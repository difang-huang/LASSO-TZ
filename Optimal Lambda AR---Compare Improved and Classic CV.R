##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Compare Improved CV and Classic CV by Density and Mutual Proportion             ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





setwd("D:/Data/National Stock Market Indeces")
LogReturns = read.csv(file = "Cleaned_Log_Returns_LongPeriod.csv", header = FALSE)
LogReturns = as.matrix(LogReturns)


BeginT = 34946   # 1995/9/4
EndT   = 44890   # 2022/11/25

Alldays = BeginT:EndT
Alldays_m = c( Alldays, rep(0, 7-length(Alldays)%%7) )
Weekdays_m = matrix(Alldays_m, nrow = 7, byrow = FALSE)
Num_Weekdays = 5*(ncol(Weekdays_m)-1) + length(Alldays)%%7
Weekdays_v = as.vector(Weekdays_m[1:5, ])[1:Num_Weekdays]

# Country Names
Countries = c("Australia","Malaysia","Indonesia","Korea","Japan","Singapore",
              "NewZealand","Philippines","Thailand","China","HK",
              
              "Netherlands","Greece","Belgium","France","Germany","Finland","Spain",
              "Ireland","Italy","Denmark","Norway","Sweden","Portugal","Russia",
              "Switzerland","UK","Poland","Turkey","Austria",
              
              "Brazil","Chile","Argentina","Mexico","Canada","US")
Country.Abbr = c("AU","MY","ID","KR","JP","SG","NZ","PH","TH","CN","HK",
                 
                 "NL","GR","BE","FR","DE","FI","ES","IE","IT","DK",
                 "NO","SE","PT","RU","CH","UK","PL","TR","AT",
                 
                 "BR","CL","AR","MX","CA","US")





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Generate Adjacency and Weighted Matrices                              ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Before Subprime => During Subprime => After Subprime
# => During European Debt => After European Debt
Period_Range = c(38930,39294, 39295,39903, 39904,40147,   40148,40624, 41625,42369)
# Before Subprime = 38930~39294 = 2006/8/1~2007/7/31
# During Subprime = 39295~39903 = 2007/8/1~2009/3/31
# After Subprime  = 39904~40147 = 2009/4/1~2009/11/30
# During European Debt = 40148~40624 = 2009/12/1 ~ 2013/12/16
# After European Debt 41625~42369 = 2013/12/17~ 2015/12/31---Optimal Choice
for(phase in 1:(length(Period_Range)/2)){
  Phase_Start = which(Weekdays_v == Period_Range[2*phase-1])
  Phase_End   = which(Weekdays_v == Period_Range[2*phase])
  XX = LogReturns[Phase_Start:Phase_End, ]
  Period = nrow(XX)
  print(paste("Phase ", phase, " has ", Period, " observations.", sep = ""))
}


N = length(Countries)
Phase.Num = length(Period_Range)/2
Rep.Times = 300 # 48
# Set.M = 10
Total.Link_classic  = matrix(0, nrow = Rep.Times, ncol = Phase.Num)
Total.Link_improve  = matrix(0, nrow = Rep.Times, ncol = Phase.Num)
Mutual.Link_classic = matrix(0, nrow = Rep.Times, ncol = Phase.Num)
Mutual.Link_improve = matrix(0, nrow = Rep.Times, ncol = Phase.Num)
StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/"


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggthemes)


for(phase in 1:Phase.Num){

# Define Variables
{
Phase_Start  = which(Weekdays_v == Period_Range[2*phase-1])
Phase_End    = which(Weekdays_v == Period_Range[2*phase])
XX = LogReturns[Phase_Start:Phase_End, ]
N = ncol(XX)  #Define the dimension of the weight matrix
Period = nrow(XX)
print(paste("Phase ", phase, " has ", Period, " observations.", sep = ""))
  
# Centralized and Standardized
X = scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1)
#diag(t(X) %*% X)
  
# Divide Data into Different Areas
Asia_area     =  1:11;  Asia    = X[, Asia_area];    N_Asia    = ncol(Asia)
Europe_area   = 12:30;  Europe  = X[, Europe_area];  N_Europe  = ncol(Europe)
America_area  = 31:36;  America = X[, America_area]; N_America = ncol(America)
  
Num_kfolds = floor(Period/30)
auto_term_loc = rep(1, N)

Num_Cores = detectCores(logical = FALSE) #- 1
chunk.size = Rep.Times/Num_Cores
}
Total.Link_classic_Period  = matrix(0, nrow = Rep.Times, ncol = N)
Total.Link_improve_Period  = matrix(0, nrow = Rep.Times, ncol = N)
Mutual.Link_classic_Period = matrix(0, nrow = Rep.Times, ncol = N)
Mutual.Link_improve_Period = matrix(0, nrow = Rep.Times, ncol = N)


# Repeat to Generate Adajacency and Weight Matrices
pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

# Define Asian, European and American Section
if(k >= 1                 & k <= N_Asia){
  y = Asia[2:Period, k]
  x = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+1          & k <= N_Asia+N_Europe){
  y = Europe[2:Period, k-N_Asia]
  x = cbind(Asia[2:Period,], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+N_Europe+1 & k <= N){
  y = America[2:Period, k-N_Asia-N_Europe]
  x = cbind(Asia[2:Period,], Europe[2:Period,], America[1:(Period-1),])
}


# DoParallel
cl = makeCluster(Num_Cores)
registerDoParallel(cl, cores = Num_Cores)
CV_M = foreach(core = 1:Num_Cores, .combine = "rbind", .packages = "glmnet") %dopar%
{
  Lambda_Est = matrix(0, nrow = chunk.size, ncol = N+1)
  for(rr in 1:chunk.size){
    cvfit = cv.glmnet(sqrt(Period-1)*x, sqrt(Period-1)*y, nfolds = Num_kfolds, 
                       alpha = 1, type.measure = "mse", nlambda = 300, 
                       intercept = FALSE, standardize = FALSE, thresh = 1e-7, 
                       penalty.factor = auto_term_loc, parallel = FALSE)
    Lambda_Est[rr, 1] = cvfit$lambda.min
    Lambda_Est[rr, 2:(N+1)] = as.vector( coef(cvfit, s="lambda.min") )[-1]
  }
  return(Lambda_Est)
}
stopImplicitCluster(); stopCluster(cl)


# Calculate Mutual Links by the Classic CV
Adjacency_Period_k = abs(sign( CV_M[, 2:(N+1)] ))
Total.Link_classic_Period[1, k]  = sum( Adjacency_Period_k[1, ] )
Mutual.Link_classic_Period[1, k] = sum( Adjacency_Period_k[1, ] )
for(rr in 2:Rep.Times){
  Total.Link_classic_Period[rr, k]  = sum( Adjacency_Period_k[rr, ] )
  Mutual.Link_classic_Period[rr, k] = sum( apply(Adjacency_Period_k[1:rr, ], 2, prod) )
}


# Calculate Mutual Links by the NEW CV
New.Adjacency_Period_k = matrix(0, nrow = Rep.Times, ncol = N)
for(rr in 1:Rep.Times){
  if(rr <= 4){ Set.M = rr }
  if(rr >  5){ Set.M = 5 }
  
  # Generate the Freqency of Lambdas
  Rep.Lambdas = CV_M[1:rr, 1]
  Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                            as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
  Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ],
                               decreasing = TRUE)[1:min(Set.M, ncol(Frq.Lambdas))]
  Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]
  # Generate the Adjacency Matrix
  Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5)+1e-9 
                         & Rep.Lambdas > max(Lambdas_Top5)-1e-9)[1]
  New.Adjacency_Period_k[rr, ] = abs( sign( CV_M[Adjacency_loc, 2:(N+1)] ) )
}
Mutual.Link_improve_Period[1, k] = sum( New.Adjacency_Period_k[1, ] )
Total.Link_improve_Period[1, k]  = sum( New.Adjacency_Period_k[1, ] )
for(rr in 2:Rep.Times){
  Mutual.Link_improve_Period[rr, k] = sum( apply(New.Adjacency_Period_k[1:rr, ], 2, prod) )
  Total.Link_improve_Period[rr, k]  = sum( New.Adjacency_Period_k[rr, ] )
}


Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)


# Calculate Mutual and Total Links in Each Period
Mutual.Link_classic[, phase] = apply(Mutual.Link_classic_Period, 1, sum)
Total.Link_classic[, phase]  = apply(Total.Link_classic_Period, 1, sum)
Mutual.Link_improve[, phase] = apply(Mutual.Link_improve_Period, 1, sum)
Total.Link_improve[, phase]  = apply(Total.Link_improve_Period, 1, sum)


print(paste("Phase ", phase, " has completed.", sep = ""))
  
}

View( Total.Link_classic/(N*(N-1)) )
View( Total.Link_improve/(N*(N-1)) )
View( Mutual.Link_classic/(N*(N-1)) )
View( Mutual.Link_improve/(N*(N-1)) )


write.csv(
  cbind(Total.Link_improve/(N*(N-1)), Total.Link_classic/(N*(N-1))),
  paste(StoragePath, "NewFigure-CV Comparison_Density.csv", sep = ""), quote = FALSE
)
write.csv(
  cbind(Mutual.Link_improve/(N*(N-1)), Mutual.Link_classic/(N*(N-1))),
  paste(StoragePath, "NewFigure-CV Comparison_Mutual Proportion.csv", sep = ""), quote = FALSE
)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Figure: Compare Improved CV and Classic CV                            ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Output Phase Names
Phs.Func = function(phs, Rep.Num = N){
  if(phs == 1){ return(rep("Before Subprime", Rep.Num)) }
  if(phs == 2){ return(rep("During Subprime", Rep.Num)) }
  if(phs == 3){ return(rep("After Subprime",  Rep.Num)) }
  if(phs == 4){ return(rep("During Eu Debt", Rep.Num)) }
  if(phs == 5){ return(rep("After Eu Debt",  Rep.Num)) }
}


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Import data
NewFig_CV.Density = read.csv(
  file = "C:/Users/581/Documents/Optimal LASSO Networks/NewFigure-CV Comparison_Density.csv",
  header = FALSE
)
NewFig_CV.Density = as.matrix(NewFig_CV.Density)

NewFig_CV.MutualP = read.csv(
  file = "C:/Users/581/Documents/Optimal LASSO Networks/NewFigure-CV Comparison_Mutual Proportion.csv",
  header = FALSE
)
NewFig_CV.MutualP = as.matrix(NewFig_CV.MutualP)


DF_CV.Density.MutualP = data.frame(
  Measure.Type = "CV", CV.Type = "improve classic", Value = 1, RepT = 1, Phase.Now = "BDA"
)
for(phase in 1:Phase.Num){
  new.df_CV.Density = data.frame(
    Measure.Type = rep("Density", Rep.Times*2),
    CV.Type = c(rep("Improved CV", Rep.Times), rep("Classic CV", Rep.Times)),
    Value = as.vector(NewFig_CV.Density[, c(phase, Phase.Num+phase)]),
    RepT = c(1:Rep.Times, 1:Rep.Times), Phase.Now = Phs.Func(phase, Rep.Num = Rep.Times*2)
  )
  new.df_CV.MutualP = data.frame(
    Measure.Type = rep("Mutual Proportion", Rep.Times*2), 
    CV.Type = c(rep("Improved CV", Rep.Times), rep("Classic CV", Rep.Times)),
    Value = as.vector(NewFig_CV.MutualP[, c(phase, Phase.Num+phase)]),
    RepT = c(1:Rep.Times, 1:Rep.Times), Phase.Now = Phs.Func(phase, Rep.Num = Rep.Times*2)
  )
  DF_CV.Density.MutualP = rbind(DF_CV.Density.MutualP, new.df_CV.Density, new.df_CV.MutualP)
}
DF_CV.Density.MutualP = DF_CV.Density.MutualP[-1, ]
# View(DF_CV.Density.MutualP)


DF_CV.Density.MutualP$CV.Type = factor(
  DF_CV.Density.MutualP$CV.Type, levels = c("Improved CV", "Classic CV")
)
DF_CV.Density.MutualP$Measure.Type = factor(
  DF_CV.Density.MutualP$Measure.Type, levels = c("Density", "Mutual Proportion")
)
DF_CV.Density.MutualP$Phase.Now = factor(
  DF_CV.Density.MutualP$Phase.Now,
  levels = c("Before Subprime", "During Subprime", "After Subprime", "During Eu Debt", "After Eu Debt")
)
BFig_CV.Density.MutualP = ggplot(
  data = DF_CV.Density.MutualP, mapping = aes(x = RepT, y = Value, color = CV.Type)
)
Fig_CV.Density.MutualP = (
  BFig_CV.Density.MutualP + geom_line() # alpha = 0.5   + guides(color = "none")
  + scale_color_manual(values = c("darkgreen", "darkorange"), name = "CV Type")
  + labs(title = "Comparison between Improved and Classic CV",
         x = "Replicate Time", y = "Proportion")
  + facet_grid(rows = vars(Phase.Now), cols = vars(Measure.Type), scales="free")
  # + ylim(0.2, 0.4)
  + theme_calc() + theme(panel.grid=element_blank())
)
Fig_CV.Density.MutualP


ggsave(filename = "NewFig_Campare.Improved.Classic.CV.eps", device = cairo_ps,
       dpi = 300, limitsize = FALSE,
       path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Figures on Stable Network Structue                                    ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





library(ggplot2)
library(grid)


# Define Figure Function
Fig_Stable.Nets = function(Classic.CV, New.CV, Num.Rep, Name.Phase, Name.Y.Values){
  DF = data.frame(
    Replicate.Time = rep(1:Num.Rep, 2), Y.Values = c(Classic.CV, New.CV),
    Y.Types = c(rep("Classic CV", Num.Rep), rep("Improved CV", Num.Rep))
  )
  DF$Y.Types = factor(DF$Y.Types, levels = c("Classic CV", "Improved CV"))
  BFig = ggplot(
    data = DF, mapping = aes(x = Replicate.Time, y = Y.Values, color = Y.Types)
  )
  Fig = (
    BFig + geom_line(size = 0.5) + guides(color = "none")
    + scale_color_manual(values = c("red", "blue"))
    + scale_x_continuous(breaks = seq(from=0,to=300,by=50),
                         labels = as.character(seq(from=0,to=300,by=50)))
    + labs(title = Name.Phase, x = "Replicate Time", y = Name.Y.Values)
  )
  return(Fig)
}


# Output Figure
Fig_Density.Bef = Fig_Stable.Nets(
  Total.Link_classic[,1]/(N*(N-1)), Total.Link_improve[,1]/(N*(N-1)), Rep.Times,
  Name.Phase = "Before", Name.Y.Values = "Density"
)
Fig_Density.Dur = Fig_Stable.Nets(
  Total.Link_classic[,2]/(N*(N-1)), Total.Link_improve[,2]/(N*(N-1)), Rep.Times,
  Name.Phase = "During", Name.Y.Values = "Density"
)
Fig_Density.Aft = Fig_Stable.Nets(
  Total.Link_classic[,3]/(N*(N-1)), Total.Link_improve[,3]/(N*(N-1)), Rep.Times,
  Name.Phase = "After", Name.Y.Values = "Density"
)

Fig_Mutual.Bef = Fig_Stable.Nets(
  Mutual.Link_classic[,1]/(N*(N-1)), Mutual.Link_improve[,1]/(N*(N-1)), Rep.Times,
  Name.Phase = "Before", Name.Y.Values = "Mutual Proportion"
)
Fig_Mutual.Dur = Fig_Stable.Nets(
  Mutual.Link_classic[,2]/(N*(N-1)), Mutual.Link_improve[,2]/(N*(N-1)), Rep.Times,
  Name.Phase = "During", Name.Y.Values = "Mutual Proportion"
)
Fig_Mutual.Aft = Fig_Stable.Nets(
  Mutual.Link_classic[,3]/(N*(N-1)), Mutual.Link_improve[,3]/(N*(N-1)), Rep.Times,
  Name.Phase = "After", Name.Y.Values = "Mutual Proportion"
)


grid.newpage()  ###?½?ͼ??????
pushViewport(viewport(layout = grid.layout(3,2))) ####???????ֳ?3*2????
vplayout = function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}

print(Fig_Density.Bef, vp = vplayout(1, 1))
print(Fig_Density.Dur, vp = vplayout(2, 1))
print(Fig_Density.Aft, vp = vplayout(3, 1))
print(Fig_Mutual.Bef,  vp = vplayout(1, 2))
print(Fig_Mutual.Dur,  vp = vplayout(2, 2))
print(Fig_Mutual.Aft,  vp = vplayout(3, 2))












plot(Mutual.Link_classic[,1]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.15,0.35))
lines(Mutual.Link_improve[,1]/(N*(N-1)), lty = 1, col = "red")
plot(Mutual.Link_classic[,2]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.3,0.45))
lines(Mutual.Link_improve[,2]/(N*(N-1)), lty = 1, col = "red")
plot(Mutual.Link_classic[,3]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.1,0.38))
lines(Mutual.Link_improve[,3]/(N*(N-1)), lty = 1, col = "red")


plot(Total.Link_classic[,1]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.15,0.35))
lines(Total.Link_improve[,1]/(N*(N-1)), lty = 1, col = "red")
plot(Total.Link_classic[,2]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.3,0.45))
lines(Total.Link_improve[,2]/(N*(N-1)), lty = 1, col = "red")
plot(Total.Link_classic[,3]/(N*(N-1)), type = "l", lty = 1, ylim = c(0.15,0.35))
lines(Total.Link_improve[,3]/(N*(N-1)), lty = 1, col = "red")












StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/2008 Subprime Crisis by"
StoragePath_FA = paste(StoragePath, "AR Opitmal Lambdas-FA.csv")
write.csv(Final.Adjacency,  StoragePath_FA,  quote = FALSE)
StoragePath_FW = paste(StoragePath, "AR Opitmal Lambdas-FW.csv")
write.csv(Final.Weight,  StoragePath_FW,  quote = FALSE)
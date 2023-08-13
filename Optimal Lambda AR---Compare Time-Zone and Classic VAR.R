##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Compare Time-Zone and Classic VAR Model by In-Sample R-Square                   ##
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
nlinks_Unsign = rep(0, Phase.Num)
nlinks_counter = rep(0, Phase.Num)
Rep.Times = 3000
Set.M = 5 # 10
StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/"

Type.VAR = "both"  # Type.VAR = c("time zone", "classic", "mix", "all", "both")
# both = time-zone model + classic model
# all  = time-zone model + classic model + mix model
# mix model = put time-zone terms and classic lag terms together


Final.Adj_tz  = matrix(0, nrow = Phase.Num*N, ncol = N)
Final.Wet_tz  = matrix(0, nrow = Phase.Num*N, ncol = N)
Final.Adj_old = matrix(0, nrow = Phase.Num*N, ncol = N)
Final.Wet_old = matrix(0, nrow = Phase.Num*N, ncol = N)
Final.Adj_mix = matrix(0, nrow = Phase.Num*(36+30), ncol = N)
Final.Wet_mix = matrix(0, nrow = Phase.Num*(36+30), ncol = N)
CV.Error_tz  = matrix(0, nrow = Phase.Num*3, ncol = N)
CV.Error_old = matrix(0, nrow = Phase.Num*3, ncol = N)


mse_tz  = matrix(0, nrow = N, ncol = Phase.Num)
mse_old = matrix(0, nrow = N, ncol = Phase.Num)


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
# auto_term_loc = rep(1, N)

Num_Cores = detectCores(logical = FALSE) # - 1
chunk.size = Rep.Times/Num_Cores
}


##=====================================================================================##
##=====================================================================================##  
##=====================================================================================##


# Repeat to Generate Adajacency and Weight Matrices===Time-Zone VAR
if(Type.VAR == "time zone" || Type.VAR == "all" || Type.VAR == "both"){

Adj_Period_tz = matrix(0, nrow = N, ncol = N)
Wet_Period_tz = matrix(0, nrow = N, ncol = N)
auto_term_loc = rep(1, N)

pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

# Define Asia, Europe and America Section
if(k >= 1                 & k <= N_Asia){
  y = Asia[2:Period, k]
  x_tz = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+1          & k <= N_Asia+N_Europe){
  y = Europe[2:Period, k-N_Asia]
  x_tz = cbind(Asia[2:Period,], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+N_Europe+1 & k <= N){
  y = America[2:Period, k-N_Asia-N_Europe]
  x_tz = cbind(Asia[2:Period,], Europe[2:Period,], America[1:(Period-1),])
}


# DoParallel
cl = makeCluster(Num_Cores)
registerDoParallel(cl, cores = Num_Cores)
CV_M = foreach(core = 1:Num_Cores, .combine = "rbind", .packages = "glmnet") %dopar%
{
  Lambda_Est = matrix(0, nrow = chunk.size, ncol = 1+N+2)
  for(rr in 1:chunk.size){
    cvfit = cv.glmnet(sqrt(Period-1)*x_tz, sqrt(Period-1)*y, nfolds = Num_kfolds, 
                       alpha = 1, type.measure = "mse", nlambda = 300, 
                       intercept = FALSE, standardize = FALSE, thresh = 1e-7, 
                       penalty.factor = auto_term_loc, parallel = FALSE)
    Lambda_Est[rr, 1] = cvfit$lambda.min
    Lambda_Est[rr, 2:(N+1)] = as.vector( coef(cvfit, s="lambda.min") )[-1]
    min.CV.Error_loc = which.min(cvfit$cvm)
    # cvm = mean cross-validated error
    Lambda_Est[rr, N+2] = cvfit$cvm[min.CV.Error_loc]
    # cvsd = estimate of standard error of cvm
    Lambda_Est[rr, N+3] = cvfit$cvsd[min.CV.Error_loc] 
  }
  return(Lambda_Est)
}
stopImplicitCluster(); stopCluster(cl)

# Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]

# Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5, na.rm = TRUE)+1e-9 
                       & Rep.Lambdas > max(Lambdas_Top5, na.rm = TRUE)-1e-9)[1]
Adj_Period_tz[, k] = sign( CV_M[Adjacency_loc, 2:(N+1)] )
# upper curve = cvm + cvsd
CV.Error_tz[phase*3-2, k] = CV_M[Adjacency_loc, N+2] + CV_M[Adjacency_loc, N+3]
CV.Error_tz[phase*3-1, k] = CV_M[Adjacency_loc, N+2]
# lower curve = cvm - cvsd
CV.Error_tz[phase*3,   k] = CV_M[Adjacency_loc, N+2] - CV_M[Adjacency_loc, N+3]


# mse of time-zone VAR model
mse_tz[k, phase] = sum( (y - x_tz %*% CV_M[Adjacency_loc, 2:(N+1)])^2 )


# Generate the Weight Matrix
Betas_Top5_w = Frq.Lambdas[2, Lambdas_Top5_loc]/sum(Frq.Lambdas[2, Lambdas_Top5_loc])
Betas_Top5_coefs = rep(0, N)
for(rr in 1:Set.M){
  Betas_loc = which(Rep.Lambdas < Lambdas_Top5[rr]+1e-9 
                     & Rep.Lambdas > Lambdas_Top5[rr]-1e-9)[1]
  Betas_Top5_coefs = Betas_Top5_coefs + Betas_Top5_w[rr] * CV_M[Betas_loc, 2:(N+1)]
}
Wet_Period_tz[, k] = Betas_Top5_coefs * abs(Adj_Period_tz[, k])


Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)

Final.Adj_tz[(phase*N-N+1):(phase*N), ] = Adj_Period_tz
Final.Wet_tz[(phase*N-N+1):(phase*N), ] = Wet_Period_tz

print(paste("Time-Zone VAR: Phase ", phase, " has completed.", sep = ""))

}


##=====================================================================================##
##=====================================================================================##  
##=====================================================================================##


# Repeat to Generate Adajacency and Weight Matrices===Classic VAR
if(Type.VAR == "classic" || Type.VAR == "all" || Type.VAR == "both"){

Adj_Period_old = matrix(0, nrow = N, ncol = N)
Wet_Period_old = matrix(0, nrow = N, ncol = N)
auto_term_loc = rep(1, N)

pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){
    
# Define Asian, Europe and America Section
if(k >= 1                 & k <= N_Asia){
  y = Asia[2:Period, k]
  x_old = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+1          & k <= N_Asia+N_Europe){
  y = Europe[2:Period, k-N_Asia]
  x_old = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+N_Europe+1 & k <= N){
  y = America[2:Period, k-N_Asia-N_Europe]
  x_old = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}


# DoParallel
cl = makeCluster(Num_Cores)
registerDoParallel(cl, cores = Num_Cores)
CV_M = foreach(core = 1:Num_Cores, .combine = "rbind", .packages = "glmnet") %dopar%
{
  Lambda_Est = matrix(0, nrow = chunk.size, ncol = 1+N+2)
  for(rr in 1:chunk.size){
    cvfit = cv.glmnet(sqrt(Period-1)*x_old, sqrt(Period-1)*y, nfolds = Num_kfolds, 
                       alpha = 1, type.measure = "mse", nlambda = 300, 
                       intercept = FALSE, standardize = FALSE, thresh = 1e-7, 
                       penalty.factor = auto_term_loc, parallel = FALSE)
    Lambda_Est[rr, 1] = cvfit$lambda.min
    Lambda_Est[rr, 2:(N+1)] = as.vector( coef(cvfit, s="lambda.min") )[-1]
    min.CV.Error_loc = which.min(cvfit$cvm)
    Lambda_Est[rr, N+2] = cvfit$cvm[min.CV.Error_loc]
    Lambda_Est[rr, N+3] = cvfit$cvsd[min.CV.Error_loc]
  }
  return(Lambda_Est)
}
stopImplicitCluster(); stopCluster(cl)


# Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]


# Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5, na.rm = TRUE)+1e-9 
                       & Rep.Lambdas > max(Lambdas_Top5, na.rm = TRUE)-1e-9)[1]
Adj_Period_old[, k] = sign( CV_M[Adjacency_loc, 2:(N+1)] )
CV.Error_old[phase*3-2, k] = CV_M[Adjacency_loc, N+2] + CV_M[Adjacency_loc, N+3]
CV.Error_old[phase*3-1, k] = CV_M[Adjacency_loc, N+2]
CV.Error_old[phase*3,   k] = CV_M[Adjacency_loc, N+2] - CV_M[Adjacency_loc, N+3]


# mse of classic VAR model
mse_old[k, phase] = sum( (y - x_old %*% CV_M[Adjacency_loc, 2:(N+1)])^2 )


# Generate the Weight Matrix
Betas_Top5_w = Frq.Lambdas[2, Lambdas_Top5_loc]/sum(Frq.Lambdas[2, Lambdas_Top5_loc])
Betas_Top5_coefs = rep(0, N)
for(rr in 1:Set.M){
  Betas_loc = which(Rep.Lambdas < Lambdas_Top5[rr]+1e-9 
                     & Rep.Lambdas > Lambdas_Top5[rr]-1e-9)[1]
  Betas_Top5_coefs = Betas_Top5_coefs + Betas_Top5_w[rr] * CV_M[Betas_loc, 2:(N+1)]
}
Wet_Period_old[, k] = Betas_Top5_coefs * abs(Adj_Period_old[, k])


Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)

Final.Adj_old[(phase*N-N+1):(phase*N), ] = Adj_Period_old
Final.Wet_old[(phase*N-N+1):(phase*N), ] = Wet_Period_old

print(paste("Classic VAR: Phase ", phase, " has completed.", sep = ""))

}


##=====================================================================================##
##=====================================================================================##  
##=====================================================================================##


# Repeat to Generate Adajacency and Weight Matrices===Mix VAR
if(Type.VAR == "mix" || Type.VAR == "all"){

Adj_Period_mix = matrix(0, nrow = N+30, ncol = N)
Wet_Period_mix = matrix(0, nrow = N+30, ncol = N)  
  
pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){
    
# Define Asian, Europe and America Section
if(k >= 1                 & k <= N_Asia){
  y = Asia[2:Period, k]
  x_mix = cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
}
if(k >= N_Asia+1          & k <= N_Asia+N_Europe){
  y = Europe[2:Period, k-N_Asia]
  x_mix = cbind(Asia[2:Period,], Europe[1:(Period-1),], America[1:(Period-1),],
                 Asia[1:(Period-1),])
  x_mix = matrix(as.numeric(x_mix), nrow = Period-1, byrow = FALSE)
}
if(k >= N_Asia+N_Europe+1 & k <= N){
  y = America[2:Period, k-N_Asia-N_Europe]
  x_mix = cbind(Asia[2:Period,], Europe[2:Period,], America[1:(Period-1),],
                 Asia[1:(Period-1),], Europe[1:(Period-1),])
}
auto_term_loc = rep(1, ncol(x_mix))


# DoParallel
cl = makeCluster(Num_Cores)
registerDoParallel(cl, cores = Num_Cores)
CV_M = foreach(core = 1:Num_Cores, .combine = "rbind", .packages = "glmnet") %dopar%
{
  Lambda_Est = matrix(0, nrow = chunk.size, ncol = 1+ncol(x_mix))
  for(rr in 1:chunk.size){
    cvfit = cv.glmnet(sqrt(Period-1)*x_mix, sqrt(Period-1)*y, nfolds = Num_kfolds, 
                       alpha = 1, type.measure = "mse", nlambda = 300, 
                       intercept = FALSE, standardize = FALSE, thresh = 1e-7, 
                       penalty.factor = auto_term_loc, parallel = FALSE)
    Lambda_Est[rr, 1] = cvfit$lambda.min
    Lambda_Est[rr, 2:(1+ncol(x_mix))] = as.vector( coef(cvfit, s="lambda.min") )[-1]
  }
  return(Lambda_Est)
}
stopImplicitCluster(); stopCluster(cl)


# Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]


# Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5, na.rm = TRUE)+1e-9 
                       & Rep.Lambdas > max(Lambdas_Top5)-1e-9, na.rm = TRUE)[1]
Adj_Period_mix[1:ncol(x_mix), k] = sign( CV_M[Adjacency_loc, 2:(1+ncol(x_mix))] )


# Generate the Weight Matrix
Betas_Top5_w = Frq.Lambdas[2, Lambdas_Top5_loc]/sum(Frq.Lambdas[2, Lambdas_Top5_loc])
Betas_Top5_coefs = rep(0, ncol(x_mix))
for(rr in 1:Set.M){
  Betas_loc = which(Rep.Lambdas < Lambdas_Top5[rr]+1e-9 
                     & Rep.Lambdas > Lambdas_Top5[rr]-1e-9)[1]
  Betas_Top5_coefs = Betas_Top5_coefs+Betas_Top5_w[rr]*CV_M[Betas_loc, 2:(1+ncol(x_mix))]
}
Wet_Period_mix[1:ncol(x_mix), k] = Betas_Top5_coefs*abs(Adj_Period_mix[1:ncol(x_mix), k])


Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)

Final.Adj_mix[(phase*66-66+1):(phase*66), ] = Adj_Period_mix
Final.Wet_mix[(phase*66-66+1):(phase*66), ] = Wet_Period_mix

print(paste("Mix VAR: Phase ", phase, " has completed.", sep = ""))

}

}


NewFig_R.Square.IS = 1 - mse_tz/mse_old
View(NewFig_R.Square.IS)

write.csv(NewFig_R.Square.IS,
          paste(StoragePath, "NewTable-R2 In-Sample Performance.csv", sep = ""),
          quote = FALSE)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Figure: In-Sample R-Sequare Performance                                    ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Output Phase Names
Phs.Func = function(phs, Rep.Num = N){
  if(phs == 1){ return(rep("Before 2008 Subprime Crisis", Rep.Num)) }
  if(phs == 2){ return(rep("During 2008 Subprime Crisis", Rep.Num)) }
  if(phs == 3){ return(rep("After 2008 Subprime Crisis",  Rep.Num)) }
  if(phs == 4){ return(rep("During European Debt Crisis", Rep.Num)) }
  if(phs == 5){ return(rep("After European Debt Crisis",  Rep.Num)) }
}


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Generate ggplot-type data
NewFig_R.Square.IS = read.csv(
  file = "C:/Users/581/Documents/Optimal LASSO Networks/NewTable-R2 In-Sample Performance.csv",
  header = FALSE
)
NewFig_R.Square.IS = as.matrix(NewFig_R.Square.IS)

DF_IS.R2 = data.frame(Market = "US", R2.Value = 1, Continent = "AsEuAm", Phase.Now = "BDA")
for(phase in 1:Phase.Num){
  new.df_R2 = data.frame(
    Market = Country.Abbr, R2.Value = as.vector( as.matrix(NewFig_R.Square.IS) ), 
    Continent = c(rep("Asia",11), rep("Europe",19), rep("America",6)),
    Phase.Now = switch(phase, "Before Subprime", "During Subprime", "After Subprime",
                       "During Eu Debt", "After Eu Debt")
  )
  DF_IS.R2 = rbind(DF_IS.R2, new.df_R2)
}
DF_IS.R2 = DF_IS.R2[-1, ]
# View(DF_IS.R2)


RBG = c("red", "blue", "green")
RBG = setNames(RBG, c("Asia", "Europe", "America"))
DF_IS.R2$Market = factor(DF_IS.R2$Market, levels = Country.Abbr)
DF_IS.R2$Continent = factor(DF_IS.R2$Continent, levels = c("Asia", "Europe", "America"))
DF_IS.R2$Phase.Now = factor(DF_IS.R2$Phase.Now,
  levels = c("Before Subprime", "During Subprime", "After Subprime", "During Eu Debt",
             "After Eu Debt")
)
BFig_IS.R2 = ggplot(
  data = DF_IS.R2,
  mapping = aes(x = Market, y = R2.Value, fill = Continent, shape = Continent)
)
Fig_IS.R2 = (
  BFig_IS.R2 + geom_point(size = 2.5, alpha = 0.5)
  + scale_color_manual(values = RBG) + guides(color = "none")
  + scale_shape_manual(values = c(16,18,17))
  + labs(title = "In-Sample Performance Evaluated by R-Square", y = "R-Square")
  + facet_grid(rows = vars(Phase.Now), cols = vars(Continent), scales="free", space = "free")
  + ylim(-0.05, 0.8)
  + theme_calc() + theme(panel.grid=element_blank())
)
Fig_IS.R2


ggsave(filename = "NewFig_In-Sample.Performance.R2.eps", device = cairo_ps,
       dpi = 300, limitsize = FALSE,
       path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)







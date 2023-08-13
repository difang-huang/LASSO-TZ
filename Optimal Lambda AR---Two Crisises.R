##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Analysis on Two Crisises by Optimal Lambda AR                                   ##
##     Two Crisises: 2008 Subprime + Europe Debt                                       ##
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
{
Countries = c("Australia","Malaysia","Indonesia","Korea","Japan","Singapore",
              "NewZealand","Philippines","Thailand","China","HK",
             
              "Netherlands","Greece","Belgium","France","Germany","Finland","Spain",
              "Ireland","Italy","Denmark","Norway","Sweden","Portugal","Russia",
              "Switzerland","UK","Poland","Turkey","Austria",
             
              "Brazil","Chile","Argentina","Mexico","Canada","US")
}



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
# After  European Debt = 41625~42369 = 2013/12/17~ 2015/12/31---Optimal Choice
for(phase in 1:(length(Period_Range)/2)){
  Phase_Start = which(Weekdays_v == Period_Range[2*phase-1])
  Phase_End   = which(Weekdays_v == Period_Range[2*phase])
  XX = LogReturns[Phase_Start:Phase_End, ]
  Period = nrow(XX)
  print(paste("Phase ", phase, " has ", Period, " observations.", sep = ""))
}


Phase_Num = length(Period_Range)/2
nlinks_Unsign = rep(0, Phase_Num)
nlinks_counter = rep(0, Phase_Num)
AR_Terms = matrix(0, nrow = length(Countries), ncol = Phase_Num+1)
AR_Terms[, 1] = Countries
Rep.Times = 10000  #3000
Set.M = 10 # 5
Final.Adjacency = matrix(0, nrow = Phase_Num*length(Countries), ncol = length(Countries))
Final.Weight    = matrix(0, nrow = Phase_Num*length(Countries), ncol = length(Countries))
StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/"


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)


for(phase in 1:Phase_Num){

# Define Variables
{
Phase_Start = which(Weekdays_v == Period_Range[2*phase-1])
Phase_End   = which(Weekdays_v == Period_Range[2*phase])
XX = LogReturns[Phase_Start:Phase_End, ]
N = ncol(XX)  #Define the dimension of the weight matrix
Period = nrow(XX)
MaxLinks = N*(N-1)
print(paste("Phase ", phase, " has ", Period, " observations.", sep = ""))
  
# Centralized and Standardized
X = scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1)
#diag(t(X) %*% X)
  
# Divide Data into Different Areas
Asia_area    =  1:11;  Asia    = X[, Asia_area];    N_Asia    = ncol(Asia)
Europe_area  = 12:30;  Europe  = X[, Europe_area];  N_Europe  = ncol(Europe)
America_area = 31:36;  America = X[, America_area]; N_America = ncol(America)

Num_kfolds = max(floor(Period/30), 3) # floor(Period/30)
auto_term_loc =  rep(1, N)
Adjacency_Period = matrix(0, nrow = N, ncol = N)
Weight_Period    = matrix(0, nrow = N, ncol = N)

Num_Cores = detectCores(logical = FALSE) # - 1
chunk.size = Rep.Times/Num_Cores
}


# Repeat to Generate Adajacency and Weight Matrices
pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

# Define Asian, Europe and America Section
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

# Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]

# Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5)+1e-9 
                       & Rep.Lambdas > max(Lambdas_Top5)-1e-9)[1]
Adjacency_Period[, k] = sign( CV_M[Adjacency_loc, 2:(N+1)] )

# Generate the Weight Matrix
Betas_Top5_w = Frq.Lambdas[2, Lambdas_Top5_loc]/sum(Frq.Lambdas[2, Lambdas_Top5_loc])
Betas_Top5_coefs = rep(0, N)
for(rr in 1:Set.M){
  Betas_loc = which(Rep.Lambdas < Lambdas_Top5[rr]+1e-9 
                     & Rep.Lambdas > Lambdas_Top5[rr]-1e-9)[1]
  Betas_Top5_coefs = Betas_Top5_coefs + Betas_Top5_w[rr] * CV_M[Betas_loc, 2:(N+1)]
}
Weight_Period[, k] = Betas_Top5_coefs * abs(Adjacency_Period[, k])

Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)


Final.Adjacency[(phase*N-N+1):(phase*N), ] = Adjacency_Period
Final.Weight[(phase*N-N+1):(phase*N), ]    = Weight_Period


# Analysis on Autogressive Terms and Link Numbers
AR_Terms[, phase+1] = diag(Weight_Period)
# diag(Adjacency_Period) = rep(0, N)
nlinks_Unsign[phase] = sum( abs(Adjacency_Period) )
nlinks_counter[phase] = sum( Adjacency_Period )


print(paste("The Phase ", phase, " has completed.", sep = ""))

}
nlinks_Unsign; nlinks_counter
View(AR_Terms)
# StoragePath_AR.Terms = paste(StoragePath, "Two Sequent Crisis-AR Terms.csv", sep = "")
# write.csv(AR_Terms, StoragePath_AR.Terms,  quote = FALSE)


StoragePath_FA = paste(StoragePath, "Two Sequent Crisis-FA_1e-4 Replications.csv", sep = "")
StoragePath_FW = paste(StoragePath, "Two Sequent Crisis-FW_1e-4 Replications.csv", sep = "")
write.csv(Final.Adjacency, StoragePath_FA,  quote = FALSE)
write.csv(Final.Weight,    StoragePath_FW,  quote = FALSE)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Output Nodes and Links Lists for Gephi                                ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##



# Define the Continent Detection Function
AsEuAm.Detect = function(k){
  if(k <= N_Asia){ return("Asia") }
  if(k >= N_Asia+1 && k <= N_Asia+N_Europe){ return("Europe") }
  if(k >= N_Europe+1){ return("America") }
}


# Output Links Lists for Gephi
for(phase in 1:Phase_Num){
Adjacency = Final.Adjacency[(phase*N-N+1):(phase*N), ]
Weights   = Final.Weight[(phase*N-N+1):(phase*N), ]
diag(Adjacency) = rep(0, N); diag(Weights) = rep(0, N)

Links_T = c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
Links_P = c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
Links_N = c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
for(i in 1:N){
for(j in 1:N){
  judgement = Adjacency[i, j]
  if(judgement == 1){
    New.Link = c(Countries[i], Countries[j],"Directed", "P", 
                  abs(Weights[i, j]), AsEuAm.Detect(i))
    Links_T = rbind(Links_T, New.Link)
    Links_P = rbind(Links_P, New.Link)
  }
  if(judgement == -1){
    New.Link = c(Countries[i], Countries[j],"Directed", "N", 
                  abs(Weights[i, j]), AsEuAm.Detect(i))
    Links_T = rbind(Links_T, New.Link)
    Links_N = rbind(Links_N, New.Link)
  }
}
}

StPath_Gephi_T = paste(StoragePath,
                       "New Data_Lambda-Full Graph Phase ", phase, ".csv", sep = "")
write.csv(Links_T, StPath_Gephi_T, quote = FALSE)

StPath_Gephi_P = paste(StoragePath,
                       "New Data_Lambda-Sub Graph Phase ", phase, "-P", ".csv", sep = "")
write.csv(Links_P, StPath_Gephi_P, quote = FALSE)

StPath_Gephi_N = paste(StoragePath,
                       "New Data_Lambda-Sub Graph Phase ", phase, "-N", ".csv", sep = "")
write.csv(Links_N, StPath_Gephi_N, quote = FALSE)
}










##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##



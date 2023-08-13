##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Rolling Window for Time-Zone VAR Model with LASSO                               ##
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
##               Continental Adjacency Weight Function                                 ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Define the Continental Adjacency Function
Cont.Summary = function(adjacency, weight, type = c("U", "P", "N"),
                        output.adjacency = TRUE, Asia_area = 1:11,
                        Europe_area = 12:30, America_area = 31:36){

Cont.Adjacency = matrix(0, nrow = 3, ncol = 3)
Cont.Weight    = matrix(0, nrow = 3, ncol = 3)
if(type == "U"){ ct.adjacency = abs(adjacency) }               # Unsign
if(type == "P"){ ct.adjacency = (abs(adjacency)+adjacency)/2 } # Positive
if(type == "N"){ ct.adjacency = (abs(adjacency)-adjacency)/2 } # Negative
ct.weight = ct.adjacency * abs(weight)


if(output.adjacency == TRUE){

Cont.Adjacency[1, 1] = sum( ct.adjacency[Asia_area, Asia_area] )
Cont.Adjacency[1, 2] = sum( ct.adjacency[Asia_area, Europe_area] )
Cont.Adjacency[1, 3] = sum( ct.adjacency[Asia_area, America_area] )
Cont.Adjacency[2, 1] = sum( ct.adjacency[Europe_area, Asia_area] )
Cont.Adjacency[2, 2] = sum( ct.adjacency[Europe_area, Europe_area] )
Cont.Adjacency[2, 3] = sum( ct.adjacency[Europe_area, America_area] )
Cont.Adjacency[3, 1] = sum( ct.adjacency[America_area, Asia_area] )
Cont.Adjacency[3, 2] = sum( ct.adjacency[America_area, Europe_area] )
Cont.Adjacency[3, 3] = sum( ct.adjacency[America_area, America_area] )

return(Cont.Adjacency)
}


if(output.adjacency == FALSE){

Cont.Weight[1, 1] = sum( ct.weight[Asia_area, Asia_area] )
Cont.Weight[1, 2] = sum( ct.weight[Asia_area, Europe_area] )
Cont.Weight[1, 3] = sum( ct.weight[Asia_area, America_area] )
Cont.Weight[2, 1] = sum( ct.weight[Europe_area, Asia_area] )
Cont.Weight[2, 2] = sum( ct.weight[Europe_area, Europe_area] )
Cont.Weight[2, 3] = sum( ct.weight[Europe_area, America_area] )
Cont.Weight[3, 1] = sum( ct.weight[America_area, Asia_area] )
Cont.Weight[3, 2] = sum( ct.weight[America_area, Europe_area] )
Cont.Weight[3, 3] = sum( ct.weight[America_area, America_area] )

return(Cont.Weight)
}

}





##=====================================================================================##
##=====================================================================================##
##               Generate Adjacency and Weighted Matrices                              ##
##=====================================================================================##
##=====================================================================================##





# Data Input
Period_Range = c(36899, 44890)   # From = 36899 = 2001/1/8, To = 44890 = 2022/11/25
Period_Start = which(Weekdays_v == Period_Range[1])
Period_End   = which(Weekdays_v == Period_Range[2])
Roll.Win = 150 #100
New.h = 5 # 40 10
Num_kfolds = 10 # Roll.Win/New.h
N = length(Countries)

# Phase.Num = (Period_End - Period_Start - Roll.Win + 1)/New.h + 1
Phase.Num = floor( (Period_End - Period_Start - Roll.Win + 1)/New.h ) + 1
StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/"

Density = matrix(0, nrow = Phase.Num, ncol = 3)
Continent.Weight_U = matrix(0, nrow = Phase.Num, ncol = 3*3)
Continent.Weight_P = matrix(0, nrow = Phase.Num, ncol = 3*3)
Continent.Weight_N = matrix(0, nrow = Phase.Num, ncol = 3*3)
Final.Adjacency_rollwin = matrix(0, nrow = Phase.Num*N, ncol = N)
Final.Weight_rollwin    = matrix(0, nrow = Phase.Num*N, ncol = N)


Rep.Times = 1600 # 960 20
PreSet.M = 5


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggthemes)
library(grid)


time.start = Sys.time()
for(phase in 1:Phase.Num){ # 44:46


# Define Variables
XX = LogReturns[(Period_Start+(phase-1)*New.h):(Period_Start-1+Roll.Win+(phase-1)*New.h), ]
Period = nrow(XX); Max.Link = N*(N-1)
X = scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1) # Centralized and Standardized
auto_term_loc = rep(1, N)
Adjacency_Period = matrix(0, nrow = N, ncol = N)
Weight_Period    = matrix(0, nrow = N, ncol = N)

Num_Cores = detectCores(logical = TRUE) # - 1
chunk.size = Rep.Times/Num_Cores


# Divide Data into Different Areas
Asia_area     =  1:11;  Asia    = X[, Asia_area];    N_Asia    = ncol(Asia)
Europe_area   = 12:30;  Europe  = X[, Europe_area];  N_Europe  = ncol(Europe)
America_area  = 31:36;  America = X[, America_area]; N_America = ncol(America)


# Repeat to Generate Adajacency and Weight Matrices
pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

## Define Asian, Europe and America Section
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

## DoParallel
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

## Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
if(ncol(Frq.Lambdas) <= PreSet.M) {Set.M = ncol(Frq.Lambdas)}
if(ncol(Frq.Lambdas) >  PreSet.M) {Set.M = PreSet.M}
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]

## Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5)+1e-9 
                      & Rep.Lambdas > max(Lambdas_Top5)-1e-9)[1]
Adjacency_Period[, k] = sign( CV_M[Adjacency_loc, 2:(N+1)] )

## Generate the Weight Matrix
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


Density[phase, 1] = sum( abs(Adjacency_Period) )/Max.Link
Density[phase, 2] = sum( (abs(Adjacency_Period)+Adjacency_Period)/2 )/Max.Link
Density[phase, 3] = sum( (abs(Adjacency_Period)-Adjacency_Period)/2 )/Max.Link


Continent.Weight_U[phase, ] = as.vector( t(Cont.Summary(
  adjacency = Adjacency_Period, weight = Weight_Period, type = "U", output.adjacency = FALSE
)) )
Continent.Weight_P[phase, ] = as.vector( t(Cont.Summary(
  adjacency = Adjacency_Period, weight = Weight_Period, type = "P", output.adjacency = FALSE
)) )
Continent.Weight_N[phase, ] = as.vector( t(Cont.Summary(
  adjacency = Adjacency_Period, weight = Weight_Period, type = "N", output.adjacency = FALSE
)) )


Final.Adjacency_rollwin[(phase*N-N+1):(phase*N), 1:N] = Adjacency_Period
Final.Weight_rollwin[(phase*N-N+1):(phase*N), 1:N]    = Weight_Period


print(paste("Phase ", phase, " (", 100*round(phase/Phase.Num, digit =3), 
            "%) has completed (", Phase.Num, " phases involved).", sep = ""))


}
time.end = Sys.time()
print(time.end - time.start)


View(Density)
View(Continent.Weight_U)
View(Continent.Weight_P)
View(Continent.Weight_N)



StoragePath = "D:/Research/2019  Financial Risk Networks by LASSO/New Data/"
write.csv(Density, paste(StoragePath, "NewFigure-Roll.Window_Desity.csv", sep = ""), quote = FALSE)

write.csv(Continent.Weight_U,
          paste(StoragePath, "NewFigure-Roll.Window_Continent.Weight_Unsigned.csv", sep = ""),
          quote = FALSE)
write.csv(Continent.Weight_P,
          paste(StoragePath, "NewFigure-Roll.Window_Continent.Weight_Positive.csv", sep = ""),
          quote = FALSE)
write.csv(Continent.Weight_N,
          paste(StoragePath, "NewFigure-Roll.Window_Continent.Weight_Negative.csv", sep = ""),
          quote = FALSE)



##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Figure: Dynamic Density                                               ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##



# Figure for Density
Density = read.csv(
  file = "C:/Users/581/Documents/Optimal LASSO Networks/Desity_Roll.Window---60+10.csv",
  header = TRUE)
Density = as.matrix(Density)[, -1]


DF_Density = data.frame(
  Time = rep(1:Phase.Num, ncol(Density)),
  Value = c(Density[, 1], Density[, 2], Density[, 3]),
  Type = c(rep("Unsigned", Phase.Num), rep("Positive", Phase.Num), 
           rep("Negative", Phase.Num))
)
DF_Density$Type = factor(DF_Density$Type,
                         levels = c("Unsigned", "Positive", "Negative"))
BFig_Density = ggplot(
  data = DF_Density, mapping = aes(x = Time, y = Value, color = Type)
)
Fig_Density = (
  BFig_Density + geom_line(linetype = "solid", size = 1) #+ guides(color = "none")
  + scale_color_manual(values = c("black", "red", "blue"))
  # + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
  #                      labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Density")
  + theme_calc() + theme(panel.grid=element_blank())
)
Fig_Density





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Figure for Continent Weights (3*3 Format)                             ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Function to Output Country Names
Cont.Name_Fig = function(ii, jj){
  source.country = switch(ii, "As", "Eu", "Am")
  target.country = switch(jj, "As", "Eu", "Am")
  return(paste("From ", source.country, " to ", target.country, sep = ""))
}


# Function for Multi Graphs
vplayout = function(x,y){ viewport(layout.pos.row = x, layout.pos.col = y) }


# Input Data
# setwd("C:/Users/581/Documents/Optimal LASSO Networks")
# Continent.Weight_U = read.csv(file = "Continent.Weight.Unsigned_Roll.Window.csv",
#                               header = TRUE)
# Continent.Weight_U = as.matrix(Continent.Weight_U)[, -1]
# 
# Continent.Weight_P = read.csv(file = "Continent.Weight.Positive_Roll.Window.csv",
#                               header = TRUE)
# Continent.Weight_P = as.matrix(Continent.Weight_P)[, -1]
# 
# Continent.Weight_N = read.csv(file = "Continent.Weight.Negative_Roll.Window.csv",
#                               header = TRUE)
# Continent.Weight_N = as.matrix(Continent.Weight_N)[, -1]


# Figure for Continent Weights (3*3 Format)
grid.newpage()
pushViewport( viewport(layout = grid.layout(3, 3)) )
for(ii in 1:3){ for(jj in 1:3){

kk = 3*(ii - 1) + jj

DF_Cont.Strength = data.frame(
  Time = rep(1:Phase.Num, 3),
  Value = c(Continent.Weight_U[, kk], Continent.Weight_P[, kk], Continent.Weight_N[, kk]),
  Link.Type = c(rep("Unsigned", Phase.Num), rep("Positive", Phase.Num),
                rep("Negative", Phase.Num))
)
DF_Cont.Strength$Link.Type = factor(DF_Cont.Strength$Link.Type,
                                    levels = c("Unsigned", "Positive", "Negative"))
BFig_Cont.Strength = ggplot(
  data = DF_Cont.Strength,
  mapping = aes(x = Time, y = Value, color = Link.Type)
)
Fig_Cont.Strength = (
  BFig_Cont.Strength + geom_line(linetype = "solid", size = 1)
  + guides(color = "none") #+ ylim(0, 20)
  + scale_color_manual(values = c("black", "red", "blue"))
  # + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
  #                      labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Strength", title = Cont.Name_Fig(ii, jj))
)

print(Fig_Cont.Strength, vp = vplayout(ii, jj))

} }


# Figure for Continent Weights with Positive Links only (3*3 Format)
grid.newpage()
pushViewport( viewport(layout = grid.layout(3, 3)) )
for(ii in 1:3){ for(jj in 1:3){

kk = 3*(ii - 1) + jj

DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, kk])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1)
  #+ guides(color = "none")
  + ylim(0, 16)
  # + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
  #                      labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)

print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))

} }



# Figure for Continent Weights with Negative Links only (3*3 Format)
grid.newpage()
pushViewport( viewport(layout = grid.layout(3, 3)) )
for(ii in 1:3){ for(jj in 1:3){

kk = 3*(ii - 1) + jj

DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, kk])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1)
  #+ guides(color = "none") + ylim(0, 20)
  # + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
  #                      labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)

print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))

} }





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Figure for Continent Weights with Positive Links only (3*3 Format)              ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





grid.newpage()
pushViewport( viewport(layout = grid.layout(3, 3)) )

{
# As to As
ii = 1
jj = 1
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# As to Eu
ii = 1
jj = 2
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 15)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# As to Eu
ii = 1
jj = 3
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Eu to As
ii = 2
jj = 1
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Eu to Eu
ii = 2
jj = 2
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Eu to Am
ii = 2
jj = 3
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Am to As
ii = 3
jj = 1
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Am to As
ii = 3
jj = 1
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Am to Eu
ii = 3
jj = 2
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))


# Am to Am
ii = 3
jj = 3
DF_Cont.Strength_P = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_P[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_P = ggplot(
  data = DF_Cont.Strength_P, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_P = (
  BFig_Cont.Strength_P + geom_line(color = "red", linetype = "solid", size = 1.5)
  + ylim(0, 8)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Positive Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_P, vp = vplayout(ii, jj))

}





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Figure for Continent Weights with Negative Links only (3*3 Format)              ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





grid.newpage()
pushViewport( viewport(layout = grid.layout(3, 3)) )

{
# As to As
ii = 1
jj = 1
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# As to Eu
ii = 1
jj = 2
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# As to Am
ii = 1
jj = 3
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Eu to As
ii = 2
jj = 1
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Eu to Eu
ii = 2
jj = 2
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 15)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Eu to Am
ii = 2
jj = 3
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Am to As
ii = 3
jj = 1
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Am to Eu
ii = 3
jj = 2
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))


# Am to Am
ii = 3
jj = 3
DF_Cont.Strength_N = data.frame(
  Time = 1:Phase.Num, Value = c(Continent.Weight_N[, 3*(ii - 1) + jj])
)
BFig_Cont.Strength_N = ggplot(
  data = DF_Cont.Strength_N, mapping = aes(x = Time, y = Value)
)
Fig_Cont.Strength_N = (
  BFig_Cont.Strength_N + geom_line(color = "blue", linetype = "solid", size = 1.5)
  + ylim(0, 4)
  + scale_x_continuous(breaks = round(nrow(Density)/40*c(2.5,13,23,34.5)),
                       labels = c("2006", "2007", "2008", "2009"))
  + labs(x = "Date", y = "Negative Strength", title = Cont.Name_Fig(ii, jj))
)
print(Fig_Cont.Strength_N, vp = vplayout(ii, jj))

}







##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##               Trail                                                                 ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





time.start = Sys.time()

# Define Variables
Roll.Win = 100
New.h = 10
XX = LogReturns[(Period_Start+(44-1)*New.h):(Period_Start-1+Roll.Win+(46-1)*New.h), ]
N = ncol(XX); Period = nrow(XX); Max.Link = N*(N-1)
X = scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1) #Centralized and Standardized
auto_term_loc = rep(1, N)
Adjacency_Period = matrix(0, nrow = N, ncol = N)
Weight_Period    = matrix(0, nrow = N, ncol = N)

Num_Cores  = detectCores(logical = TRUE) # - 1
chunk.size = Rep.Times/Num_Cores


# Divide Data into Different Areas
Asia_area     =  1:11;  Asia    = X[, Asia_area];    N_Asia    = ncol(Asia)
Europe_area   = 12:30;  Europe  = X[, Europe_area];  N_Europe  = ncol(Europe)
America_area  = 31:36;  America = X[, America_area]; N_America = ncol(America)


# Repeat to Generate Adajacency and Weight Matrices
pb = txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

## Define Asian, Europe and America Section
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

## DoParallel
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

## Generate the Freqency of Lambdas
Rep.Lambdas = CV_M[, 1]
Frq.Lambdas = matrix( c( as.numeric(names(table(Rep.Lambdas))),
                         as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
if(ncol(Frq.Lambdas) <= PreSet.M) {Set.M = ncol(Frq.Lambdas)}
if(ncol(Frq.Lambdas) >  PreSet.M) {Set.M = PreSet.M}
Lambdas_Top5_loc = sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 = Frq.Lambdas[1, Lambdas_Top5_loc]

## Generate the Adjacency Matrix
Adjacency_loc = which(Rep.Lambdas < max(Lambdas_Top5)+1e-9 
                      & Rep.Lambdas > max(Lambdas_Top5)-1e-9)[1]
Adjacency_Period[, k] = sign( CV_M[Adjacency_loc, 2:(N+1)] )

## Generate the Weight Matrix
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

time.end = Sys.time()
print(time.end - time.start)


sum( abs(Adjacency_Period) )/Max.Link
sum( abs(Adjacency_Period) - diag(diag(abs(Adjacency_Period))) )/Max.Link
sum( (abs(Adjacency_Period)+Adjacency_Period)/2 )/Max.Link
sum( (abs(Adjacency_Period)-Adjacency_Period)/2 )/Max.Link





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##






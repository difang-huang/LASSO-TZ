##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Tables in Time-Zone VAR Model by Optimal Lambda AR                             ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





setwd("C:/Users/581/Documents/Optimal LASSO Networks")

# Input adjacency matrix over full period
FULL.Adj.M = read.csv(file = "Two Sequent Crisis-FA_1e-4 Replications.csv", header = FALSE)
FULL.Adj.M = as.matrix(FULL.Adj.M)

# Input weight matrix over full period
FULL.Wet.M = read.csv(file = "Two Sequent Crisis-FW_1e-4 Replications.csv", header = FALSE)
FULL.Wet.M = as.matrix(FULL.Wet.M)


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


N = length(Countries)
Phase.Num = nrow(FULL.Adj.M)/N


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggthemes)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Table: AR Coefficients over Five Periods                                   ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# AR Coefficients
NewTab_AR.Coefs = matrix(0, nrow = N, ncol = Phase.Num)
for(phase in 1:Phase.Num){
  Weight = FULL.Wet.M[(phase*N-N+1):(phase*N), ]
  NewTab_AR.Coefs[, phase] = diag(Weight)
}
NewTab_AR.Coefs = cbind(Countries, NewTab_AR.Coefs)
write.csv(NewTab_AR.Coefs,
          paste(StoragePath, "NewTable-AR Coefficients.csv", sep = ""),
          quote = FALSE)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Table: Continent and Degree Assortativity                                  ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Assortativity and Disassortativity
Assort.Continent = matrix(0, nrow = 3, ncol = Phase.Num)
Assort.Degree    = matrix(0, nrow = 3, ncol = Phase.Num)
for(phase in 1:Phase.Num){
  Continents_Lables = c(rep(1, N_Asia), rep(2, N_Europe), rep(3, N_America))
  Adjacency = FULL.Adj.M[(phase*N-N+1):(phase*N), ]
  diag(Adjacency) = rep(0, N); # diag(Weights) = rep(0, N)
  Adj.M_P   = (abs(Adjacency) + Adjacency)/2
  Adj.M_N   = (abs(Adjacency) - Adjacency)/2
  
  A_U = graph_from_adjacency_matrix(abs(Adjacency), mode = "directed")
  Assort.Continent[1, phase] = assortativity(A_U, Continents_Lables, directed = TRUE)
  Assort.Degree[1, phase]    = assortativity_degree(A_U, directed = TRUE)
  
  A_P = graph_from_adjacency_matrix(Adj.M_P, mode = "directed")
  Assort.Continent[2, phase] = assortativity(A_P, Continents_Lables, directed = TRUE)
  Assort.Degree[2, phase]    = assortativity_degree(A_P, directed = TRUE)
  
  A_N = graph_from_adjacency_matrix(Adj.M_N, mode = "directed")
  Assort.Continent[3, phase] = assortativity(A_N, Continents_Lables)
  Assort.Degree[3, phase]    = assortativity_degree(A_N, directed = TRUE)
}
NewTab_Assort = rbind(Assort.Continent, Assort.Degree)
View(NewTab_Assort)
write.csv(NewTab_Assort,
          paste(StoragePath, "NewTable-Assortativity of Continent and Degree.csv", sep = ""),
          quote = FALSE)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Table: Continent Degrees and Strengths                                     ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Function to calculate continent degrees and strengths
Continent.Sum = function(adjacency, weight, type = c("U", "P", "N"),
                         output.type = c("adjacency","weight"), Asia_area = 1:11,
                         Europe_area = 12:30, America_area = 31:36){

if(type == "U"){ ct.adjacency = abs(adjacency) }                # Unsign
if(type == "P"){ ct.adjacency = (abs(adjacency)+adjacency)/2 }  # Positive
if(type == "N"){ ct.adjacency = (abs(adjacency)-adjacency)/2 }  # Negative
ct.weight = ct.adjacency * abs(weight)
# diag(ct.adjacency) = rep(0, nrow(adjacency)) ????????ǰ???Խ???Ԫ??????Ϊ0

if(output.type == "adjacency"){
  Cont.Adjacency = matrix(0, nrow = 3, ncol = 3)
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

if(output.type == "weight"){
  Cont.Weight    = matrix(0, nrow = 3, ncol = 3)
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
##=====================================================================================##


# Degree and strength analysis
NewTab_Cont.Adj = matrix(0, nrow = Phase.Num*3, ncol = 3*3)
NewTab_Cont.Wet = matrix(0, nrow = Phase.Num*3, ncol = 3*3)
for(phase in 1:Phase.Num){

Adjacency = FULL.Adj.M[(phase*N-N+1):(phase*N), ]
Weight    = FULL.Wet.M[(phase*N-N+1):(phase*N), ]
diag(Adjacency) = rep(0, N); diag(Weight) = rep(0, N)

# Analysis on Continent Degrees
CT.A_U = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "U",
                       output.type = "adjacency")
CT.A_P = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "P",
                       output.type = "adjacency")
CT.A_N = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "N",
                       output.type = "adjacency")
NewTab_Cont.Adj[(3*phase-2):(3*phase), ] = cbind(CT.A_U, CT.A_P, CT.A_N)

# Analysis on Continent Strengths
CT.W_U = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "U",
                       output.type = "weight")
CT.W_P = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "P",
                       output.type = "weight")
CT.W_N = Continent.Sum(adjacency = Adjacency, weight = Weight, type = "N",
                       output.type = "weight")
NewTab_Cont.Wet[(3*phase-2):(3*phase), ] = cbind(CT.W_U, CT.W_P, CT.W_N)

}
View(NewTab_Cont.Wet)
write.csv(NewTab_Cont.Adj,
          paste(StoragePath, "NewTable-Continent Degree Analysis.csv", sep = ""), 
          quote = FALSE)
write.csv(NewTab_Cont.Wet,
          paste(StoragePath, "NewTable-Continent Strength Analysis.csv", sep = ""),
          quote = FALSE)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Table: Summary statistics of national equity market indexes                ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Output Phase Names
Phs.Func = function(phs, Rep.Num = N){
  if(phs == 1){ return(rep("Before Subprime Crisis", Rep.Num)) }
  if(phs == 2){ return(rep("During Subprime Crisis", Rep.Num)) }
  if(phs == 3){ return(rep("After Subprime Crisis",  Rep.Num)) }
  if(phs == 4){ return(rep("During European Debts",  Rep.Num)) }
  if(phs == 5){ return(rep("After European Debts",   Rep.Num)) }
}


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Input Return Data
LogReturns = read.csv(
  file = "D:/Data/National Stock Market Indeces/Cleaned_Log_Returns_LongPeriod.csv",
  header = FALSE
  )
LogReturns = as.matrix(LogReturns)
BeginT = 34946   # 1995/9/4
EndT   = 44890   # 2022/11/25
Alldays = BeginT:EndT
Alldays_m = c( Alldays, rep(0, 7-length(Alldays)%%7) )
Weekdays_m = matrix(Alldays_m, nrow = 7, byrow = FALSE)
Num_Weekdays = 5*(ncol(Weekdays_m)-1) + length(Alldays)%%7
Weekdays_v = as.vector(Weekdays_m[1:5, ])[1:Num_Weekdays]
# ?δ?ǰ->?δ???->?δ???->ŷծ??->ŷծ??
Period_Range = c(38930,39294, 39295,39903, 39904,40147,   40148,40624, 41625,42369)


# Summary Statistics
NewTab_Summary.Statistics = matrix(0, nrow = Phase.Num*N, ncol = 7)
for(phase in 1:Phase.Num){

Phase_Start = which(Weekdays_v == Period_Range[2*phase-1])
Phase_End   = which(Weekdays_v == Period_Range[2*phase])
XX = LogReturns[Phase_Start:Phase_End, ]

Summary.Stats = matrix(0, nrow = N, ncol = 7)
Summary.Stats[, 1] = apply(XX, 2, mean)
Summary.Stats[, 2] = apply(XX, 2, sd)
Summary.Stats[, 3:7] = apply(XX, 2, fivenum)
NewTab_Summary.Statistics[(phase*N-N+1):(phase*N), ] = Summary.Stats

}
write.csv(NewTab_Summary.Statistics,
          paste(StoragePath, "NewTable-Summary Statistics of Index Returns.csv", sep = ""), 
          quote = FALSE)








##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##                              Draft                                                  ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##



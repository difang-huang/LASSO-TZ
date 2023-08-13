setwd("D:/Data/National Stock Market Indeces")
LogReturns <- 
  read.csv(file = "D:/Data/National Stock Market Indeces/Cleaned_Log_Returns.csv", 
           header = FALSE)
LogReturns <- as.matrix(LogReturns)


BeginT <- 34946   # 1995/9/4
EndT   <- 42587   # 2016/8/5

Alldays <- BeginT:EndT
Alldays_m <- c( Alldays, rep(0, 7-length(Alldays)%%7) )
Weekdays_m <- matrix(Alldays_m, nrow = 7, byrow = FALSE)
Num_Weekdays <- 5*(ncol(Weekdays_m)-1) + length(Alldays)%%7
Weekdays_v <- as.vector(Weekdays_m[1:5, ])[1:Num_Weekdays]

# Country Names
{
Countries <- c("Australia","Malaysia","Indonesia","Korea","Japan","Singapore",
               "NewZealand","Philippines","Thailand","China","HK",
               
               "Holland","Greece","Belgium","France","Germany","Finland","Spain",
               "Ireland","Italy","Denmark","Norway","Sweden","Portugal","Russia",
               "Switzerland","UK","Poland","Turkey","Austria",
               
               "Brazil","Chile","Argentina","Mexico","Canada","US")
}



##=====================================================================================##
##=====================================================================================##
##               Generate Adjacency and Weighted Matrices                              ##
##=====================================================================================##
##=====================================================================================##



# 次贷前->次贷中->次贷后->欧债
Period_Range <-c(38930, 39294, 39295, 39903, 39904, 40147, 40148, 40543)
# 次贷前->次贷中1->次贷中2->次贷后->欧债
#Period_Range <-c(38930, 39294, 39295, 39706, 39707, 39903, 39904, 40147, 40148, 40543)

Phase_Num <- length(Period_Range)/2
nlinks_Unsign <- rep(0, Phase_Num)
nlinks_counter <- rep(0, Phase_Num)
AR_Terms <- matrix(0, nrow = length(Countries), ncol = Phase_Num+1)
AR_Terms[, 1] <- Countries
Adaptive.AR <- TRUE
Rep.Times <- 3000#0
Set.M <- 10
Final.Adjacency <- matrix(0, nrow = Phase_Num*length(Countries), ncol = length(Countries))
Final.Weight    <- matrix(0, nrow = Phase_Num*length(Countries), ncol = length(Countries))
StoragePath <- "C:/Users/581/Documents/Optimal LASSO Networks/2008 Subprime Crisis by"


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)

for(phase in 1:Phase_Num){

# Define Variables
{
Phase_Start  <- which(Weekdays_v == Period_Range[2*phase-1])
Phase_End    <- which(Weekdays_v == Period_Range[2*phase])
XX <- LogReturns[Phase_Start:Phase_End, ]
N <- ncol(XX)  #Define the dimension of the weight matrix
Period <- nrow(XX)
MaxLinks <- N*(N-1)
print(paste("Phase ", phase, " has ", Period, " observations.", sep = ""))
  
# Centralized and Standardized
X <- scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1)
#diag(t(X) %*% X)
  
# Divide Data into Different Areas
Asia_area     <-  1:11;  Asia    <- X[, Asia_area];    N_Asia    <- ncol(Asia)
Europe_area   <- 12:30;  Europe  <- X[, Europe_area];  N_Europe  <- ncol(Europe)
America_area  <- 31:36;  America <- X[, America_area]; N_America <- ncol(America)

Num_kfolds <- floor(Period/30)
auto_term_loc <- rep(1, N)
Adjacency_Period <- matrix(0, nrow = N, ncol = N)
Weight_Period    <- matrix(0, nrow = N, ncol = N)

Num_Cores <- detectCores(logical = FALSE) # - 1
chunk.size <- Rep.Times/Num_Cores
}


# Repeat to Generate Adajacency and Weight Matrices
pb <- txtProgressBar(min = 0, max = N, style = 3)
for(k in 1:N){

# Define Asian, Europe and America Section
if(k >= 1                 & k <= N_Asia){
    y <- Asia[2:Period, k]
    x <- cbind(Asia[1:(Period-1),], Europe[1:(Period-1),], America[1:(Period-1),])
  }
if(k >= N_Asia+1          & k <= N_Asia+N_Europe){
    y <- Europe[2:Period, k-N_Asia]
    x <- cbind(Asia[2:Period,], Europe[1:(Period-1),], America[1:(Period-1),])
  }
if(k >= N_Asia+N_Europe+1 & k <= N){
    y <- America[2:Period, k-N_Asia-N_Europe]
    x <- cbind(Asia[2:Period,], Europe[2:Period,], America[1:(Period-1),])
  }
if(Adaptive.AR == FALSE){auto_term_loc[k] <- 0}

# DoParallel
cl <- makeCluster(Num_Cores)
registerDoParallel(cl, cores = Num_Cores)
CV_M <- foreach(core = 1:Num_Cores, .combine = "rbind", .packages = "glmnet") %dopar%
{
  Lambda_Est <- matrix(0, nrow = chunk.size, ncol = N+1)
  for(rr in 1:chunk.size){
    cvfit <- cv.glmnet(sqrt(Period-1)*x, sqrt(Period-1)*y, nfolds = Num_kfolds, 
                       alpha = 1, type.measure = "mse", nlambda = 300, 
                       intercept = FALSE, standardize = FALSE, thresh = 1e-7, 
                       penalty.factor = auto_term_loc, parallel = FALSE)
    Lambda_Est[rr, 1] <- cvfit$lambda.min
    Lambda_Est[rr, 2:(N+1)] <- as.vector( coef(cvfit, s="lambda.min") )[-1]
  }
  return(Lambda_Est)
}
stopImplicitCluster(); stopCluster(cl)

# Generate the Freqency of Lambdas
Rep.Lambdas <- CV_M[, 1]
Frq.Lambdas <- matrix( c( as.numeric(names(table(Rep.Lambdas))),
                          as.vector(table(Rep.Lambdas)) ), nrow = 2, byrow = TRUE )
Lambdas_Top5_loc <- sort.list(Frq.Lambdas[2, ], decreasing = TRUE)[1:Set.M]
Lambdas_Top5 <- Frq.Lambdas[1, Lambdas_Top5_loc]

# Generate the Adjacency Matrix
Adjacency_loc <- which(Rep.Lambdas < max(Lambdas_Top5)+1e-9 
                       & Rep.Lambdas > max(Lambdas_Top5)-1e-9)[1]
Adjacency_Period[, k] <- sign( CV_M[Adjacency_loc, 2:(N+1)] )

# Generate the Weight Matrix
Betas_Top5_w <- Frq.Lambdas[2, Lambdas_Top5_loc]/sum(Frq.Lambdas[2, Lambdas_Top5_loc])
Betas_Top5_coefs <- rep(0, N)
for(rr in 1:Set.M){
  Betas_loc <- which(Rep.Lambdas < Lambdas_Top5[rr]+1e-9 
                     & Rep.Lambdas > Lambdas_Top5[rr]-1e-9)[1]
  Betas_Top5_coefs <- Betas_Top5_coefs + Betas_Top5_w[rr] * CV_M[Betas_loc, 2:(N+1)]
}
Weight_Period[, k] <- Betas_Top5_coefs * abs(Adjacency_Period[, k])

Sys.sleep(0.1)
setTxtProgressBar(pb, k)
}
close(pb)

Final.Adjacency[(phase*N-N+1):(phase*N), ] <- Adjacency_Period
Final.Weight[(phase*N-N+1):(phase*N), ]    <- Weight_Period


# Analysis on Autogressive Terms and Link Numbers
AR_Terms[, phase+1] <- diag(Weight_Period)
# diag(Adjacency_Period) <- rep(0, N)
nlinks_Unsign[phase] <- sum( abs(Adjacency_Period) )
nlinks_counter[phase] <- sum( Adjacency_Period )


print(paste("The Phase ", phase, " has completed.", sep = ""))

}
nlinks_Unsign; nlinks_counter
View(AR_Terms)


StoragePath <- "C:/Users/581/Documents/Optimal LASSO Networks/2008 Subprime Crisis by"
StoragePath_FA <- paste(StoragePath, "AR Opitmal Lambdas-FA.csv")
write.csv(Final.Adjacency,  StoragePath_FA,  quote = FALSE)
StoragePath_FW <- paste(StoragePath, "AR Opitmal Lambdas-FW.csv")
write.csv(Final.Weight,  StoragePath_FW,  quote = FALSE)



##=====================================================================================##
##=====================================================================================##
##               Degrees and Strengths Analysis                                        ##
##=====================================================================================##
##=====================================================================================##



# Define the Degree Output Function
Ranks_Output <- function(DG, LastPosition = 36){
  Nations <- Countries[sort.list(DG, decreasing = TRUE)][1:LastPosition]
  Values  <- sort(DG, decreasing = TRUE)[1:LastPosition]
  return( matrix(c(Nations,Values), nrow = LastPosition, ncol = 2, byrow = FALSE) )
}

# Define the Continental Adjacency Function
Continents_Sum <- function(adjacency, weight, type = c("U", "P", "N"),
                           output.adjacency = TRUE, Asia_area = 1:11,
                           Europe_area = 12:30, America_area = 31:36){
  Cont.Adjacency = matrix(0, nrow = 3, ncol = 3)
  Cont.Weight    = matrix(0, nrow = 3, ncol = 3)
  if(type == "U"){ ct.adjacency = abs(adjacency) }               # Unsign
  if(type == "P"){ ct.adjacency = (abs(adjacency)+adjacency)/2 } # Positive
  if(type == "N"){ ct.adjacency = (abs(adjacency)-adjacency)/2 } # Negative
  ct.weight = ct.adjacency * abs(weight)
  # diag(ct.adjacency) <- rep(0, nrow(adjacency))
  
  Cont.Adjacency[1, 1] = sum( ct.adjacency[Asia_area, Asia_area] )
  Cont.Adjacency[1, 2] = sum( ct.adjacency[Asia_area, Europe_area] )
  Cont.Adjacency[1, 3] = sum( ct.adjacency[Asia_area, America_area] )
  Cont.Adjacency[2, 1] = sum( ct.adjacency[Europe_area, Asia_area] )
  Cont.Adjacency[2, 2] = sum( ct.adjacency[Europe_area, Europe_area] )
  Cont.Adjacency[2, 3] = sum( ct.adjacency[Europe_area, America_area] )
  Cont.Adjacency[3, 1] = sum( ct.adjacency[America_area, Asia_area] )
  Cont.Adjacency[3, 2] = sum( ct.adjacency[America_area, Europe_area] )
  Cont.Adjacency[3, 3] = sum( ct.adjacency[America_area, America_area] )
  
  Cont.Weight[1, 1] = sum( ct.weight[Asia_area, Asia_area] )
  Cont.Weight[1, 2] = sum( ct.weight[Asia_area, Europe_area] )
  Cont.Weight[1, 3] = sum( ct.weight[Asia_area, America_area] )
  Cont.Weight[2, 1] = sum( ct.weight[Europe_area, Asia_area] )
  Cont.Weight[2, 2] = sum( ct.weight[Europe_area, Europe_area] )
  Cont.Weight[2, 3] = sum( ct.weight[Europe_area, America_area] )
  Cont.Weight[3, 1] = sum( ct.weight[America_area, Asia_area] )
  Cont.Weight[3, 2] = sum( ct.weight[America_area, Europe_area] )
  Cont.Weight[3, 3] = sum( ct.weight[America_area, America_area] )
  
  if(output.adjacency == TRUE){  return(Cont.Adjacency) }
  if(output.adjacency == FALSE){ return(Cont.Weight) }
}

# Output the Continental Ranges
CT.Ranges <- function(m){
  if(m == 1){return( 1:N_Asia )}
  if(m == 2){return( (N_Asia+1):(N_Asia+N_Europe) )}
  if(m == 3){return( (N_Asia+N_Europe+1):N )}
}

# Output Phase Names
Phs.Func <- function(phs, Rep.Num = N){
  if(phs == 1){ return(rep("Before", Rep.Num)) }
  if(phs == 2){ return(rep("During", Rep.Num)) }
  if(phs == 3){ return(rep("After",  Rep.Num)) }
  if(phs == 4){ return(rep("Europe Debts", Rep.Num)) }
}

##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Degree Analysis
Degree_Unsign <- matrix(0, nrow = 30, ncol = Phase_Num*6)
Degree_Offset <- matrix(0, nrow = 30, ncol = Phase_Num*6)
CT.Adjacency = matrix(0, nrow = Phase_Num*3, ncol = 3*3)
CT.Weight    = matrix(0, nrow = Phase_Num*3, ncol = 3*3)
for(phase in 1:Phase_Num){
  Adjacency = Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weight    = Final.Weight[(phase*N-N+1):(phase*N), ]
  diag(Adjacency) = rep(0, N); diag(Weight) = rep(0, N)
  
  # Unsigns
  InD_U  <- apply(abs(Adjacency), 2, sum)
  OutD_U <- apply(abs(Adjacency), 1, sum)
  InD_Unsign  <- Ranks_Output(InD_U, 30)
  OutD_Unsign <- Ranks_Output(OutD_U, 30)
  TotalD_Unsign <- Ranks_Output(InD_U+OutD_U, 30)
  Degree_Unsign[, (6*phase-5):(6*phase)] <- cbind(InD_Unsign, OutD_Unsign, TotalD_Unsign)
  
  # Offset
  InD_O  <- apply(Adjacency, 2, sum)
  OutD_O <- apply(Adjacency, 1, sum)
  InD_Offset    <- Ranks_Output(InD_O, 30)
  OutD_Offset   <- Ranks_Output(OutD_O, 30)
  TotalD_Offset <- Ranks_Output(InD_O+OutD_O, 30)
  Degree_Offset[, (6*phase-5):(6*phase)] <- cbind(InD_Offset, OutD_Offset,TotalD_Offset)
  
  # Analysis on Continent Degrees
  CT.A_U = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "U",
                          output.adjacency = TRUE)
  CT.A_P = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "P",
                          output.adjacency = TRUE)
  CT.A_N = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "N",
                          output.adjacency = TRUE)
  CT.Adjacency[(3*phase-2):(3*phase), ] <- cbind(CT.A_U, CT.A_P, CT.A_N)
  
  # Analysis on Continent Strengths
  CT.W_U = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "U",
                          output.adjacency = FALSE)
  CT.W_P = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "P",
                          output.adjacency = FALSE)
  CT.W_N = Continents_Sum(adjacency = Adjacency, weight = Weight, type = "N",
                          output.adjacency = FALSE)
  CT.Weight[(3*phase-2):(3*phase), ] <- cbind(CT.W_U, CT.W_P, CT.W_N)
}
View(Degree_Unsign)
View(Degree_Offset)
View(CT.Adjacency)
View(CT.Weight)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Strength Summary Analysis
library(ggplot2)
RBG <- c("red", "blue", "green")
RBG <- setNames(RBG, c("Asian", "Europe", "American"))
CT_List <- c(rep("Asian",N_Asia), rep("Europe",N_Europe), rep("American",N_America))

Strength_Sumry  <- matrix(0, nrow = 36, ncol = Phase_Num*4)
S <- matrix(0, nrow = N, ncol = 4*Phase_Num)
S_Spain  <- matrix(0, nrow = N, ncol = 4*Phase_Num)
S_Turkey <- matrix(0, nrow = N, ncol = 4*Phase_Num)
for(phase in 1:Phase_Num){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  diag(Weights) <- rep(0, N)
  Adj.M_P   <- (abs(Adjacency) + Adjacency)/2; diag(Adj.M_P) <- rep(0, N)
  Adj.M_N   <- (abs(Adjacency) - Adjacency)/2; diag(Adj.M_N) <- rep(0, N)
  
  # Spain
  S_Spain[, 4*phase-3] <- (Adj.M_P*Weights)[, 18]
  S_Spain[, 4*phase-2] <- (Adj.M_N*Weights)[, 18]
  S_Spain[, 4*phase-1] <- (Adj.M_P*Weights)[18, ]
  S_Spain[, 4*phase  ] <- (Adj.M_N*Weights)[18, ]
  
  # Turkey
  S_Turkey[, 4*phase-3] <- (Adj.M_P*Weights)[, 29]
  S_Turkey[, 4*phase-2] <- (Adj.M_N*Weights)[, 29]
  S_Turkey[, 4*phase-1] <- (Adj.M_P*Weights)[29, ]
  S_Turkey[, 4*phase  ] <- (Adj.M_N*Weights)[29, ]
  
  # Strength Summary
  InS_P  <- apply(abs(Adj.M_P*Weights), 2, sum) 
  InS_N  <- apply(abs(Adj.M_N*Weights), 2, sum)
  OutS_P <- apply(abs(Adj.M_P*Weights), 1, sum)
  OutS_N <- apply(abs(Adj.M_N*Weights), 1, sum)
  Strength_Sumry[, 4*phase-3] <- InS_P;   Strength_Sumry[, 4*phase-2] <- InS_N
  Strength_Sumry[, 4*phase-1] <- OutS_P;  Strength_Sumry[, 4*phase  ] <- OutS_N
  
  # Draw Strengths Scattter Plots Phase by Phase
  {
    # Net In-Strengths and Out-Strengths
    Net.Strengths <- data.frame(Agent = Countries, Continent = CT_List,
                                Net.InS  = InS_P-InS_N, Net.OutS = OutS_P-OutS_N)
    BFig_Net.S <- ggplot(data = Net.Strengths,
                         mapping = aes(x = Net.OutS, y = Net.InS, color = Continent))
    if(phase == 1){
      Fig_Net.S_1 <- ( 
        BFig_Net.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        # + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        # + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-2,5) + ylim(-0.1,0.9)
        + labs(title = "Before", x = "Net Out-Strengths", y = "Net In-Strengths") )
    }
    if(phase == 2){
      Fig_Net.S_2 <- (
        BFig_Net.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        # + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        # + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-2,5) + ylim(-0.1,0.9)
        + labs(title = "During", x = "Net Out-Strengths", y = "Net In-Strengths") )
    }
    if(phase == 3){
      Fig_Net.S_3 <- (
        BFig_Net.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        # + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        # + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-2,5) + ylim(-0.1,0.9)
        + labs(title = "After", x = "Net Out-Strengths", y = "Net In-Strengths") )
    }
    if(phase == 4){
      Fig_Net.S_4 <- (
        BFig_Net.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        # + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        # + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-2,5) + ylim(-0.1,0.9)
        + labs(title = "Europe Debts", x = "Net Out-Strengths", y = "Net In-Strengths") )
    }
    # color = NA
    
    # Strengths in Positive and Negative Sub Networks
    PN.Strengths <- data.frame(Agent = Countries, Continent = CT_List,
                               PoS  = OutS_P-InS_P, NeS = OutS_N-InS_N)
    BFig_PN.S <- ggplot(data = PN.Strengths,
                        mapping = aes(x = PoS, y = NeS, color = Continent))
    if(phase == 1){
      Fig_PN.S_1 <- ( 
        BFig_PN.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-1,4) + ylim(-1,1.5)
        + labs(title = "Before", x = "Positive Strengths", y = "Negative Strengths") )
    }
    if(phase == 2){
      Fig_PN.S_2 <- ( 
        BFig_PN.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-1,4) + ylim(-1,1.5)
        + labs(title = "During", x = "Positive Strengths", y = "Negative Strengths") )
    }
    if(phase == 3){
      Fig_PN.S_3 <- ( 
        BFig_PN.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-1,4) + ylim(-1,1.5)
        + labs(title = "After", x = "Positive Strengths", y = "Negative Strengths") )
    }
    if(phase == 4){
      Fig_PN.S_4 <- ( 
        BFig_PN.S + geom_point(size = 1)
        + scale_color_manual(values = RBG) + guides(color = "none")
        + geom_text(aes(label=Countries), size=3)
        + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
        + geom_abline(slope = -1, intercept = 0, linetype="dashed")
        + xlim(-1,4) + ylim(-1,1.5)
        + labs(title = "Europe Debts", x = "Positive Strengths", y = "Negative Strengths") )
    }
  }
}
View(Strength_Sumry)

View(S_Spain); View(S_Turkey)
Fig_Net.S_1; Fig_Net.S_2; Fig_Net.S_3; # Fig_Net.S_4
Fig_PN.S_1; Fig_PN.S_2; Fig_PN.S_3; # Fig_PN.S_4
Strength_Sumry_Path <- paste(StoragePath, "Opitmal AR Strength Summary.csv")
write.csv(Strength_Sumry, Strength_Sumry_Path, quote = FALSE)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Draw Scatter Plots for Strengths in 3*1 Forms
Net.Strengths <- data.frame(Agent = "aaa", Continent = "bbb", Net.InS  = 1, 
                            Net.OutS = 1, Phase.Now = "ccc")
PN.Strengths  <- data.frame(Agent = "aaa", Continent = "bbb", PoS  = 1, NeS = 1, 
                            Phase.Now = "ccc")
for(phase in 1:(Phase_Num-1)){
  # Strength Summary
  InS_P  <- Strength_Sumry[, 4*phase-3]
  InS_N  <- Strength_Sumry[, 4*phase-2]
  OutS_P <- Strength_Sumry[, 4*phase-1]
  OutS_N <- Strength_Sumry[, 4*phase  ]
  
  # Net In-Strengths and Out-Strengths
  Net.New <- data.frame(Agent = Countries, Continent = CT_List, Net.InS  = InS_P-InS_N,
                        Net.OutS = OutS_P-OutS_N, Phase.Now = Phs.Func(phase))
  Net.Strengths <- rbind(Net.Strengths, Net.New)

  # Strengths in Positive and Negative Sub Networks
  PN.New <- data.frame(Agent = Countries, Continent = CT_List, PoS  = OutS_P-InS_P, 
                       NeS = OutS_N-InS_N, Phase.Now = Phs.Func(phase))
  PN.Strengths <- rbind(PN.Strengths, PN.New)
}

# Net In- and Out- Strengths
Net.Strengths <- Net.Strengths[-1, ]
BFig_Net.S <- ggplot(data = Net.Strengths,
                     mapping = aes(x = Net.OutS, y = Net.InS, color = Continent))
Fig_Net.S <- ( 
  BFig_Net.S + geom_point(size = 1)
  + scale_color_manual(values = RBG) + guides(color = "none")
  + geom_text(aes(label=rep(Countries,(Phase_Num-1))), size=2)
  + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  # + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
  # + geom_abline(slope = -1, intercept = 0, linetype="dashed")
  + xlim(-2,5) + ylim(-0.1,0.9) + labs(title = "Net In- and Out- Strengths", 
         x = "Net Out-Strengths", y = "Net In-Strengths")
  )
Fig_Net.S + facet_grid(Phase.Now~.)

# Positive and Negative Strengths
PN.Strengths  <- PN.Strengths[-1, ]
BFig_PN.S <- ggplot(data = PN.Strengths,
                    mapping = aes(x = PoS, y = NeS, color = Continent))
Fig_PN.S <- ( 
  BFig_PN.S + geom_point(size = 1)
  + scale_color_manual(values = RBG) + guides(color = "none")
  + geom_text(aes(label=rep(Countries,(Phase_Num-1))), size=2)
  + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  + geom_abline(slope = 1,  intercept = 0, linetype="dashed")
  + geom_abline(slope = -1, intercept = 0, linetype="dashed")
  + xlim(-1,4) + ylim(-1,1.5) + labs(title = "Positive and Negative Strengths",
        x = "Positive Strengths", y = "Negative Strengths") 
  )
Fig_PN.S + facet_grid(Phase.Now~.)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Roles Detection
Roles.Detect <- function(Var.X, Var.Y, Details = FALSE){
  if(Var.X > 0 & Var.Y > 0){
    if(Details == FALSE){ return(1) }
    else{
      if(abs(Var.X) > abs(Var.Y)){ return(11) }
      if(abs(Var.X) < abs(Var.Y)){ return(12) }
    }
  }
  if(Var.X < 0 & Var.Y > 0){
    if(Details == FALSE){ return(2) }
    else{
      if(abs(Var.X) < abs(Var.Y)){ return(23) }
      if(abs(Var.X) > abs(Var.Y)){ return(24) }
    }
  }
  if(Var.X < 0 & Var.Y < 0){
    if(Details == FALSE){ return(3) }
    else{
      if(abs(Var.X) > abs(Var.Y)){ return(35) }
      if(abs(Var.X) < abs(Var.Y)){ return(36) }
    }
  }
  if(Var.X > 0 & Var.Y < 0){
    if(Details == FALSE){ return(4) }
    else{
      if(abs(Var.X) > abs(Var.Y)){ return(47) }
      if(abs(Var.X) < abs(Var.Y)){ return(48) }
    }
  }
  
  if(Var.X >  0 & Var.Y == 0){ return(1148) }
  if(Var.X == 0 & Var.Y >  0){ return(1223) }
  if(Var.X <  0 & Var.Y == 0){ return(2435) }
  if(Var.X == 0 & Var.Y <  0){ return(3647) }
}

Roles_Net <- matrix(0, nrow = length(Countries), ncol = 1+Phase_Num)
Roles_PN  <- matrix(0, nrow = length(Countries), ncol = 1+Phase_Num)
Roles_Net[, 1] <- Countries; Roles_PN[, 1] <- Countries

for(phase in 1:Phase_Num){
  InS_P  <- Strength_Sumry[, 4*phase-3]
  InS_N  <- Strength_Sumry[, 4*phase-2]
  OutS_P <- Strength_Sumry[, 4*phase-1]
  OutS_N <- Strength_Sumry[, 4*phase  ]
  
  Net.InS <- InS_P - InS_N;  Net.OutS <- OutS_P - OutS_N
  Net.PoS <- OutS_P - InS_P; Net.NeS  <- OutS_N - InS_N
  
  for(kk in 1:length(Countries)){
    Roles_Net[kk, 1+phase] <- Roles.Detect(Net.OutS[kk], Net.InS[kk], Details = FALSE)
    Roles_PN[kk,  1+phase] <- Roles.Detect(Net.PoS[kk],  Net.NeS[kk], Details = TRUE)
  }
}
View(Roles_Net); View(Roles_PN)



##=====================================================================================##
##=====================================================================================##
##               Output Nodes and Links Lists for Gephi                                ##
##=====================================================================================##
##=====================================================================================##



# Define the Continent Detection Function
AsEuAm.Detect <- function(k){
  if(k <= N_Asia){ return("Asia") }
  if(k >= N_Asia+1 && k <= N_Asia+N_Europe){ return("Europe") }
  if(k >= N_Europe+1){ return("America") }
}


StPath_Gephi_root <- "C:/Users/581/Documents/Optimal LASSO Networks/AR Opitmal"


# Output Links Lists for Gephi
for(phase in 1:Phase_Num){
Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
diag(Adjacency) <- rep(0, N); diag(Weights) <- rep(0, N)

Links_T <- c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
Links_P <- c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
Links_N <- c("Source", "Target", "Type", "Signs", "Com.Weights", "Continent.From")
for(i in 1:N){
for(j in 1:N){
  judgement <- Adjacency[i, j]
  if(judgement == 1){
    New.Link <- c(Countries[i], Countries[j],"Directed", "P", 
                  abs(Weights[i, j]), AsEuAm.Detect(i))
    Links_T <- rbind(Links_T, New.Link); Links_P <- rbind(Links_P, New.Link)
  }
  if(judgement == -1){
    New.Link <- c(Countries[i], Countries[j],"Directed", "N", 
                  abs(Weights[i, j]), AsEuAm.Detect(i))
    Links_T <- rbind(Links_T, New.Link); Links_N <- rbind(Links_N, New.Link)
  }
}
}

StPath_Gephi_T <- paste(StPath_Gephi_root,
                        " Lambda-Full Graph Phase ", phase, ".csv", sep = "")
write.csv(Links_T, StPath_Gephi_T, quote = FALSE)
StPath_Gephi_P <- paste(StPath_Gephi_root,
                        " Lambda-Sub Graph Phase ", phase, "-P", ".csv", sep = "")
write.csv(Links_P, StPath_Gephi_P, quote = FALSE)
StPath_Gephi_N <- paste(StPath_Gephi_root,
                        " Lambda-Sub Graph Phase ", phase, "-N", ".csv", sep = "")
write.csv(Links_N, StPath_Gephi_N, quote = FALSE)
}


# Output Nodes Lists for Gephi
for(phase in 1:Phase_Num){
InS_P  <- Strength_Sumry[, 4*phase-3]
InS_N  <- Strength_Sumry[, 4*phase-2]
OutS_P <- Strength_Sumry[, 4*phase-1]
OutS_N <- Strength_Sumry[, 4*phase  ]

Nodes_List <- c("Id", "Label", "Continent", "U-Strength", "P-Strength", "N-Strength")
for(i in 1:N){
  New.Node <- c(Countries[i], Countries[i], AsEuAm.Detect(i),
                (InS_P+InS_N+OutS_P+OutS_N)[i], (InS_P+OutS_P)[i], (InS_N+OutS_N)[i])
  Nodes_List <- rbind(Nodes_List, New.Node)
}
  
StPath_Gephi_Nodes <- paste(StPath_Gephi_root,
                            " Lambda-Agents Info-Phase ", phase, ".csv", sep = "")
write.csv(Nodes_List, StPath_Gephi_Nodes, quote = FALSE)
}


# Output Continental Graphs for Gephi
for(phase in 1:Phase_Num){
Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]

for(ii in 1:3){
for(jj in 1:3){
  Com.Weights <- "Com.Weights"
  Links_List <- c("Source", "Target", "Type", "Signs")
  for(i in CT.Ranges(ii)){
  for(j in CT.Ranges(jj)){
    judgement <- Adjacency[i, j]
    if(judgement == 1){
      Links_List<-rbind(Links_List,c(Countries[i],Countries[j],"Directed","P"))
      Com.Weights <- c(Com.Weights, Weights[i, j])
    }
    if(judgement == -1){
      Links_List<-rbind(Links_List,c(Countries[i],Countries[j],"Directed","N"))
      Com.Weights <- c(Com.Weights, Weights[i, j])
    }
  }
  }
  StoragePath_SubGephi <- paste(StoragePath_Gephi_root, " Lambdas-SubGraph Phase ",
                                phase, " Part ", ii, jj, ".csv", sep = "")
  Gephi_Results <- cbind(Links_List[-1, ], Com.Weights[-1])
  write.csv(Gephi_Results, StoragePath_SubGephi, quote = FALSE)
}
}

}



##=====================================================================================##
##=====================================================================================##
##               Topological Properties Analysis                                       ##
##=====================================================================================##
##=====================================================================================##



# Assortativity and Disassortativity
Mix_Continents <- matrix(0, nrow = 3, ncol = Phase_Num)
Mix_Degrees    <- matrix(0, nrow = 3, ncol = Phase_Num)
for(phase in 1:Phase_Num){
  Continents_Lables <- c(rep(1, N_Asia), rep(2, N_Europe), rep(3, N_America))
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  mmP <- (abs(Adjacency) + Adjacency)/2
  mmN <- (abs(Adjacency) - Adjacency)/2
  
  A_U <- graph_from_adjacency_matrix(abs(Adjacency), mode = "directed")
  Mix_Continents[1, phase] <- assortativity(A_U, Continents_Lables, directed = TRUE)
  Mix_Degrees[1, phase]    <- assortativity_degree(A_U, directed = TRUE)
  
  A_P <- graph_from_adjacency_matrix(mmP, mode = "directed")
  Mix_Continents[2, phase] <- assortativity(A_P, Continents_Lables, directed = TRUE)
  Mix_Degrees[2, phase]    <- assortativity_degree(A_P, directed = TRUE)
  
  A_N <- graph_from_adjacency_matrix(mmN, mode = "directed")
  Mix_Continents[3, phase] <- assortativity(A_N, Continents_Lables)
  Mix_Degrees[3, phase]    <- assortativity_degree(A_N, directed = TRUE)
}
View(Mix_Continents); View(Mix_Degrees)
# U.Degrees_Labels <- apply(abs(Adjacency), 2, sum) + apply(abs(Adjacency), 1, sum)
# Mix_Degrees[1, phase]    <- assortativity(A_U, U.Degrees_Labels, directed = TRUE)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Densities
TPA_Density <- data.frame(Phase.Now = "BDA", Net.Density = 1, Type = "PN")
for(phase in 1:(Phase_Num-1)){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  # Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  diag(Adjacency) <- rep(0, N); # diag(Weights) <- rep(0, N)
  Adj.M_P   <- (abs(Adjacency) + Adjacency)/2
  Adj.M_N   <- (abs(Adjacency) - Adjacency)/2
  
  PN_Densities <- c(sum(Adj.M_P)/(N^2-N), sum(Adj.M_N)/(N^2-N))
  TPA_Density.New <- data.frame(Phase.Now = Phs.Func(phase, Rep.Num = 2),
                                Net.Density = PN_Densities, Type = c("P","N"))
  TPA_Density <- rbind(TPA_Density, TPA_Density.New)
}
TPA_Density$Type <- factor(TPA_Density$Type, levels = c("N","P"))
BFig_TPA_Density <- ggplot(data = TPA_Density[-1,],
                           mapping = aes(x=Phase.Now, y=Net.Density, fill=Type))
Fig_TPA_Density <- (
  BFig_TPA_Density
  + geom_bar(stat = "identity", position = "stack")
  + scale_fill_manual(values = c("blue","red"),
                      labels = c("Negative Links","Positive Links"))
  + geom_text(mapping = aes(label = round(Net.Density,2)), colour = 'white',
              vjust = 1.8, hjust = 0.5, position = position_stack())
  + ylim(0,0.35) + labs(x = "", y = "Densities")
  + theme(legend.position=c(0.85,0.85), legend.title=element_blank(),
          legend.background = element_blank())
  )
Fig_TPA_Density





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





Before_Weight <- read.delim("clipboard", header = FALSE)
Before_Weight <- as.matrix(Before_Weight)
During_Weight <- read.delim("clipboard", header = FALSE)
During_Weight <- as.matrix(During_Weight)
After_Weight  <- read.delim("clipboard", header = FALSE)
After_Weight  <- as.matrix(After_Weight)


G7 <- c(5, 15,16,20,27, 35,36)
Rest.Asia    <- c(1:4,6:11)
Rest.Europe  <- c(12:14,17:19,21:26,28:30)
Rest.America <- 31:34
View(Before_Weight[G7,G7])
View(During_Weight[G7,G7])
View(After_Weight[G7,G7])


Interactions_G7.Continents <- function(weight, Signs = c("P","N")){
  G7.Others <- matrix(0, nrow = 4, ncol = 4)
  if(Signs == "P"){ Weight.M <- (abs(weight) + weight)/2 }
  if(Signs == "N"){ Weight.M <- (abs(weight) - weight)/2 }
  if(Signs == "Offset"){ Weight.M <- weight }
  
  G7.Others[1,1] <- sum(Weight.M[G7, G7])
  G7.Others[1,2] <- sum(Weight.M[G7, Rest.Asia])
  G7.Others[1,3] <- sum(Weight.M[G7, Rest.Europe])
  G7.Others[1,4] <- sum(Weight.M[G7, Rest.America])
  
  G7.Others[2,1] <- sum(Weight.M[Rest.Asia, G7])
  G7.Others[2,2] <- sum(Weight.M[Rest.Asia, Rest.Asia])
  G7.Others[2,3] <- sum(Weight.M[Rest.Asia, Rest.Europe])
  G7.Others[2,4] <- sum(Weight.M[Rest.Asia, Rest.America])
  
  G7.Others[3,1] <- sum(Weight.M[Rest.Europe, G7])
  G7.Others[3,2] <- sum(Weight.M[Rest.Europe, Rest.Asia])
  G7.Others[3,3] <- sum(Weight.M[Rest.Europe, Rest.Europe])
  G7.Others[3,4] <- sum(Weight.M[Rest.Europe, Rest.America])
  
  G7.Others[4,1] <- sum(Weight.M[Rest.America, G7])
  G7.Others[4,2] <- sum(Weight.M[Rest.America, Rest.Asia])
  G7.Others[4,3] <- sum(Weight.M[Rest.America, Rest.Europe])
  G7.Others[4,4] <- sum(Weight.M[Rest.America, Rest.America])
  
  return(G7.Others)
}

View(Interactions_G7.Continents(Before_Weight, Signs = "P"))
View(Interactions_G7.Continents(Before_Weight, Signs = "N"))

View(Interactions_G7.Continents(During_Weight, Signs = "P"))
View(Interactions_G7.Continents(During_Weight, Signs = "N"))

View(Interactions_G7.Continents(After_Weight, Signs = "P"))
View(Interactions_G7.Continents(After_Weight, Signs = "N"))









##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##







# Impulse Response Functions
Acc.IR_Unsign <- matrix(0, nrow = N, ncol = Phase_Num*2)
Acc.IR_Offset <- matrix(0, nrow = N, ncol = Phase_Num*2)
# CT_Acc.IR <- matrix(0, nrow = Phase_Num*3, ncol = 3*2)
H <- 30
Imp.Res    <- matrix(0, nrow = N*H, ncol = Phase_Num*N)
Ulti.Influ <- matrix(0, nrow = Phase_Num*N, ncol = N)
for(phase in 1:Phase_Num){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  
  BB <- t(abs(Adjacency) * Weights)  # 转置之后是系数矩阵
  BB_11 <- BB[1:N_Asia, 1:N_Asia]
  BB_12 <- BB[1:N_Asia, (N_Asia+1):(N_Asia+N_Europe)]
  BB_13 <- BB[1:N_Asia, (N_Asia+N_Europe+1):N]
  BB_21 <- BB[(N_Asia+1):(N_Asia+N_Europe), 1:N_Asia]
  BB_22 <- BB[(N_Asia+1):(N_Asia+N_Europe), (N_Asia+1):(N_Asia+N_Europe)]
  BB_23 <- BB[(N_Asia+1):(N_Asia+N_Europe), (N_Asia+N_Europe+1):N]
  BB_31 <- BB[(N_Asia+N_Europe+1):N, 1:N_Asia]
  BB_32 <- BB[(N_Asia+N_Europe+1):N, (N_Asia+1):(N_Asia+N_Europe)]
  BB_33 <- BB[(N_Asia+N_Europe+1):N, (N_Asia+N_Europe+1):N]
  
  Ulti.Influ_phase <- matrix(0, nrow = N, ncol = N)
  for(k in 1:N){
    Imp.Res_k <- matrix(0, nrow = H, ncol = N)
    X_1 <- Imp.Res_k[1, 1:N_Asia]
    Y_1 <- Imp.Res_k[1, (N_Asia+1):(N_Asia+N_Europe)]
    Z_1 <- Imp.Res_k[1, (N_Asia+N_Europe+1):N]
    if(k >= 1 & k <= N_Asia){
      X_1[k] <- 1;  Y_1 <- BB_21 %*% X_1;  Z_1 <- BB_31 %*% X_1 + BB_32 %*% Y_1
    }
    if(k >= N_Asia+1 & k <= N_Asia+N_Europe){ Y_1[k-N_Asia] <- 1;  Z_1 <- BB_32 %*% Y_1  }
    if(k >= N_Asia+N_Europe+1 & k <= N){ Z_1[k-N_Asia-N_Europe] <- 1 }
    Imp.Res_k[1, ] <- c(X_1, Y_1, Z_1)
    
    for(h in 2:H){
      X_t_1 <- Imp.Res_k[h-1, 1:N_Asia]
      Y_t_1 <- Imp.Res_k[h-1, (N_Asia+1):(N_Asia+N_Europe)]
      Z_t_1 <- Imp.Res_k[h-1, (N_Asia+N_Europe+1):N]
      
      X_t <- BB_11 %*% X_t_1 + BB_12 %*% Y_t_1 + BB_13 %*% Z_t_1
      Y_t <- BB_21 %*% X_t   + BB_22 %*% Y_t_1 + BB_23 %*% Z_t_1
      Z_t <- BB_31 %*% X_t   + BB_32 %*% Y_t   + BB_33 %*% Z_t_1
      
      Imp.Res_k[h, ] <- c(X_t, Y_t, Z_t)
    }
    Imp.Res[(k*H-H+1):(k*H), (phase*N-N+1):(phase*N)] <- Imp.Res_k
    Ulti.Influ_phase[, k] <- apply(Imp.Res_k, 2, sum)
  }
  Ulti.Influ[(phase*N-N+1):(phase*N),] <- Ulti.Influ_phase
  
  print(diag(Ulti.Influ_phase))
  diag(Ulti.Influ_phase) <- rep(0, N)
  # Acc.IR_U  <- apply(abs(Ulti.Influ_phase), 2, sum)
  # Acc.IR_Unsign[,(2*phase-1):(2*phase)] <- Ranks_Output(Acc.IR_U, N)
  Acc.IR_O <- apply(Ulti.Influ_phase, 2, sum) # BB.infty - diag(1,N)
  Acc.IR_Offset[,(2*phase-1):(2*phase)] <- Ranks_Output(Acc.IR_O, N)
}
# View(Acc.IR_Unsign)
View(Acc.IR_Offset)
View(Ulti.Influ)   


for(phase in 1:Phase_Num){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  BB        <- t(abs(Adjacency) * Weights)  # 转置之后是系数矩阵
  BB.infty  <- solve( diag(1,N) - BB )
  
  # Offset
  Acc.IR_O  <- apply(BB.infty, 2, sum) # BB.infty - diag(1,N)
  Acc.IR_Offset[,(2*phase-1):(2*phase)] <- Ranks_Output(Acc.IR_O, 30)

  # Analysis on Continents
  CT_Acc.IR_U <- Continents_Sum(BB.infty, type = "U")
  CT_Acc.IR_P <- Continents_Sum(BB.infty, type = "P")
  CT_Acc.IR_N <- Continents_Sum(BB.infty, type = "N")
  CT_Acc.IR[(3*phase-2):(3*phase), ] <- cbind(CT_Acc.IR_U, CT_Acc.IR_P-CT_Acc.IR_N)
}
View(CT_Acc.IR)






StoragePath_Unsign <- paste(StoragePath, "AR Opitmal Lambdas-Unsign.csv")
write.csv(Degree_Analysis_Unsign, StoragePath_Unsign, quote = FALSE)
StoragePath_Offset <- paste(StoragePath, "AR Opitmal Lambdas-Offset.csv")
write.csv(Degree_Analysis_Offset,  StoragePath_Offset,  quote = FALSE)
StoragePath_CT     <- paste(StoragePath, "AR Opitmal Lambdas-Continent.csv")
write.csv(CT.Adjacency,  StoragePath_CT,  quote = FALSE)



# mmP <- (abs(Adjacency_Period) + Adjacency_Period)/2
# mmN <- (abs(Adjacency_Period) - Adjacency_Period)/2
# InDegree <- apply(mmP, 2, sum) - apply(mmN, 2, sum)
# InD   <- Ranks_Output(InDegree, 30)
# OutDegree <- apply(mmP, 1, sum) - apply(mmN, 1, sum)
# OutD   <- Ranks_Output(OutDegree, 30)
# TotalDegree <- InDegree + OutDegree
# TotalD   <- Ranks_Output(TotalDegree, 30)
# Degree_Analysis[, (6*phase-5):(6*phase)] <- cbind(InD, OutD, TotalD)


# Adjusted Degree Analysis
# InDegree.adj <- (apply(abs(adjacency[Asia_area, ]), 2, sum)/N_Asia
#                  + apply(abs(adjacency[Europe_area, ]), 2, sum)/N_Europe 
#                  + apply(abs(adjacency[America_area ,]),2, sum)/N_America)
# InD.adj   <- Ranks_Output(InDegree.adj, 30)
# # OutDegree <- apply(Adjacency, 1, sum)
# OutDegree.adj <- (apply(abs(adjacency[, Asia_area]), 1, sum)/N_Asia
#                   + apply(abs(adjacency[, Europe_area]),  1, sum)/N_Europe
#                   + apply(abs(adjacency[, America_area]), 1, sum)/N_America)
# OutD.adj   <- Ranks_Output(OutDegree.adj, 30)
# TotalDegree.adj <- InDegree.adj + OutDegree.adj
# TotalD.adj   <- Ranks_Output(TotalDegree.adj, 30)
# Degree_Analysis.adjust[, (6*phase-5):(6*phase)] <- cbind(InD.adj, OutD.adj, TotalD.adj)


Acc.IR_Unsign <- matrix(0, nrow = 30, ncol = Phase_Num*2)
Acc.IR_Offset <- matrix(0, nrow = 30, ncol = Phase_Num*2)
CT_Acc.IR <- matrix(0, nrow = Phase_Num*3, ncol = 3*2)
for(phase in 1:Phase_Num){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  BB        <- t(abs(Adjacency) * Weights)  # 转置之后是系数矩阵
  BB.infty  <- solve( diag(1,N) - BB )
  
  # Unsigs
  Acc.IR_U  <- apply(abs(BB.infty), 2, sum) # abs(BB.infty - diag(1,N))
  Acc.IR_Unsign[,(2*phase-1):(2*phase)] <- Ranks_Output(Acc.IR_U, 30)
  
  # Offset
  Acc.IR_O  <- apply(BB.infty, 2, sum) # BB.infty - diag(1,N)
  Acc.IR_Offset[,(2*phase-1):(2*phase)] <- Ranks_Output(Acc.IR_O, 30)
  
  # Analysis on Continents
  CT_Acc.IR_U <- Continents_Sum(BB.infty, type = "U")
  CT_Acc.IR_P <- Continents_Sum(BB.infty, type = "P")
  CT_Acc.IR_N <- Continents_Sum(BB.infty, type = "N")
  CT_Acc.IR[(3*phase-2):(3*phase), ] <- cbind(CT_Acc.IR_U, CT_Acc.IR_P-CT_Acc.IR_N)
}
View(Acc.IR_Unsign)
View(Acc.IR_Offset)
View(CT_Acc.IR)




# Strength Analysis
Strength_Unsign <- matrix(0, nrow = 30, ncol = Phase_Num*6)
Strength_Offset <- matrix(0, nrow = 30, ncol = Phase_Num*6)
Strength_Sumry  <- matrix(0, nrow = 36, ncol = Phase_Num*4)
# Strength_Sumry[, 1] <- Countries
CT.Strength <- matrix(0, nrow = Phase_Num*3, ncol = 3*2)
Totl.Resonance <- matrix(0, nrow = 36, ncol = Phase_Num*2)
Incs.Resonance <- matrix(0, nrow = 36, ncol = Phase_Num*2)
Decs.Resonance <- matrix(0, nrow = 36, ncol = Phase_Num*2)
for(phase in 1:Phase_Num){
  Adjacency <- Final.Adjacency[(phase*N-N+1):(phase*N), ]
  Weights   <- Final.Weight[(phase*N-N+1):(phase*N), ]
  Strength  <- abs(Adjacency) * Weights
  # diag(Strength) <- rep(0, N)
  Adj.M_P   <- (abs(Adjacency) + Adjacency)/2; diag(Adj.M_P) <- rep(0, N)
  Adj.M_N   <- (abs(Adjacency) - Adjacency)/2; diag(Adj.M_N) <- rep(0, N)
  
  # Unsigs
  InS_U  <- apply(abs(Strength), 2, sum)
  OutS_U <- apply(abs(Strength), 1, sum)
  InS_Unsign  <- Ranks_Output(InS_U, 30)
  OutS_Unsign <- Ranks_Output(OutS_U, 30)
  TotalS_Unsign <- Ranks_Output(InS_U+OutS_U, 30)
  Strength_Unsign[,(6*phase-5):(6*phase)] <- cbind(InS_Unsign,OutS_Unsign,TotalS_Unsign)
  
  # Offset
  InS_O  <- apply(Strength, 2, sum)
  OutS_O <- apply(Strength, 1, sum)
  InS_Offset    <- Ranks_Output(InS_O, 30)
  OutS_Offset   <- Ranks_Output(OutS_O, 30)
  TotalS_Offset <- Ranks_Output(InS_O+OutS_O, 30)
  Strength_Offset[,(6*phase-5):(6*phase)] <- cbind(InS_Offset,OutS_Offset,TotalS_Offset)
  
  # Analysis on Continents
  CT.S_U <- Continents_Sum(Strength, type = "U")
  CT.S_P <- Continents_Sum(Strength, type = "P")
  CT.S_N <- Continents_Sum(Strength, type = "N")
  CT.Strength[(3*phase-2):(3*phase), ] <- cbind(CT.S_U, CT.S_P-CT.S_N)
  
  # Resonance Effects
  OutIn_P <- apply(Adj.M_P*Weights, 1, sum) - apply(Adj.M_P*Weights, 2, sum)
  OutIn_N <- abs(apply(Adj.M_N*Weights, 1, sum)) - abs(apply(Adj.M_N*Weights, 2, sum))
  Totl.Resonance[, (2*phase-1):(2*phase)] <- Ranks_Output(OutIn_P-OutIn_N, 36)
  Incs.Resonance[, (2*phase-1):(2*phase)] <- Ranks_Output(OutIn_P, 36)
  Decs.Resonance[, (2*phase-1):(2*phase)] <- Ranks_Output(OutIn_N, 36)
  
  # Strength Summary
  Strength_Sumry[, 4*phase-3] <- apply(abs(Adj.M_P*Weights), 2, sum) # InS+
  Strength_Sumry[, 4*phase-2] <- apply(abs(Adj.M_N*Weights), 2, sum) # InS-
  Strength_Sumry[, 4*phase-1] <- apply(abs(Adj.M_P*Weights), 1, sum) # OutS+
  Strength_Sumry[, 4*phase  ] <- apply(abs(Adj.M_N*Weights), 1, sum) # OutS-
}

View(Strength_Unsign)
View(Strength_Offset)
View(Strength_Sumry)
View(CT.Strength)
View(Totl.Resonance); View(Incs.Resonance); View(Decs.Resonance)

Strength_Unsign_Path <- paste(StoragePath, "Opitmal AR Strength_Unsign.csv")
write.csv(Strength_Unsign, Strength_Unsign_Path, quote = FALSE)
Strength_Offset_Path <- paste(StoragePath, "Opitmal AR Strength_Offset.csv")
write.csv(Strength_Offset, Strength_Offset_Path,  quote = FALSE)
StoragePath_CT       <- paste(StoragePath, "AR Opitmal Lambdas-Continent.csv") 
write.csv(CT.Adjacency,  StoragePath_CT,  quote = FALSE)




##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Classic VAR
{
# 次贷前->次贷中->次贷后->欧债
Period_Range <-c(38930, 39294, 39295, 39903, 39904, 40147, 40148, 40543)
# 次贷前->次贷中1->次贷中2->次贷后->欧债
#Period_Range <-c(38930, 39294, 39295, 39706, 39707, 39903, 39904, 40147, 40148, 40543)

Phase_Num <- length(Period_Range)/2
nlinks_Unsign <- rep(0, Phase_Num)
nlinks_counter <- rep(0, Phase_Num)
Degree_Analysis <- matrix(0, nrow = 30, ncol = Phase_Num*6)
AR_Terms <- matrix(0, nrow = length(Countries), ncol = Phase_Num+1)
AR_Terms[, 1] <- Countries


library(glmnet)
for(phase in 1:Phase_Num){
  
  # Define Variables
  {
    Phase_Start  <- which(Weekdays_v == Period_Range[2*phase-1])
    Phase_End    <- which(Weekdays_v == Period_Range[2*phase])
    XX <- LogReturns[Phase_Start:Phase_End, ]
    N <- ncol(XX)  #Define the dimension of the weight matrix
    Period <- nrow(XX)
    MaxLinks <- N*(N-1)
    
    # Centralized and Standardized
    X <- scale(XX, center = TRUE, scale = TRUE) / sqrt(Period-1)
    #diag(t(X) %*% X)
  }
  
  
  
  Weight <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    y <- X[2:Period, i]
    x <- X[1:(Period-1), ]
    
    auto_term_loc <- rep(1, N); auto_term_loc[i] <- 0
    cvfit <- cv.glmnet(sqrt(Period-1)*x, sqrt(Period-1)*y, nfolds = 10, 
                       type.measure = "mse", alpha = 1, nlambda = 300,
                       intercept = FALSE, standardize = FALSE, 
                       penalty.factor = auto_term_loc)
    
    best_est <- as.vector( coef(cvfit, s="lambda.min") )[-1]
    Weight[i, ] <- best_est
  }
  AR_Terms[, phase+1] <- diag(Weight)
  
  # Generate the adjacency matrix
  adjacency <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      if(i == j) adjacency[i, j] <- 0
      else{
        w <- Weight[i, j]
        if(w > 0)  adjacency[j, i] <-  1
        if(w == 0) adjacency[j, i] <-  0
        if(w < 0)  adjacency[j, i] <- -1
      }
    }
  }
  nlinks_Unsign[phase] <- sum( abs(adjacency) )
  nlinks_counter[phase] <- sum( adjacency )
  
  
  
  # Degree Analysis
  nlinks_Unsign[phase] <- sum( abs(adjacency) )
  nlinks_counter[phase] <- sum( adjacency )
  
  mmP <- (abs(adjacency) + adjacency)/2
  mmN <- (abs(adjacency) - adjacency)/2
  
  InDegree <- apply(mmP, 2, sum) - apply(mmN, 2, sum)
  InD   <- Ranks_Output(InDegree, 30)
  OutDegree <- apply(mmP, 1, sum) - apply(mmN, 1, sum)
  OutD   <- Ranks_Output(OutDegree, 30)
  TotalDegree <- InDegree + OutDegree
  TotalD   <- Ranks_Output(TotalDegree, 30)
  
  Degree_Analysis[, (6*phase-5):(6*phase)] <- cbind(InD, OutD, TotalD)
  
}

nlinks_Unsign; #nlinks_counter
View(AR_Terms)
View(Degree_Analysis)


write.csv(Degree_Analysis,
          "C:/Users/581/Documents/2008 Subprime Crisis by AR Opitmal Lambdas.csv",
          quote = FALSE)
}

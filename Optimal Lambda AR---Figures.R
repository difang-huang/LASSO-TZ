##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##     Figures in Time-Zone VAR Model by Optimal Lambda AR                             ##
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
StoragePath = "C:/Users/581/Documents/Optimal LASSO Networks/"


library(glmnet)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggthemes)





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Figure: Network Densities over Five Periods                                ##
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


# Calculations
NewFig_Density = matrix(0, nrow = Phase.Num, ncol = 3)
TPA_Density = data.frame(Phase.Now = "BDA", Net.Density = 1, Type = "PN")
for(phase in 1:(Phase.Num)){

# Input DATA
Adjacency = FULL.Adj.M[(phase*N-N+1):(phase*N), ]
# Weights = FULL.Wet.M[(phase*N-N+1):(phase*N), ]
diag(Adjacency) = rep(0, N); # diag(Weights) = rep(0, N)
Adj.M_P   = (abs(Adjacency) + Adjacency)/2
Adj.M_N   = (abs(Adjacency) - Adjacency)/2

# Calculate unsigned, positive and negative densities
NewFig_Density[phase, ] = c(sum(abs(Adjacency))/(N^2-N), sum(Adj.M_P)/(N^2-N), sum(Adj.M_N)/(N^2-N))

# Output ggplot-type data
TPA_Density.New = data.frame(Phase.Now = Phs.Func(phase, Rep.Num = 3),
                             Net.Density = NewFig_Density[phase, ], Type = c("U","P","N"))
TPA_Density = rbind(TPA_Density, TPA_Density.New)

}
write.csv(NewFig_Density, paste(StoragePath, "NewFigure-Network Density.csv", sep = ""), 
          quote = FALSE)


# Draw ggplot figure of Network Density
TPA_Density$Type = factor(TPA_Density$Type, levels = c("N","P","U"))
TPA_Density$Phase.Now = factor(TPA_Density$Phase.Now,
   levels = c("Before 2008 Subprime Crisis", "During 2008 Subprime Crisis",
              "After 2008 Subprime Crisis", "During European Debt Crisis",
              "After European Debt Crisis")
)
BFig_TPA_Density = ggplot(data = TPA_Density[-1,],
                          mapping = aes(x=Phase.Now, y=Net.Density, fill=Type))
Fig_TPA_Density = (
  BFig_TPA_Density
  + geom_bar(stat = "identity", position = "dodge", alpha=0.8) # position = 'stack'
  + scale_fill_manual(values = c("blue","red","grey40"),
                      labels = c("Negative Links","Positive Links","Unsigned Links"))
  + geom_text(mapping = aes(label = round(Net.Density,3)), colour = 'white', vjust = 1.8,
              position = position_dodge(0.9))
  # , hjust = 0.5
  + ylim(0,0.35) + labs(x = "", y = "Network Density")
  + theme_calc()
  + theme(legend.position=c(0.9,0.9), legend.title=element_blank(),
          legend.background = element_blank())
)
Fig_TPA_Density


ggsave(filename = "NewFig_Network_Density.eps", device = cairo_ps, dpi = 300,
       path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data")
ggsave(filename = "NewFig_Network_Density.tiff", device = "tiff", dpi = 300,
       width = 15, height = 10.6, units = "cm",
       path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data")





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Figure: Heat Map of VAR Coefficients                                       ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Define basic functions
Title_HeatMap.VAR.Coefs = function(phs){
  switch(phs,
         "Before Subprime Crisis (1 August 2006 to 31 July 2007)",
         "During Subprime Crisis (1 August 2007 to 31 March 2009)",
         "After Subprime Crisis (1 April 2009 to 30 November 2009)",
         "During European Debt Crisis (1 December 2009 to 16 December 2013)",
         "After European Debt Crisis (17 December 2013 to 31 December 2015)")
}


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Heat Map of VAR Coefficients
for(phase in 1:Phase.Num){

Adj.M = FULL.Adj.M[(phase*N-N+1):(phase*N), ]
Wet.M = FULL.Wet.M[(phase*N-N+1):(phase*N), ]
rownames(Adj.M) = Countries
colnames(Adj.M) = Countries


VAR.Coefs.M = matrix(0, nrow = N*N, ncol = 6)
for(ii in 1:N){for(jj in 1:N){
  loc = (ii-1)*N + jj
  VAR.Coefs.M[loc, 1] = Country.Abbr[ii]
  VAR.Coefs.M[loc, 2] = switch((ii>=1)+(ii>=12)+(ii>=30), "Asia", "Europe", "Americas")
  VAR.Coefs.M[loc, 3] = Country.Abbr[jj]
  VAR.Coefs.M[loc, 4] = switch((jj>=1)+(jj>=12)+(jj>=30), "Asia", "Europe", "Americas")
  VAR.Coefs.M[loc, 5] = Adj.M[ii, jj]
  VAR.Coefs.M[loc, 6] = Wet.M[ii, jj]
}}


DF_AdjWet.M = data.frame(From = VAR.Coefs.M[, 1], From_ConT = VAR.Coefs.M[, 2],
                         To   = VAR.Coefs.M[, 3], To_ConT   = VAR.Coefs.M[, 4],
                         Adj = VAR.Coefs.M[, 5], Wet = VAR.Coefs.M[, 6])
DF_AdjWet.M$From = factor(DF_AdjWet.M$From, levels = rev(Country.Abbr))
DF_AdjWet.M$To   = factor(DF_AdjWet.M$To,   levels = Country.Abbr)
DF_AdjWet.M$From_ConT = factor(DF_AdjWet.M$From_ConT, levels = c("Asia", "Europe", "Americas"))
DF_AdjWet.M$To_ConT   = factor(DF_AdjWet.M$To_ConT,   levels = c("Asia", "Europe", "Americas"))
DF_AdjWet.M$Adj  = factor(DF_AdjWet.M$Adj,  levels = c(1,0,-1))
BFig_AdjWet.M = ggplot(data = DF_AdjWet.M,
                       mapping = aes(x = To, y = From, fill = Adj))
Fig_AdjWet.M = (
  BFig_AdjWet.M + geom_tile(alpha=0.6)   # alpha=透明度
  + scale_fill_manual(name = "Link Sign", breaks = c(1,0,-1),
                      values = c("red","grey85","blue"),
                      labels = c("Positive", "No Exist", "Negative"))
  + facet_grid(rows = vars(From_ConT), cols = vars(To_ConT), scales="free", space = "free")
  # scale, space = ???????ݶԿ̶ȡ???Ԫ??????????
  + labs(x = "To", y = "From", title =Title_HeatMap.VAR.Coefs(phs = phase))
  + theme_calc()   # bw calc pander // stata hc
  + theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black"),
          legend.background = element_blank())        # ȥ??????
  # + theme(panel.border = element_rect(fill"black", size=1, linetype="solid"))
  # + geom_vline(xintercept=c(11.5, 30.5), size=0.8)
  # + geom_hline(yintercept=c(6.5, 25.5), size=0.8)
)
plot(Fig_AdjWet.M)


ggsave(
  filename = switch (phase,
    "NewFig_HeatMap.VAR.Coefs_Subprime.Before.eps", "NewFig_HeatMap.VAR.Coefs_Subprime.During.eps",
    "NewFig_HeatMap.VAR.Coefs_Subprime.After.eps", "NewFig_HeatMap.VAR.Coefs_EuDebt.During.eps",
    "NewFig_HeatMap.VAR.Coefs_EuDebt.After.eps"
  ),
  device = cairo_ps, dpi = 300, width = 25, height = 15, units = "cm",
  path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)


ggsave(
  filename = switch (phase,
    "NewFig_HeatMap.VAR.Coefs_Subprime.Before.png", "NewFig_HeatMap.VAR.Coefs_Subprime.During.png",
    "NewFig_HeatMap.VAR.Coefs_Subprime.After.png", "NewFig_HeatMap.VAR.Coefs_EuDebt.During.png",
    "NewFig_HeatMap.VAR.Coefs_EuDebt.After.png"
  ),
  device = "png", dpi = 300, width = 25, height = 15, units = "cm",
  path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)

}





##=====================================================================================##
##=====================================================================================##
##=====================================================================================##
##          Figure: National Market Statuses Based on Net Strength                     ##
##=====================================================================================##
##=====================================================================================##
##=====================================================================================##





# Basic calculations for ggplot
NewTab_Strength = matrix(0, nrow = N, ncol = 6*Phase.Num)
Net.Strength = data.frame(Agent = "aaa", Continent = "bbb", Net.InS  = 1,
                           Net.OutS = 1, Phase.Now = "ccc")
for(phase in 1:Phase.Num){

# Define basic variables
Adj.M = FULL.Adj.M[(phase*N-N+1):(phase*N), ]
Wet.M = FULL.Wet.M[(phase*N-N+1):(phase*N), ]
diag(Wet.M) = rep(0, N)
Adj.M_P   = (abs(Adj.M) + Adj.M)/2; diag(Adj.M_P) = rep(0, N)
Adj.M_N   = (abs(Adj.M) - Adj.M)/2; diag(Adj.M_N) = rep(0, N)


# Strength summary
InS_P  = apply(abs(Adj.M_P*Wet.M), 2, sum) 
InS_N  = apply(abs(Adj.M_N*Wet.M), 2, sum)
OutS_P = apply(abs(Adj.M_P*Wet.M), 1, sum)
OutS_N = apply(abs(Adj.M_N*Wet.M), 1, sum)
NewTab_Strength[, 6*phase-5] = InS_P;  NewTab_Strength[, 6*phase-4] = InS_N
NewTab_Strength[, 6*phase-3] = OutS_P; NewTab_Strength[, 6*phase-2] = OutS_N
NewTab_Strength[, 6*phase-1] = OutS_P - OutS_N
NewTab_Strength[, 6*phase]   = InS_P  - InS_N


# Output ggplot-type data
Net.New = data.frame(Agent = Country.Abbr,
                     Continent = c(rep("Asia",11), rep("Europe",19), rep("Americas",6)),
                     Net.InS  = InS_P-InS_N, Net.OutS = OutS_P-OutS_N,
                     Phase.Now = Phs.Func(phase))
Net.Strength = rbind(Net.Strength, Net.New)

}
Net.Strength = Net.Strength[-1, ]
Net.Strength$Continent = factor(Net.Strength$Continent, levels = c("Asia", "Europe", "Americas"))
View(Net.Strength)

write.csv(NewTab_Strength,
          paste(StoragePath, "NewFigure-Four Strengths and Net Strengths.csv", sep = ""),
          quote = FALSE)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Draw Scatter Plots of Net Strengths in 2008 Subprime Crisis---3*1 Form
RBG = c("red", "blue", "green")
RBG = setNames(RBG, c("Asia", "Europe", "America"))

Net.Strength_Subprime = Net.Strength[1:(3*N), ]
Net.Strength_Subprime$Phase.Now = c(rep("Before",N), rep("During",N), rep("After",N))
Net.Strength_Subprime$Phase.Now = factor(Net.Strength_Subprime$Phase.Now,
                                         levels = c("Before", "During", "After"))
BFig_Net.S_Subprime = ggplot(
  data = Net.Strength_Subprime,
  mapping = aes(x = Net.OutS, y = Net.InS, color = Continent, shape = Continent)
)
Fig_Net.S_Subprime = (
  BFig_Net.S_Subprime + geom_point(size = 2.5, alpha = 0.5)
  + scale_color_manual(values = RBG) + guides(color = "none")
  + scale_shape_manual(values = c(16,18,17))
  + geom_text(aes(label=rep(Country.Abbr,3)), size=5)
  + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  + xlim(-2,5) + ylim(-0.1,1)
  + labs(title = "Net In- and Out- Strengths in 2008 Subprime Crisis",
         x = "Net Out-Strength", y = "Net In-Strength")
  + theme_calc() + theme(panel.grid=element_blank())
)
Fig_Net.S_Subprime + facet_grid(Phase.Now~.)


ggsave(filename = "NewFig_Net.Strengths_Subprime.eps", device = cairo_ps, dpi = 300,
  path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Draw Scatter Plots of Net Strengths in European Debt Crisis---2*1 Form
RBG = c("red", "blue", "green")
RBG = setNames(RBG, c("Asia", "Europe", "America"))

Net.Strength_EuDebt = Net.Strength[(3*N+1):(5*N), ]
Net.Strength_EuDebt$Phase.Now = c(rep("During",N), rep("After",N))
Net.Strength_EuDebt$Phase.Now = factor(Net.Strength_EuDebt$Phase.Now,
                                       levels = c("During", "After"))
BFig_Net.S_EuDebt = ggplot(
  data = Net.Strength_EuDebt,
  mapping = aes(x = Net.OutS, y = Net.InS, color = Continent, shape = Continent)
)
Fig_Net.S_EuDebt = (
  BFig_Net.S_EuDebt + geom_point(size = 2.5, alpha = 0.5)
  + scale_color_manual(values = RBG) + guides(color = "none")
  + scale_shape_manual(values = c(16,18,17))
  + geom_text(aes(label=rep(Country.Abbr,2)), size=5)
  + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  + xlim(-2,5) + ylim(-0.1, 1)
  + labs(title = "Net In- and Out- Strengths in European Debt Crisis",
         x = "Net Out-Strength", y = "Net In-Strength")
  + theme_calc() + theme(panel.grid=element_blank())
)
Fig_Net.S_EuDebt + facet_grid(Phase.Now~.)


ggsave(filename = "NewFig_Net.Strengths_EuDebt.eps", device = cairo_ps, dpi = 300,
       path = "D:/Research/2019  Financial Risk Networks by LASSO/New Data"
)


##=====================================================================================##
##=====================================================================================##
##=====================================================================================##


# Roles Detection
Roles.Detect = function(Var.X, Var.Y){
  if(Var.X > 0 & Var.Y > 0){ return(1) }
  if(Var.X < 0 & Var.Y > 0){ return(2) }
  if(Var.X < 0 & Var.Y < 0){ return(3) }
  if(Var.X > 0 & Var.Y < 0){ return(4) }
  
  if(Var.X >  0 & Var.Y == 0){ return(14) }
  if(Var.X == 0 & Var.Y >  0){ return(12) }
  if(Var.X <  0 & Var.Y == 0){ return(23) }
  if(Var.X == 0 & Var.Y <  0){ return(34) }
}


# Country Status Detection
NewFig_Country.Roles = matrix(0, nrow = N, ncol = Phase.Num)
for(phase in 1:Phase.Num){
  Net.OutS = NewTab_Strength[, 6*phase-1]
  Net.InS  = NewTab_Strength[, 6*phase]
  
  for(kk in 1:N){
    NewFig_Country.Roles[kk, phase] = Roles.Detect(Net.OutS[kk], Net.InS[kk])
  }
}
View(NewFig_Country.Roles)

write.csv(NewFig_Country.Roles,
          paste(StoragePath, "NewTable-Country Statuses Based on Net Strengths.csv", sep = ""),
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





DF_AdjWet.M = data.frame(From = VAR.Coefs.M[, 1], From_ConT = VAR.Coefs.M[, 2],
                         To   = VAR.Coefs.M[, 3], To_ConT   = VAR.Coefs.M[, 4],
                         Adj = VAR.Coefs.M[, 5], Wet = VAR.Coefs.M[, 6])
DF_AdjWet.M$From = factor(DF_AdjWet.M$From, levels = rev(Country.Abbr))
DF_AdjWet.M$To   = factor(DF_AdjWet.M$To,   levels = Country.Abbr)
DF_AdjWet.M$From_ConT = factor(DF_AdjWet.M$From_ConT, levels = c("America", "Europe", "Asia"))
DF_AdjWet.M$To_ConT   = factor(DF_AdjWet.M$To_ConT,   levels = c("Asia", "Europe", "America"))
DF_AdjWet.M$Adj  = factor(DF_AdjWet.M$Adj,  levels = c(1,0,-1))
BFig_AdjWet.M = ggplot(data = DF_AdjWet.M,
                       mapping = aes(x = To, y = From, fill = Adj))
Fig_AdjWet.M = (
  BFig_AdjWet.M + geom_tile(alpha=0.5)   # alpha=透明度
  + scale_fill_manual(values = c("red","lightgrey","blue"))
  + geom_vline(xintercept=c(11.5, 30.5), size=0.8)
  + geom_hline(yintercept=c(6.5, 25.5), size=0.8)
  + facet_grid(From_ConT~To_ConT)
  # + facet_grid(rows = From_ConT, cols = To_ConT)
  + theme_bw()   # bw stata calc hc pander
  + theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))        # ȥ??????
  # + theme(panel.border = element_rect(fill"black", size=1, linetype="solid"))
)
Fig_AdjWet.M
#2016 hopping performance analyses
# Author: Rory Telemeco

rm(list=ls(all = TRUE))

library(stringr)
library(foreach)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)
library(nls2)
library(foreach)

source('Functions.R')

####### Input files ##########

setwd("./data/")

HopRaw <- read.csv(file = "Pdodge_PerformanceData1_06202016_ALL.csv", header = TRUE, stringsAsFactors = FALSE)

##### Bookkeeping ########
#HopRaw <- HopRaw[is.na(HopRaw$TEST.TEMP1) == FALSE, ] #remove individuals that we couldn't measure

Labels <- str_split_fixed(HopRaw[,1], "_", n = 4) #collect Temp and Light treatments
HopRaw$Temp <- Labels[,2]
HopRaw$Light <- Labels[,3]

#stack TEST1, TEST2, and TEST3 in pseudo-longform (HopPL)
Hop1 <- HopRaw[, c(1:6, 60:61, 7:19)]
colnames(Hop1) <- c("ID", "Batch", "CC", "Site", "Species", "Sex", "Temp", "Light", "TestTemp", "Y0", "X0", "Y1", "X1", "Y2", "X2", "Y3", "X3", "Y4", "X4", "Y5", "X5")
Hop2 <- HopRaw[, c(1:6, 60:61, 20:32)]
colnames(Hop2) <- c("ID", "Batch", "CC", "Site", "Species", "Sex", "Temp", "Light", "TestTemp", "Y0", "X0", "Y1", "X1", "Y2", "X2", "Y3", "X3", "Y4", "X4", "Y5", "X5")
Hop3 <- HopRaw[, c(1:6, 60:61, 33:45)]
colnames(Hop3) <- c("ID", "Batch", "CC", "Site", "Species", "Sex", "Temp", "Light", "TestTemp", "Y0", "X0", "Y1", "X1", "Y2", "X2", "Y3", "X3", "Y4", "X4", "Y5", "X5")
Hop4 <- HopRaw[, c(1:6, 60:61, 46:58)]
colnames(Hop4) <- c("ID", "Batch", "CC", "Site", "Species", "Sex", "Temp", "Light", "TestTemp", "Y0", "X0", "Y1", "X1", "Y2", "X2", "Y3", "X3", "Y4", "X4", "Y5", "X5")

HopPL <- rbind(Hop1, Hop2, Hop3, Hop4)

#remove rows with no values
HopPL <- HopPL[is.na(HopPL$TestTemp) == FALSE, ] #remove individuals that we couldn't measure


###### Convert measures from Feet_Inches to cm #######

#Function to convert a single vector (input is a Feet_Inches vector and output is a CM vector)
HopsToCM <- function(Vec) {
  Div <- data.frame(str_split_fixed(Vec, "_", n = 2), stringsAsFactors = FALSE) #divide components
  DivI <- 12 * as.numeric(Div[,1]) #convert feet to inches
  Inches <- DivI + as.numeric(Div[,2]) #total inches
  CM <- Inches * 2.5 #convert inches to cm
  return(CM)
}

#Function to convert an entire dataframe (uses above function) (input is a DF of Feet_Inches vectors and output is a DF of CM vectors)
HopsToCM_ALL <- function(DF) {
  NewPoints <- foreach(i = 1:ncol(DF), .combine = cbind) %do%
    HopsToCM(DF[,i])
  NewPoints <- data.frame(NewPoints)
  colnames(NewPoints) <- colnames(DF)
  return(NewPoints)
}

#convert all measures to CM using above functions
CM <- HopsToCM_ALL(HopPL[,10:21])

#reattach metadata
HopCM <- cbind(HopPL[,1:9], CM)


##### Calculate Euclidian Distances ######

EucDist <- function (X1, Y1, X2, Y2) {sqrt((X2-X1)^2 + (Y2-Y1)^2)} #Input coordinates and output Euclidean Distance

HopDist <- function(DF) {
  Hop1 <- EucDist(DF$X0, DF$Y0, DF$X1, DF$Y1)
  Hop2 <- EucDist(DF$X1, DF$Y1, DF$X2, DF$Y2)
  Hop3 <- EucDist(DF$X2, DF$Y2, DF$X3, DF$Y3)
  Hop4 <- EucDist(DF$X3, DF$Y3, DF$X4, DF$Y4)
  Hop5 <- EucDist(DF$X4, DF$Y4, DF$X5, DF$Y5)
  Hops <- data.frame(cbind(Hop1, Hop2, Hop3, Hop4, Hop5))
  return(Hops)
} #Input DF with numeric columns labeled "X0" "Y0" ... "X5" "Y5" and get dataframe of hop distances as output.  Uses EucDist above

HopLengths <- HopDist(HopCM)

#Add distances to full dataframe
Hop <- cbind(HopCM, HopLengths)


##### Calculate summary statistics for hops (use "apply" to use function across rows rather than columns) #####
Hop$Mean <- apply(HopLengths, 1, mean, na.rm = TRUE)
Hop$SD <- apply(HopLengths, 1, sd, na.rm = TRUE)
Hop$Min <- apply(HopLengths, 1, min, na.rm = TRUE)
Hop$Max <- apply(HopLengths, 1, max, na.rm = TRUE)
Hop$Median <- apply(HopLengths, 1, median, na.rm = TRUE)


####### Assign Factors ######

Hop$Batch <- as.factor(Hop$Batch)
Hop$Site <- as.factor(Hop$Site)
Hop$Species <- as.factor(Hop$Species)
Hop$Sex <- as.factor(Hop$Sex)
Hop$Temp <- as.factor(Hop$Temp)
Hop$Light <- as.factor(Hop$Light)
Hop$TestTempF <- as.factor(Hop$TestTemp)
Hop$Treat <- paste(Hop$Temp, Hop$Light, sep = "_")

######## assess sample sizes ##########

table(Hop[, c(4,6,9,33)]) #results collated in google drive file "Performance Sample Sizes"


##### Plots of the data #####

#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}
Give.N <- function(x){return(c(y = MaxSE(x)*1.1, label = length(x))) }


#averaging
ggplot(data = Hop, aes(x = TestTemp, y = Mean, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 5)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 5)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = 5), size = 5) +
  ylab("Mean hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))

ggplot(data = Hop, aes(x = TestTemp, y = Mean, col = Site, linetype = Sex)) +
  rory_theme +
  facet_grid(Light~Temp) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 3)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 3)) +
  ylab("Mean hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))


#every individual
ggplot(data = Hop, aes(x = TestTemp, y = Mean, col = Site, group = ID)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  geom_line()+
  ylab("Mean hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))


#plot based on individual medians
ggplot(data = Hop, aes(x = TestTemp, y = Median, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1) +
  stat_summary(fun.y = mean, geom = "line", size = 2) +
  ylab("Median hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))

ggplot(data = Hop, aes(x = TestTemp, y = Median, col = Site, linetype = Sex)) +
  rory_theme +
  facet_grid(Light~Temp) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1) +
  stat_summary(fun.y = mean, geom = "line", size = 2) +
  ylab("Median hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))


#plot based on individual maxima 

ggplot(data = Hop, aes(x = TestTemp, y = Max, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1) +
  stat_summary(fun.y = mean, geom = "line", size = 2) +
  ylab("Max hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))

ggplot(data = Hop, aes(x = TestTemp, y = Max, col = Site, linetype = Sex)) +
  rory_theme +
  facet_grid(Light~Temp) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1) +
  stat_summary(fun.y = mean, geom = "line", size = 2) +
  ylab("Maximum hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))


#every individual
ggplot(data = Hop, aes(x = TestTemp, y = Max, col = Site, group = ID)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  geom_line()+
  ylab("Max hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))

#looking at median instead of mean does not make qualitative changes
#looking at the max makes some difference



##### Analyze performance with linear mixed models #####


Hop_2 <- merge(Hop, Instars[,c(1,3,54)], by.x = "ID", by.y = "Individual") #Instars from the "GrowthRate_Script_2016"
hist(Hop_2$Mean) #looks pretty normal

poly(x, 3, raw=TRUE)

LmerH0 <- lmer(Mean ~ poly(TestTemp, 2, raw = TRUE) * Sex * Site * Temp * Light + 
                 (1|Female2), 
               #REML = FALSE,
               na.action = 'na.omit', data = Hop_2)

LmH0 <- lm(Mean ~ poly(TestTemp, 2, raw = TRUE) * Sex * Site * Temp * Light, data = Hop_2)

step(LmH0)

#preferred model from step:
LmH1 <- lm(Mean ~ poly(TestTemp, 2, raw = TRUE) + Sex + Site + Temp + Light + 
             poly(TestTemp, 2, raw = TRUE):Sex + poly(TestTemp, 2, raw = TRUE):Temp + Sex:Temp + Site:Temp + poly(TestTemp, 2, raw = TRUE):Light + Sex:Light + Site:Light + Temp:Light + 
             poly(TestTemp, 2, raw = TRUE):Sex:Temp + poly(TestTemp, 2, raw = TRUE):Sex:Light + poly(TestTemp, 2, raw = TRUE):Temp:Light + Sex:Temp:Light + Site:Temp:Light + 
             poly(TestTemp, 2, raw = TRUE):Sex:Temp:Light, 
           data = Hop_2)

LmerH1 <- lmer(Mean ~ poly(TestTemp, 2, raw = TRUE) + Sex + Site + Temp + Light + 
                 poly(TestTemp, 2, raw = TRUE):Sex + poly(TestTemp, 2, raw = TRUE):Temp + Sex:Temp + Site:Temp + poly(TestTemp, 2, raw = TRUE):Light + Sex:Light + Site:Light + Temp:Light + 
                 poly(TestTemp, 2, raw = TRUE):Sex:Temp + poly(TestTemp, 2, raw = TRUE):Sex:Light + poly(TestTemp, 2, raw = TRUE):Temp:Light + Sex:Temp:Light + Site:Temp:Light + 
                 poly(TestTemp, 2, raw = TRUE):Sex:Temp:Light +
                 (1|Female2), 
               #REML = FALSE,
               na.action = 'na.omit', data = Hop_2)

AIC(LmH1)
AIC(LmerH1) #preferred and looks good after graphical validation

Anova(LmerH1, type = 3) #The 4-way interaction between the polynomial and sex, temp, and light is significant.  Everything is having effects, pretty much

lsmeans(LmerH1, pairwise ~ Site|Temp+Sex+Light, adjust = "none") #A1 and C1 differ from one another in HV_S, LV_L, and LV_S
lsmeans(LmerH1, pairwise ~ Light|Temp+Sex, adjust = "none") #Light generally has an effect in LV but not HV
lsmeans(LmerH1, pairwise ~ Temp|Light+Sex) #Temp regime had an effect under every light and site, with LV always greater than HV

#plot the LSmeans
ObjH<-data.frame(summary(lsmeans(LmerH1, ~ TestTemp*Sex*Site*Temp*Light, at = list(TestTemp = c(25,17,10,35)), type = "response")))
Obj0<-data.frame(summary(lsmeans(LmerH0, ~ TestTemp*Sex*Site*Temp*Light, at = list(TestTemp = c(25,17,10,35)), type = "response")))


ggplot(data = ObjH, aes(y = lsmean, x = TestTemp, col = Site, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), position = position_dodge(width = 3)) +
  geom_line(aes(x = TestTemp), size = 1.5, position = position_dodge(width = 3)) +
  geom_point(size = 6, position = position_dodge(width = 3)) +
  ylab("Mean hopping distance\n(cm, LS-means)") + xlab("Site") +
  scale_x_continuous(breaks = c(10, 17, 25, 35)) +
  scale_shape_manual(values = c(19, 1))

##### What if analyze all measures, rather than the means #####

#make data long form
Hop_Long <- melt(Hop_2[, c(1:9, 22:26, 33:35)], id.vars = c("ID", "Batch", "CC", "Site", "Species", "Sex", "Temp", "Light", "Treat", "Female", "Female2", "TestTemp"))

#full model (huge)

LmerL_H0 <- lmer(value ~ poly(TestTemp, 2, raw = TRUE) * Sex * Site * Temp * Light + 
                   (1|variable/ID/Female2), 
                 #REML = FALSE,
                 na.action = 'na.omit', data = Hop_Long)
#the 5-way interaction is significant, so should work from the full model

ObjL0<-data.frame(summary(lsmeans(LmerL_H0, ~ TestTemp*Sex*Site*Temp*Light, at = list(TestTemp = c(25,17,10,35)), type = "response")))


ggplot(data = ObjL0, aes(y = lsmean, x = TestTemp, col = Site, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), position = position_dodge(width = 3)) +
  geom_line(aes(x = TestTemp), size = 1.5, position = position_dodge(width = 3)) +
  geom_point(size = 6, position = position_dodge(width = 3)) +
  ylab("Mean hopping distance\n(cm, LS-means)") + xlab("Site") +
  scale_x_continuous(breaks = c(10, 17, 25, 35)) +
  scale_shape_manual(values = c(19, 1))

##### Compare hopping and optimal performance (COARSE) ######

#Function used to collect row with maximum mean hopping distance
MaxHopRow <- function (DF) {DF[DF$Mean == max(DF$Mean),]}

#Use foreach to collect a data row for each individual with temp of highest perforance
Opt_DF <- foreach(i = unique(Hop_2$ID), .combine = rbind) %do%
  MaxHopRow(Hop_2[Hop_2$ID == i, c(1,4,5,6,7,8,9,27)])
colnames(Opt_DF) <- c("ID", "Site", "Species", "Sex", "Temp", "Light", "To_Hop", "PTo_Hop" )

#merge with PBT data (from other script) to have PBT and T with optimal performance
Opt_DF <- merge(Opt_DF, Obj_PBT_Ind[, c(1:3)], by.x = "ID", by.y = "Individual")

plot(Opt_DF$PBT, Opt_DF$To_Hop) #no obvious correlation
plot(Opt_DF$PBT, Opt_DF$PTo_Hop) #no obvious correlation
plot(Opt_DF$To_Hop, Opt_DF$PTo_Hop) #no obvious correlation

HopCor1 <- lm(To_Hop ~ PBT*Sex*Temp*Light*Site, data = Opt_DF)
step(HopCor1)

HopCor1.1 <- lm(To_Hop ~ PBT + Sex + Temp + Light + Site + PBT:Sex + 
                  PBT:Temp + Sex:Temp + PBT:Light + Sex:Light + Temp:Light + 
                  PBT:Site + Sex:Site + Light:Site + PBT:Sex:Temp + PBT:Sex:Light + 
                  PBT:Temp:Light + Sex:Temp:Light + Sex:Light:Site + PBT:Sex:Temp:Light, data = Opt_DF)
Anova(HopCor1.1, type = 3) #potentially an interaction between PBT and Sex on the optimal hopping temperature

HopCor1.2 <- lm(To_Hop ~ PBT * Sex, data = Opt_DF)
Anova(HopCor1.2, type = 3) #not significant outside of the more complex model though


##### Estimate thermal performance curves #####

#write formular for performance curve from Buckley and Nufio 2014 (product of Gaussian and Gompertz functions)
Performance <- function (Te, Pmax, To, rho, sig) {
  Pmax * (exp(1)^(-exp((rho*(Te-To))-6) - (sig*((Te-To)^2))))
}
#Te = temperature at which performance is assessed
#Pmax = performance at optimal temperature (in units of performance)
#To = optimal temperature (C)
#rho  = sensitivity of performance above To
#sig = sensitivity of performance below To


#Loop to build a dataframe with TPC parameters for each individual based on hopping performance, uses PBT information from that script
#Uses nls2 to calculate the performance curves

TPC_param <- data.frame(foreach(i = unique(Opt_DF$ID), .combine = rbind) %do% {
  #PBT <- Opt_DF[Opt_DF$ID == i , "PBT"]
  PBT <- 40
  Start <- data.frame(Pmax = c(10,70), To = c(PBT-15,PBT+15), sig = c(0, 0.7))
  NLS <- nls2(y ~ Performance(Te, Pmax, To, rho = 0.7, sig), start = Start, data = TPC_data[TPC_data$ID == i ,], algorithm = "brute-force") 
  c(i, as.vector(coef(NLS)), 0.7) })

colnames(TPC_param) <- c("ID", "Pmax", "To", "sig", "rho")

#Add To calculated from TPC's to hopping data to look for correlations
Opt_DF2 <- merge(Opt_DF, TPC_param[,c(1,3)], by = "ID")

# see if anything jumps out
ggplot(data = Opt_DF2, aes(x = PBT, y = as.numeric(as.character(To)), col = Site)) +
  rory_theme +
  geom_point(size = 6)

#I'm unsure what is causing the big breaks...

HopCor2 <- lm(as.numeric(as.character(To)) ~ PBT*Sex*Temp*Light*Site, data = Opt_DF2)
step(HopCor2) #nearly the full model...

HopCor2.1 <- lm(as.numeric(as.character(To)) ~ PBT*(Sex+Temp*Light+Site), data = Opt_DF2)
Anova(HopCor2.1) #Light and PBT may be important for affecting To, but that's about it


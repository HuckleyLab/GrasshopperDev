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
library(lme4)
library(MuMIn)

source('Functions.R')

####### Input files ##########

setwd("./data/")

FeedRaw <- read.csv(file = "FeedingRate_Data_All_03152017.csv", header = TRUE, stringsAsFactors = FALSE)

##### Bookkeeping ########

#bring in hopping data which have sex
HopRaw <- read.csv(file = "Pdodge_PerformanceData1_06202016_ALL.csv", header = TRUE, stringsAsFactors = FALSE) #remove individuals that we couldn't measure
#add sex with merge
FeedRaw <- merge(FeedRaw, HopRaw[,c(1,6)], by.x = "Individual.ID..EGGPOD..HATCHING._TEMP_LIGHT_SPECIES_SITE.", by.y = "Individual", all.x = TRUE)
#divide super label into informative treatments
Labels <- str_split_fixed(FeedRaw[,1], "_", n = 4) #collect Temp and Light treatments
FeedRaw$Temp <- Labels[,2]
FeedRaw$Light <- Labels[,3]
FeedRaw$Population <- Labels[,4]
FeedRaw$ID <- Labels[,1]
ID_Labels <- str_split_fixed(FeedRaw$ID, "[.]", n = 2)
FeedRaw$Mom <- ID_Labels[,1]
FeedRaw$Species <- "dodg"

#stack Feed1, Feed2, and Feed3 in pseudo-longform (FeedPL)
Feed1 <- FeedRaw[, c(31:33, 27:30, 1:2, 3:10)]
colnames(Feed1) <- c("ID", "Mom", "Species", "Sex", "Temp", "Light", "Site", "FULL_ID", "Batch", "Date", "TestTemp", "SubBatch", "Position", "Area1", "Area2", "Area3", "Area4")
Feed2 <- FeedRaw[, c(31:33, 27:30, 1:2, 11:18)]
colnames(Feed2) <- c("ID", "Mom", "Species", "Sex", "Temp", "Light", "Site", "FULL_ID", "Batch", "Date", "TestTemp", "SubBatch", "Position", "Area1", "Area2", "Area3", "Area4")
Feed3 <- FeedRaw[, c(31:33, 27:30, 1:2, 19:26)]
colnames(Feed3) <- c("ID", "Mom", "Species", "Sex", "Temp", "Light", "Site", "FULL_ID", "Batch", "Date", "TestTemp", "SubBatch", "Position", "Area1", "Area2", "Area3", "Area4")

FeedPL <- rbind(Feed1, Feed2, Feed3)

#remove rows with no values or missing sex
FeedPL <- FeedPL[is.na(FeedPL$TestTemp) == FALSE, ] #remove individuals that we couldn't measure
FeedPL <- FeedPL[is.na(FeedPL$Sex) == FALSE & FeedPL$Sex!= "", ]

#####Convert Measured Areas to Amounts Eaten#######
FeedPL$P1_Consumed <- FeedPL$Area1-FeedPL$Area2
FeedPL$P2_Consumed <- FeedPL$Area3-FeedPL$Area4
FeedPL$Tot_Consumed <- FeedPL$P1_Consumed + FeedPL$P2_Consumed

####### Assign Factors ######

FeedPL$Batch <- as.factor(FeedPL$Batch)
FeedPL$Site <- as.factor(FeedPL$Site)
FeedPL$Species <- as.factor(FeedPL$Species)
FeedPL$Sex <- as.factor(FeedPL$Sex)
FeedPL$Temp <- as.factor(FeedPL$Temp)
FeedPL$Light <- as.factor(FeedPL$Light)
FeedPL$TestTempF <- as.factor(FeedPL$TestTemp)
FeedPL$Treat <- paste(FeedPL$Temp, FeedPL$Light, sep = "_")

######## assess sample sizes ##########

table(FeedPL[, c(4,5,6, 21)]) #Sample sizes look pretty good


##### Plots of the data #####

#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}
Give.N <- function(x){return(c(y = MaxSE(x)*1.1, label = length(x))) }

#Total leaf area consumed
ggplot(data = FeedPL, aes(x = TestTemp, y = Tot_Consumed, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 5)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 5)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = 5), size = 5) +
  ylab("Leaf area consumed (cm2, 8h)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(15, 45), breaks = c(20, 30, 40))

#Leaf area consumed during first period (crop loading)
ggplot(data = FeedPL, aes(x = TestTemp, y = P1_Consumed, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 5)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 5)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = 5), size = 5) +
  ylab("Leaf area consumed (cm2, P1)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(15, 45), breaks = c(20, 30, 40))

#Leaf area consumed during second period (processing)
ggplot(data = FeedPL, aes(x = TestTemp, y = P2_Consumed, col = Site)) +
  rory_theme +
  facet_grid(Sex~Temp + Light) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 5)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 5)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = 5), size = 5) +
  ylab("Leaf area consumed (cm2, P2)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(15, 45), breaks = c(20, 30, 40))

##### Analyze performance with linear mixed models #####

hist(FeedPL$Tot_Consumed) #looks pretty normal

#scaled data to help lmer
X = scale(cbind(FeedPL$Tot_Consumed, FeedPL$TestTemp, FeedPL$Sex, FeedPL$Site, FeedPL$Temp, FeedPL$Light)) #julian date of first clutch

LmerF0 <- lmer(Tot_Consumed ~ poly(TestTemp, 2, raw = TRUE) * Sex * Site * Temp * Light + 
                 (1|ID/Mom), 
               REML = FALSE,
               na.action = 'na.omit', data = FeedPL) 

LmerF1 <- lmer(Tot_Consumed ~ TestTemp * Sex * Site * Temp * Light + 
                 (1|ID/Mom), 
               #na.action = 'na.omit',
               REML = FALSE,
               data = FeedPL) 

LmerF2 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^4 + 
                 (1|ID/Mom), 
               REML = FALSE,
               na.action = 'na.omit', data = FeedPL) 

LmerF2.1 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^4 - TestTemp:Sex:Temp:Light - TestTemp:Sex:Site:Light + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL) 

LmerF3 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^3 + 
                 (1|ID/Mom), 
               REML = FALSE,
               na.action = 'na.omit', data = FeedPL) 

LmerF3.1 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^3 +Sex:Temp:Light:Site + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF3.2 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 + TestTemp:Sex:Site + TestTemp:Sex:Temp + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF3.5 <- lmer(Tot_Consumed ~ (poly(TestTemp, 2, raw = TRUE) + Sex + Site + Temp + Light)^2 + TestTemp:Sex:Site + TestTemp:Sex:Temp + 
                   (1|ID/Mom), 
                 #REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF3.3 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 + TestTemp:Sex:Site + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF3.4 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 + TestTemp:Sex:Site - TestTemp:Light + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF4 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 + 
                 (1|ID/Mom), 
               REML = FALSE,
               na.action = 'na.omit', data = FeedPL)

LmerF4.1 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 + Site:Temp:Light + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF4.2 <- lmer(Tot_Consumed ~ (TestTemp + Sex + Site + Temp + Light)^2 - Temp:Light + 
                   (1|ID/Mom), 
                 REML = FALSE,
                 na.action = 'na.omit', data = FeedPL)

LmerF5 <- lmer(Tot_Consumed ~ TestTemp + Sex + Site + Temp + Light + 
                 (1|ID/Mom), 
               REML = FALSE,
               na.action = 'na.omit', data = FeedPL)

#AIC comparisons should be done on models with "REML = FALSE"
aictab(list(LmerF0, LmerF1, LmerF2, LmerF2.1, LmerF3, LmerF3.1, LmerF3.2, LmerF3.3, LmerF3.4, LmerF3.5, LmerF4, LmerF4.1, LmerF4.2, LmerF5), modnames = c("LmerF0", "LmerF1", "LmerF2", "LmerF2.1", "LmerF3", "Lmer3.1", "Lmer3.2", "Lmer3.3", "Lmer3.4", "Lmer3.5", "LmerF4", "LmerF4.1", "LmerF4.2", "LmerF5")) #LmerF1 is supported by far
#          K    AICc Delta_AICc AICcWt Cum.Wt       LL
#Lmer3.5  33 3542.78       0.00      1      1 -1736.30
#Lmer3.2  27 3561.44      18.66      0      1 -1752.33
#Lmer3.3  26 3562.36      19.58      0      1 -1753.89
#Lmer3.4  25 3563.92      21.14      0      1 -1755.77
#LmerF4.2 23 3564.08      21.29      0      1 -1758.03
#LmerF4   24 3566.10      23.32      0      1 -1757.95
#LmerF4.1 26 3567.44      24.66      0      1 -1756.43
#LmerF2.1 46 3570.72      27.94      0      1 -1735.23
#Lmer3.1  42 3571.17      28.38      0      1 -1740.16
#LmerF3   40 3571.31      28.53      0      1 -1742.56
#LmerF2   49 3577.61      34.83      0      1 -1735.10
#LmerF1   51 3580.20      37.42      0      1 -1733.99
#LmerF0   75 3601.19      58.40      0      1 -1714.08
#LmerF5   10 3611.82      69.04      0      1 -1795.71

#Preferred model (LmerF3.2) includes all 2-way interactions as well as interactions between TestTemp:Sex:Site and TestTemp:Sex:Temp

#P-values should be compared for models with REML = TRUE
Anova(LmerF1, type = 3) #The 4-way interaction between the TestTemp, Site, temp, and light is significant.  Everything is having effects, pretty much
Anova(LmerF2, type = 3)
Anova(LmerF2.1, type = 3)
Anova(LmerF3.2, type = 3)
Anova(LmerF3.5, type = 3)

#lsmeans (can add "adust = "none" if do not want p-value correction)
lsmeans(LmerF3.5, pairwise ~ Site|Temp+Sex) #A1 and C1 differ from one another in HV_S, LV_L, and LV_S
#B1 and C1 never differ
#A1 > C1 in all short-day females and >B1 in short-day, HV females, A1 < C1 in long-day, low-variance males 

lsmeans(LmerF3.5, pairwise ~ Light|Temp+Sex) #Light has no general effects
lsmeans(LmerF3.5, pairwise ~ Temp|Light+Sex) #Temp has no general effects
lsmeans(LmerF3.5, pairwise ~ Sex|Light+Temp) #Females are always higher

#plot the LSmeans
ObjF<-data.frame(summary(lsmeans(LmerF3.5, ~ TestTemp*Sex*Site*Temp*Light, at = list(TestTemp = c(20,30,40)), type = "response")))
#Obj0<-data.frame(summary(lsmeans(LmerH0, ~ TestTemp*Sex*Site*Temp*Light, at = list(TestTemp = c(25,17,10,35)), type = "response")))

ggplot(data = ObjF, aes(y = lsmean, x = TestTemp, col = Site, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = asymp.UCL, ymin = asymp.LCL), position = position_dodge(width = 3)) +
  geom_line(aes(x = TestTemp), size = 1.5, position = position_dodge(width = 3)) +
  geom_point(size = 6, position = position_dodge(width = 3)) +
  #ylab("Leaf area consumed\n(cm2, LS-means)") + xlab("Site") +
  xlab(expression(paste("Test temperature (", degree~C, ")"))) +
  ylab(expression(atop("Leaf area consumed", paste("(c",m^2, "/8h, LS-means)", sep="")))) +
  scale_x_continuous(limits = c(15, 45), breaks = c(20, 30, 40)) +
  scale_shape_manual(values = c(19, 1))

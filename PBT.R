rm(list=ls(all = TRUE))

library(stringr)
library(foreach)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)

source('Functions.R') #PC, can do via R studio as well

####### Input files ##########

setwd("./data/")

#upload PBT data
PBT <- read.csv(file = "PBTdata_06202016_ALL.csv", header = TRUE, stringsAsFactors = FALSE)

str(PBT)

###### Clean up the data ############

#collect treatment labels
Labels <- str_split_fixed(PBT[,1], "_", n = 4) #collect Temp and Light treatments
PBT$Temp <- Labels[,2]
PBT$Light <- Labels[,3]
PBT$Site <- Labels[,4]

HopRaw <- read.csv(file = "Pdodge_PerformanceData1_06202016_ALL.csv", header = TRUE, stringsAsFactors = FALSE)
SexData <- data.frame(HopRaw$Individual, HopRaw$sex, stringsAsFactors = FALSE)

PBT <- merge(PBT, SexData, by.x =  colnames(PBT)[1], by.y = colnames(SexData)[1])
colnames(PBT)[27] <- "Sex"
PBT$Sex <- ifelse(PBT$Sex == "", NA, PBT$Sex)


#trim out missing data (rows with no values)
PBT <- PBT[is.na(PBT$t_enter_grad) == FALSE, ]

######## calculate median and mean TB for each individual ###############

TB_df <- data.frame(PBT$Tb_1, PBT$Tb_2, PBT$Tb_3, PBT$Tb_4, PBT$Tb_5, PBT$Tb_6)

PBT$Tb_Mean <- apply(TB_df, 1, mean, na.rm = TRUE)
PBT$Tb_Median <- apply(TB_df, 1, median, na.rm = TRUE)
PBT$Tb_sd <- apply(TB_df, 1, sd, na.rm = TRUE)

######## Sample Size ##########

PBT$Treat <- paste(PBT$Temp, PBT$Light, sep = "_")
table(PBT[,c(26,27,31)]) #results collated in google drive file "Performance Sample Sizes"


######### Initial Plots ###############
#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}

Give.N <- function(x){
  return(c(y = MaxSE(x)*1.01, label = length(x))) 
}

ggplot(data = PBT, aes(x = Site, y = Tb_Mean, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Mean PBT (C)") + xlab("Site")


#Boxplot (less useful)
Give.N3 <- function(x){
  return(c(y = 48, label = length(x))) 
}

ggplot(data = PBT, aes(x = Site, y = Tb_Mean, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_boxplot() +
  #stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6), geom = "box") +
  #stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(fun.data = Give.N3, geom = "text", fun.y = median, position = position_dodge(width = .8), size = 7) +
  ylab("Mean PBT (C)") + xlab("Site")

##### Analyses #####
PBT_2 <- merge(PBT, Instars[,c(1,3,54)], by.x = "Individual", by.y = "Individual")
PBT_Long <- melt(PBT_2[,c(1,3,10,12,14,16,18,20,24:27,31,33)], id.var = c("Individual", "TprefBATCH", "Female2", "Temp", "Light", "Site", "Sex", "Treat"))
PBT_Long$Temp <- as.factor(PBT_Long$Temp)
PBT_Long$Light <- as.factor(PBT_Long$Light)
PBT_Long$Site <- as.factor(PBT_Long$Site)
PBT_Long$Sex <- as.factor(PBT_Long$Sex)
PBT_Long$Individual <- as.factor(PBT_Long$Individual)


hist(PBT_Long$value) #looks quite normal

LmerPBT.0 <- lmer(value ~ Sex * Site * Temp * Light + 
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long) #initially attempted to include batch and female, but so little variation that the scales were off and the model showed errors

LmerPBT.1 <- lmer(value ~ Sex + Site + Temp + Light + 
                    Sex:Site + Sex:Temp + Sex:Light + Site:Temp + Site:Light + Temp:Light + 
                    Site:Light:Temp +
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

LmerPBT.2 <- lmer(value ~ Sex + Site + Temp + Light + 
                    Sex:Site + Sex:Temp + Sex:Light + Site:Temp + Site:Light + Temp:Light + 
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

LmerPBT.3 <- lmer(value ~ Sex + Site + Temp + Light + 
                    Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

LmerPBT.4 <- lmer(value ~ Site + Temp + Light + 
                    Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

LmerPBT.5 <- lmer(value ~ Site + Temp + Light + 
                    Site:Temp + Site:Light + Temp:Light +
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

LmerPBT.6 <- lmer(value ~ Sex + Site + Temp + Light + 
                    (1|Individual), na.action = 'na.omit', 
                  #REML = FALSE, 
                  data = PBT_Long)

LmerPBT.7 <- lmer(value ~ Site + Temp + Light + 
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long)

aictab(list(LmerPBT.0, LmerPBT.1, LmerPBT.2, LmerPBT.3, LmerPBT.4, LmerPBT.5, LmerPBT.6, LmerPBT.7), modnames = c("LmerPBT.0", "LmerPBT.1", "LmerPBT.2", "LmerPBT.3", "LmerPBT.4", "LmerPBT.5", "LmerPBT.6", "LmerPBT.7"))
#LmerPBT.6 is the preferred model, delta AIC = 7.14 (no interactions)

summary(LmerPBT.3)
Anova(LmerPBT.6, type = 3)

lsmeans(LmerPBT.6, pairwise ~ Site|Temp+Light) #Sites never significantly differ from one another
lsmeans(LmerPBT.6, pairwise ~ Light|Temp) #L and S differ in LV only
lsmeans(LmerPBT.6, pairwise ~ Temp|Light) #HV and LV differ, only in long days and primarilly for C1

Obj_PBT<-data.frame(summary(lsmeans(LmerPBT.6, ~ Sex * Site * Temp * Light)))

ggplot(data = Obj_PBT, aes(y = lsmean, x = Site, shape = Sex, linetype = Sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), position = position_dodge(width = .6)) +
  geom_line(aes(x = as.numeric(as.factor(Site))), size = 1.5, position = position_dodge(width = .6)) +
  geom_point(size = 6, position = position_dodge(width = .6)) +
  ylab("Preferred body temperature\n(C, LS-means)") + xlab("Site")

#estimate lsmean PBT for every individual

#Unfortunately, can't figure out how to collect individual means from preferred model, so have to just use means from a model only including individual
Lm_PBT.6 <- lm(value ~ Individual, data = PBT_Long)
Obj_PBT_Ind <- data.frame(summary(lsmeans(Lm_PBT.6, ~ Individual)))
colnames(Obj_PBT_Ind) <- c("Individual", "PBT", "PBT.SE", "df", "lower.CL",  "upper.CL")  

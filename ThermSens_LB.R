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

#Process PBT
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

#calculate median and mean TB for each individual 
TB_df <- data.frame(PBT$Tb_1, PBT$Tb_2, PBT$Tb_3, PBT$Tb_4, PBT$Tb_5, PBT$Tb_6)

PBT$Tb_Mean <- apply(TB_df, 1, mean, na.rm = TRUE)
PBT$Tb_Median <- apply(TB_df, 1, median, na.rm = TRUE)
PBT$Tb_sd <- apply(TB_df, 1, sd, na.rm = TRUE)

#Sample Size 

PBT$Treat <- paste(PBT$Temp, PBT$Light, sep = "_")
table(PBT[,c(26,27,31)]) #results collated in google drive file "Performance Sample Sizes"

#Drop cases with no sex
PBT= subset(PBT, PBT$Sex %in% c("M","F"))

#----------------
#upload CT data
CT <- read.csv(file = "CT_Data_2016.csv", header = TRUE, stringsAsFactors = FALSE)

#=========================================
# PLOTS

#PBT
#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}

Give.N <- function(x){
  return(c(y = MaxSE(x)*1.01, label = length(x))) 
}

ggplot(data = PBT, aes(x = Site, y = Tb_Mean, color = Temp)) +
  rory_theme +
  facet_wrap(~Sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Mean PBT (C)") + xlab("Site")+scale_color_viridis_d()

#Analyses
PBT_2 <- merge(PBT, Instars[,c(1,3,54)], by.x = "Individual", by.y = "Individual")
PBT_Long <- melt(PBT_2[,c(1,3,10,12,14,16,18,20,24:27,31,33)], id.var = c("Individual", "TprefBATCH", "Female2", "Temp", "Light", "Site", "Sex", "Treat"))
PBT_Long$Temp <- as.factor(PBT_Long$Temp)
PBT_Long$Light <- as.factor(PBT_Long$Light)
PBT_Long$Site <- as.factor(PBT_Long$Site)
PBT_Long$Sex <- as.factor(PBT_Long$Sex)
PBT_Long$Individual <- as.factor(PBT_Long$Individual)

#Check Rory's file for model selection, check random effect
LmerPBT.0 <- lmer(value ~ Sex * Site * Temp * Light + 
                    (1|Individual), na.action = 'na.omit', 
                  REML = FALSE, 
                  data = PBT_Long) #initially attempted to include batch and female, but so little variation that the scales were off and the model showed errors

LmerPBT.6 <- lmer(value ~ Sex + Site + Temp + Light + 
                    (1|Individual), na.action = 'na.omit', 
                  #REML = FALSE, 
                  data = PBT_Long)

aictab(list(LmerPBT.0, LmerPBT.6), modnames = c("LmerPBT.0","LmerPBT.6"))

summary(LmerPBT.6)
Anova(LmerPBT.6, type = 3)

lsmeans(LmerPBT.6, pairwise ~ Site|Temp+Light) #Sites never significantly differ from one another
lsmeans(LmerPBT.6, pairwise ~ Light|Temp) #L and S differ in LV only
lsmeans(LmerPBT.6, pairwise ~ Temp|Light) #HV and LV differ, only in long days and primarilly for C1

#CT



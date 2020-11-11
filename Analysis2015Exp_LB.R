#ANALYZE 2015 data

##### Grasshopper 2015 Thermal Biology #####
# Using Sima's data to examine CTmin, CTmax, and PBT
# Initial coding by Rory Telemeco 18 May 2015

rm(list=ls(all = TRUE))

#libraries
library(ggplot2)
require(grid)
library(gridExtra)
library(lsmeans)
library(reshape2)
library(car)
library(stringr)
library(nlme)
library(lme4)

#my common custom functions/code
rory_theme <- theme_classic(base_size = 28) + theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) + #adjust axis and text size and position
  theme(axis.line=element_line(size = 1.25), plot.margin = unit(c(.3,.3,.6,.6), "cm"), axis.ticks=element_line(size = 1.25)) + #adjust plot margins and line element size
  theme(axis.ticks.length =unit(-0.3, "cm"), axis.ticks.margin = unit(0.7, "cm")) + #move tick marks to the inside
  theme(panel.margin = unit(2, units = "lines")) #+

Validate <- function(Mod, main) {
  if(missing(main)) main <- deparse(substitute(Mod))
  E<-resid(Mod)
  Fit <- fitted(Mod)
  op <- par(mfrow =c(1,2))
  plot(x=Fit, y=E,
       xlab = "Fitted values", ylab = "Residuals",
       main = paste("Residuals vs Fitted Values from", main))
  hist(E, nclass = 15, 
       xlab = "Residuals",
       main = paste("Histrogram of residuals from", main))
  par(op) 
} #graphical model validation code

##### bookkeeping for all #####
setwd("./data/")

#Thermal Data from Sima
Thermal <- read.csv("GrasshopperData_Tpref_CT_SPR15.csv", header=T, na.strings = c("", " ", "NA"), stringsAsFactors = FALSE)

#merged data from Lauren
clean <- read.csv("GH_3rdAdult_dataALL_2015Experiment.csv", header=T, na.strings = c("", " ", "NA"), stringsAsFactors = FALSE)

#create identifier in same format as what is in Thermal
Pod_ID2 <- paste(clean$Site, clean$Fem_ID, sep = "")
Pod_ID2 <- paste(Pod_ID2, clean$Treatment, sep = ".")
clean$Pod_ID2 <- Pod_ID2

#Pull out identifying data
cleanID <- clean[, c(2,3,4,6,8,47)] #just select out the identifiers
cleanID_red <- cleanID[!duplicated(cleanID$Pod_ID2),] #reduce to a single row per egg pod

#merge Thermal and clean data frames
Therm2 <- merge(Thermal, cleanID_red, by.y = "Pod_ID2", by.x = "Pod_ID", all.x = TRUE, all.y = FALSE)
Therm2$Temp <- ifelse(Therm2$Treatment == "24L" | Therm2$Treatment == "24S", "24C", "28C")
Therm2$Light <- ifelse(Therm2$Treatment == "24L" | Therm2$Treatment == "28L", "L", "S")

#generally explore sample sizes
table(Therm2[, c(9,12,11)]) #very uneven sampling...

#====================================
#Development plots

# Format dates and calculate ages 
clean$Date1_Hatch <- as.Date(clean$Date_Hatch, format = "%Y-%m-%d")
clean$Date1_Adult <- as.Date(clean$Date_Adult, format = "%Y-%m-%d")

clean$Age_Adult <- clean$Date1_Adult - clean$Date1_Hatch
clean$Age_Adult= as.numeric(clean$Age_Adult)

#split factors
clean$temp= substr(clean$Treatment, 1, 2)
clean$light= substr(clean$Treatment, 3, 4)

#restrict data
clean= subset(clean, clean$Site %in% c("A1","B1","C1"))
clean= subset(clean, clean$temp %in% c("24","28"))

#plot
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}

Give.N <- function(x){
  return(c(y = MaxSE(x)*1.01, label = length(x))) 
}

ggplot(data = clean, aes(x = Site, y = Age_Adult, color = temp)) +
  facet_grid(Species~light) + geom_point()+
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  ylab("Time to adult (d)") + xlab("Site")+scale_color_viridis_d()

#------------------
#simple temp and light models

mod1= lm(Age_Adult~temp*Site*light+Species, data=clean)
mod1= lm(Age_Adult~temp*Site*light, data=clean[clean$Species=="dodg",])

mod2 <- lmer(Age_Adult~temp*Site*light + 
                 (1|Fem_ID), na.action = 'na.omit', REML=FALSE, data=clean[clean$Species=="dodg",])
Anova(mod2)

#Discuss model more temp and light? Include supplementary figure?

#====================================
#TPC plots

Therm2L <- melt(Therm2, id.vars = c("Pod_ID", "Species.x", "Date", "CTmin_max_Note", "Individual_ID", "Site", "Fem_ID", "Species.y", "Treatment", "Temp", "Light"))
ThermDodgL <- melt(ThermDodg, id.vars = c("Pod_ID", "Species.x", "Date", "CTmin_max_Note", "Individual_ID", "Site", "Fem_ID", "Species.y", "Treatment", "Temp", "Light"))

ggplot(data = Therm2L, aes(x = Site, y = as.numeric(value), color = Temp, shape=variable)) +
  facet_grid(Light~Species.y) + geom_point()+
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  ylab("Temperature (C)") + xlab("Site")+scale_color_viridis_d()

#models
M_PBT1 <- lm(PBT ~ Site + Temp + Light, data = ThermDodg) #Temp Signficant
Anova(M_PBT1, type = 3)
#PBT higher at 28

M_MAX1 <- lm(CTmax ~ Site + Temp + Light, data = ThermDodg) #Temp Significant
Anova(M_MAX1, type = 3)
#CTmax higher at 28

M_MIN1 <- lm(CTmin ~ Site + Temp + Light, data = ThermDodg) #Nothing significant
Anova(M_MIN1, type = 3)

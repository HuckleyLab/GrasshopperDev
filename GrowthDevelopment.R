#2016 growth rate and survival analyses
# Author: Rory Telemeco

rm(list=ls(all = TRUE))

library(stringr)
library(foreach)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)
library(nlme)
library(lme4)
library(car)
library(AICcmodavg)
library(lsmeans)
library(reshape2)
library(MASS)

source('Functions.R')

####### Input files ##########

setwd("./data/")

Instars <- read.csv(file = "GrasshopperData_3rdToAdult_SPR16.csv", header = TRUE, stringsAsFactors = FALSE)
str(Instars)

Instars$sex <- ifelse(Instars$sex == "", NA, Instars$sex)

##### Format dates and calculate ages ########

Instars$Date1_Hatch <- as.Date(Instars$Date_Hatch, format = "%m/%d/%Y")
Instars$Date1_1st <- as.Date(Instars$Date_1st, format = "%m/%d/%Y")
Instars$Date1_2nd <- as.Date(Instars$Date_2nd, format = "%m/%d/%Y")
Instars$Date1_3rd <- as.Date(Instars$Date_3rd, format = "%m/%d/%Y")
Instars$Date1_4th <- as.Date(Instars$Date_4th, format = "%m/%d/%Y")
Instars$Date1_5th <- as.Date(Instars$Date_5th, format = "%m/%d/%Y")
Instars$Date1_5th. <- as.Date(Instars$Date_5th., format = "%m/%d/%Y")
Instars$Date1_Adult <- as.Date(Instars$Date_Adult, format = "%m/%d/%Y")
Instars$Date1_Death <- as.Date(Instars$Date_Death, format = "%m/%d/%Y")

Instars$Age_3rd <- Instars$Date1_3rd - Instars$Date1_Hatch
Instars$Age_4th <- Instars$Date1_4th - Instars$Date1_Hatch
Instars$Age_5th <- Instars$Date1_5th - Instars$Date1_Hatch
Instars$Age_Adult <- Instars$Date1_Adult - Instars$Date1_Hatch


#trim out pell data

Instars <- Instars[Instars$Species == "dodg",]

##### Plot time to adulthood ######

#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}

Give.N <- function(x){
  return(c(y = MaxSE(x)*1.01, label = length(x))) 
}

ggplot(data = Instars, aes(x = Site, y = Age_Adult, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Time to adult (d)") + xlab("Site")





######### Initial analysis ##################

#set factors

Instars$Site <- as.factor(Instars$Site)
Instars$sex <- as.factor(Instars$sex)
Instars$Temp <- as.factor(Instars$Temp)
Instars$Light <- as.factor(Instars$Light)
Instars$Female <- as.factor(Instars$Female)
Instars$Female2 <- as.factor(paste(Instars$Female, Instars$Site, sep = "_"))

hist(as.numeric(Instars$Age_Adult)) #not quite normal, should probably use a Poisson distribution...

Lmer0 <- glmer(as.numeric(Age_Adult) ~ sex * Site * Temp * Light + 
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

Lmer1 <- glmer(as.numeric(Age_Adult) ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                 Site:Light:Temp +
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

Lmer2 <- glmer(as.numeric(Age_Adult) ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

Lmer3 <- glmer(as.numeric(Age_Adult) ~ sex + Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

Lmer4 <- glmer(as.numeric(Age_Adult) ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

Lmer5 <- glmer(as.numeric(Age_Adult) ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light +
                 (1|Female2), family = poisson, na.action = 'na.omit', data = Instars)

aictab(list(Lmer0, Lmer1, Lmer2, Lmer3, Lmer4, Lmer5), modnames = c("Lmer0", "Lmer1", "Lmer2", "Lmer3", "Lmer4", "Lmer5"))
#Lmer3 is by far the preferred model according to AIC

summary(Lmer3)
Anova(Lmer3, type = 3)
Validate(Lmer3) #looks great

#look at individual contrasts
lsmeans(Lmer3, pairwise ~ Site|Temp+Light) #A1 and C1 differ from one another in HV_S, LV_L, and LV_S
lsmeans(Lmer3, pairwise ~ Light|Temp) #Light generally has an effect in LV but not HV
lsmeans(Lmer3, pairwise ~ Temp|Light+Site) #Temp regime had an effect under every light and site, with LV always greater than HV

#plot lsmeans
Obj1<-data.frame(summary(lsmeans(Lmer3, ~ sex*Site*Temp*Light, type = "response")))

ggplot(data = Obj1, aes(y = response, x = Site, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = asymp.UCL, ymin = asymp.LCL), position = position_dodge(width = .6)) +
  geom_line(aes(x = as.numeric(as.factor(Site))), size = 1.5, position = position_dodge(width = .6)) +
  geom_point(size = 6, position = position_dodge(width = .6)) +
  ylab("Time to adult (d, LS-means)") + xlab("Site")


#check lmer for singularity, no singulatiry with lmer, OK to assume normal?
Lmer3 <- lmer(as.numeric(Age_Adult) ~ sex + Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                 (1|Female2), na.action = 'na.omit', data = Instars)


#plot predictions
plot_model(Lmer3, type="pred", terms=c("Site","Temp","Light"), show.data=TRUE)

####### Analyze maturation rate while including all instars #########

Instars_Adult <- Instars[is.na(Instars$Age_Adult) == FALSE,]
Instars_Adult_Long <- melt(Instars_Adult[,c(1:8,50:54)], id.vars = c("Individual", "CCode", "Female", "Female2", "Site", "Species", "sex", "Temp", "Light"))

Lmer0.1 <- glmer(as.numeric(value) ~ sex * Site * Temp * Light + 
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

Lmer1.1 <- glmer(as.numeric(value) ~ sex + Site + Temp + Light + 
                   sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                   Site:Light:Temp +
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

Lmer2.1 <- glmer(as.numeric(value) ~ sex + Site + Temp + Light + 
                   sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

Lmer3.1 <- glmer(as.numeric(value) ~ sex + Site + Temp + Light + 
                   Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

Lmer4.1 <- glmer(as.numeric(value) ~ Site + Temp + Light + 
                   Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

Lmer5.1 <- glmer(as.numeric(value) ~ Site + Temp + Light + 
                   Site:Temp + Site:Light + Temp:Light +
                   (1|Individual/variable/Female2), family = poisson, na.action = 'na.omit', data = Instars_Adult_Long)

aictab(list(Lmer0.1, Lmer1.1, Lmer2.1, Lmer3.1, Lmer4.1, Lmer5.1), modnames = c("Lmer0.1", "Lmer1.1", "Lmer2.1", "Lmer3.1", "Lmer4.1", "Lmer5.1"))
#3.1 is the preferred model, suggesting that that sex has a main effect but no interaction, but that the other interactions are important

summary(Lmer3.1)
Anova(Lmer3.1)


##### Plot the instar-specific maturation data #####
levels(Instars_Adult_Long$variable) <- c("3rd", "4th", "5th", "Adult")


ggplot(data = Instars_Adult_Long, aes(x = variable, y = value, shape = sex, linetype = sex, col = Site)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(variable))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Time (d)") + xlab("Instar")



#### Mass at adult #####
ggplot(data = Instars, aes(x = Site, y = Mass_Adult, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Mass as adult (g)") + xlab("Site")
# appears to be a strong sex effect (duh), and a trend for high-elevation to be smaller


hist(Instars$Mass_Adult) #also looks a little skewed, Poisson might be best

LmerM0 <- lmer(Mass_Adult ~ sex * Site * Temp * Light + 
                 (1|Female2), 
               REML = FALSE,
               na.action = 'na.omit', data = Instars)

LmerM1 <- lmer(Mass_Adult ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                 Site:Light:Temp +
                 (1|Female2), 
               #REML = FALSE, 
               na.action = 'na.omit', data = Instars)

LmerM2 <- lmer(Mass_Adult ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                 (1|Female2), 
               #REML = FALSE,
               na.action = 'na.omit', data = Instars)

LmerM3 <- lmer(Mass_Adult ~ sex + Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                 (1|Female2), 
               REML = FALSE,
               na.action = 'na.omit', data = Instars)

LmerM4 <- lmer(Mass_Adult ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + Site:Temp:Light +
                 (1|Female2), 
               REML = FALSE, 
               na.action = 'na.omit', data = Instars)

LmerM5 <- lmer(Mass_Adult ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light +
                 (1|Female2), 
               REML = FALSE,
               na.action = 'na.omit', data = Instars)

aictab(list(LmerM0, LmerM1, LmerM2, LmerM3, LmerM4, LmerM5), modnames = c("LmerM0", "LmerM1", "LmerM2", "LmerM3", "LmerM4", "LmerM5"))
#M2 is slightly preferred over 3 and 0 which are greatly preferred over the rest
Validate(LmerM2) #looks good

Anova(LmerM2, type = 3) #sex:Temp and it's main effects are significant

lsmeans(LmerM1, pairwise ~ Site|Temp+sex+Light, adjust = "none") #A1 and C1 differ from one another in HV_S, LV_L, and LV_S
lsmeans(LmerM1, pairwise ~ Light|Temp+sex, adjust = "none") #Light generally has an effect in LV but not HV
lsmeans(LmerM1, pairwise ~ Temp|Light+sex) #Temp regime had an effect under every light and site, with LV always greater than HV

##### Body size/shape as adult #####
ggplot(data = Instars, aes(x = Site, y = Femur_Adult, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Femur length (mm)") + xlab("Site")

ggplot(data = Instars, aes(x = Site, y = Pronotum_Adult, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Pronotum length (mm)") + xlab("Site")


#plot of mass via lsmeans
Obj2<-data.frame(summary(lsmeans(LmerM1, ~ sex*Site*Temp*Light, type = "response")))

ggplot(data = Obj2, aes(y = lsmean, x = Site, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), position = position_dodge(width = .6)) +
  geom_line(aes(x = as.numeric(as.factor(Site))), size = 1.5, position = position_dodge(width = .6)) +
  geom_point(size = 6, position = position_dodge(width = .6)) +
  ylab("Adult mass (g, LS-means)") + xlab("Site")


####### Multivariate Shape (simple) ##########

AdultsOnly <- Instars[is.na(Instars$Femur_Adult) == FALSE, ]
AdultsOnly <- AdultsOnly[is.na(AdultsOnly$sex) == FALSE, ]

TS <- paste(AdultsOnly$Temp, AdultsOnly$Light, AdultsOnly$sex, AdultsOnly$Site, sep = "_")

Shape <- cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult)


#run the PCA
pca1 <- prcomp(Shape, scale. = TRUE)  #uses SVD for computations (mean-centered is default), scaled transforms data to unit variance
PC.scores1<-data.frame(pca1$x)

#Calculate means and se of PC1 and PC2 for plotting
PC1_mean <- unlist(tapply(PC.scores1[, "PC1"], TS, mean))
PC2_mean <- unlist(tapply(PC.scores1[, "PC2"], TS, mean))
PC1_se <- unlist(tapply(PC.scores1[, "PC1"], TS, std.err))
PC2_se <- unlist(tapply(PC.scores1[, "PC2"], TS, std.err))
T_L_Sex_Site <- rownames(PC1_mean)  

# Data frame for the means and errors
PCmeans <- data.frame(T_L_Sex_Site, PC1_mean, PC2_mean, PC1_se, PC2_se)
PCLabels <- str_split_fixed(PCmeans[,1], "_", n = 4) #collect Temp and Light treatments
PCmeans$Temp <- PCLabels[,1]
PCmeans$Light <- PCLabels[,2]
PCmeans$Sex <- PCLabels[,3]
PCmeans$Site <- PCLabels[,4]

error_PC1 <- aes(xmax = PC1_mean + (2*PC1_se), xmin= PC1_mean - (2*PC1_se)) #95% CI
error_PC2 <- aes(ymax = PC2_mean + (2*PC2_se), ymin= PC2_mean - (2*PC2_se)) #95% CI

#PC limits
PC1_lim <- max(abs(PC.scores1[,"PC1"]))
PC2_lim <- max(abs(PC.scores1[,"PC2"]))
PC_lim <- max(c(PC1_lim, PC2_lim))

#reduced dataframe for scatterplot
ScatDat <- data.frame(PC.scores1[, "PC1"], PC.scores1[, "PC2"])

#PC plot
ggplot(data = PCmeans, aes(x = PC1_mean, y = PC2_mean, col = Site, shape = Sex)) +
  facet_wrap(~Temp+Light, nrow = 2) +
  rory_theme +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(error_PC2, size = 1, width = 0) +
  geom_errorbarh(error_PC1, size = 1, height = 0) +
  geom_point(size = 8) +
  xlab("PC1") + ylab("PC2") +
  ylim(-PC_lim, PC_lim) + xlim(-PC_lim, PC_lim) +
  #scale_color_manual(values = c(Color1, Color2)) +
  scale_shape_manual(values = c(1, 19)) #+
#scale_fill_manual(values = c(Color1,Color2,"white","white")) +
#geom_point(data = ScatDat, aes(x = PC.scores1...xPC., y = PC.scores1...yPC.), col = colo1, shape = symb1, fill = bg1, size = 4) +
#theme(legend.position = Pos)#
# 


#### multivariate shape analysis ####
# not sure if or how to include a random effect.  For this, leaving female out of it

LmS0 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ sex * Site * Temp * Light, 
           data = AdultsOnly)

LmS1 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ sex + Site + Temp + Light + 
             sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
             Site:Light:Temp, 
           data = AdultsOnly)

LmS2 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ sex + Site + Temp + Light + 
             sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light,
           data = AdultsOnly)

LmS3 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ sex + Site + Temp + Light + 
             Site:Temp + Site:Light + Temp:Light + Site:Temp:Light, 
           data = AdultsOnly)

LmS4 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ Site + Temp + Light + 
             Site:Temp + Site:Light + Temp:Light + 
             Site:Temp:Light,
           data = AdultsOnly)

LmS5 <- lm(cbind(AdultsOnly$Mass_Adult^1/3, AdultsOnly$Femur_Adult, AdultsOnly$Pronotum_Adult) ~ Site + Temp + Light + 
             Site:Temp + Site:Light + Temp:Light,
           data = AdultsOnly)

Anova(LmS0, type = 3) #all interactions are non signficant
Anova(LmS1, type = 3) #an interaction between sex and temp
Anova(LmS2, type = 3) #interactions between sex and temp, site and temp, site and light (probably a good idea to keep all of the two-way interactions)
Anova(LmS3, type = 3) #no interactions significant
Anova(LmS4, type = 3) #no interactions significant
Anova(LmS5, type = 3) #no interactions signficant
# in all, only the effects of sex and site appear important

lsmeans(LmS2, pairwise ~ Site|Temp+sex+Light) #A1 and C1 differ from one another in HV_S, LV_L, and LV_S
lsmeans(LmS2, pairwise ~ Light|Temp+sex+Site) #Light generally has an effect in LV but not HV
lsmeans(LmS2, pairwise ~ Temp|Light+sex+Site) #Temp regime had an effect under every light and site, with LV always greater than HV


#Can't use AIC with mutlivariate so use p-values

##### What about survival? #####
Instars$Survive <- ifelse(is.na(Instars$Age_Adult) == TRUE, 0, 1)

ggplot(data = Instars, aes(x = Site, y = Survive, shape = sex, linetype = sex)) +
  rory_theme +
  facet_wrap(~Temp+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Survival frequency") + xlab("Site")


MS0 <- glmmPQL(Survive ~ sex * Site * Temp * Light, 
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

MS1 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light + 
                 Site:Light:Temp,
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

MS2 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Temp + sex:Light + Site:Temp + Site:Light + Temp:Light,
               random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS3 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + 
                 Site:Temp:Light, 
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

MS4 <- glmmPQL(Survive ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + 
                 Site:Temp:Light,
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

MS5 <- glmmPQL(Survive ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light + 
                 Site:Temp:Light,
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

MS6 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 sex:Site + sex:Light + Site:Temp + Site:Light + Temp:Light,
               random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS7 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light,
               random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS8 <- glmmPQL(Survive ~ Site + Temp + Light + 
                 Site:Temp + Site:Light + Temp:Light,
               random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS9 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                 Site:Temp + Temp:Light,
               random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS10 <- glmmPQL(Survive ~ sex + Site + Temp + Light + 
                  Site:Temp,
                random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

MS11 <- glmmPQL(Survive ~ sex + Site + Temp + 
                  Site:Temp,
                random = (~ 1|Female2), family = binomial, data = Instars, niter = 100)

Anova(MS0, type = 3)
Anova(MS1, type = 3)
Anova(MS2, type = 3)
Anova(MS3, type = 3)
Anova(MS4, type = 3)
Anova(MS5, type = 3)
Anova(MS6, type = 3)
Anova(MS7, type = 3)
Anova(MS8, type = 3)
Anova(MS9, type = 3)
Anova(MS10, type = 3)
Anova(MS11, type = 3)

#==================
#DEVELOPMENT PLOT
#load and process data in Growth Develoment.R

#try ordered factor
#Instars$Site= factor(Instars$Site, ordered=TRUE, levels=c("A1","B1","C1") )
Instars$Site <- as.factor(Instars$Site)

nLmer0 <- lmer(as.numeric(Age_Adult) ~ sex * Site * Temp * Light + 
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)


#split sexes
#nLmer0 <- lmer(as.numeric(Age_Adult) ~ Site * Temp * Light + 
#                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars[Instars$sex=="F",])

#check numeric elevation
#make numeric elevation variable
Instars$elev= 2195
Instars$elev[Instars$Site=="B1"]= 2591
Instars$elev[Instars$Site=="C1"]= 3048
#scale
Instars$elev.s= scale(Instars$elev)

nLmer0 <- lmer(as.numeric(Age_Adult) ~ sex * elev.s * Temp * Light + 
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)
#----

Anova(nLmer0, type = 3)
summary(nLmer0)

#plot marginal effects
#https://strengejacke.github.io/sjPlot/articles/plot_marginal_effects.html

# 4 way interaction
plot_model(nLmer0, type="pred", terms=c("Site","Temp","Light","sex"), show.data=FALSE)

# main effects
plot_model(nLmer0, type="pred", terms=c("sex"), show.data=FALSE)
plot_model(nLmer0, type="pred", terms=c("Temp"), show.data=FALSE)
plot_model(nLmer0, type="pred", terms=c("Light"), show.data=FALSE)

#Interpretation:
#sex:site:temp:light=
# LV slower development 
# Long days delay development for LV, speed development for HV males but not females
# high elevation (C1) develops faster, except for long day females:
# Male: main effect develops slightly faster

#-----------------
#Figure 1

#plot based on individual means
MinSE <- function(x) {mean(x, na.rm = TRUE) - std.error(x, na.rm = TRUE)}
MaxSE <- function(x) {mean(x, na.rm = TRUE) + std.error(x, na.rm = TRUE)}

Give.N <- function(x){
  return(c(y = MaxSE(x)*1.01, label = length(x))) 
}

#drop NA sex
Instars= Instars[which(Instars$sex %in% c("F","M")),]

ggplot(data = Instars, aes(x = Site, y = Age_Adult, color = Temp)) +
  rory_theme +
  facet_wrap(~sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Time to adult (d)") + xlab("Site")+scale_color_viridis_d()

#==================
#Survival

Instars$Survive <- ifelse(is.na(Instars$Age_Adult) == TRUE, 0, 1)

ggplot(data = Instars, aes(x = Site, y = Survive, color = Temp)) +
  rory_theme +
  facet_wrap(~sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Survival (%)") + xlab("Site")+scale_color_viridis_d()

mod.surv <- lmer(Survive ~ sex * Site * Temp * Light + 
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)

Anova(mod.surv, type = 3)
summary(mod.surv)

MS0 <- glmmPQL(Survive ~ sex * Site * Temp * Light, 
               random = (~ 1|Female2), 
               family = binomial, data = Instars, niter = 100)

Anova(MS0, type = 3)

#==================
#MASS

#didn't go through model selection, but may make sense to use full model
mass0 <- lmer(as.numeric(Mass_Adult) ~ sex * Site * Temp * Light + 
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)

Anova(mass0, type = 3)
summary(mass0)

#plot marginal effects
# 3 way interaction
plot_model(mass0, type="pred", terms=c("Site","Temp","sex"), show.data=FALSE)

# 4 way interaction
plot_model(mass0, type="pred", terms=c("Site","Temp","Light","sex"), show.data=FALSE)

#--------------------
#Figure 2
#Mass at adulthood
ggplot(data = Instars, aes(x = Site, y = Mass_Adult, color=Temp)) +
  rory_theme +
  facet_wrap(~sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Mass as adult (g)") + xlab("Site")
# appears to be a strong sex effect (duh), and a trend for high-elevation to be smaller

#==================
#Shape

#to long format
shape_long <- melt(Instars[,c(1:8,54,30:32)], id.vars = c("Individual", "CCode", "Female", "Female2", "Site", "Species", "sex", "Temp", "Light","Mass_Adult")) 

ggplot(data = shape_long, aes(x = Mass_Adult^0.333, y = value, color=Temp)) +
  geom_point(aes(shape=sex))+facet_grid(variable ~ Light, scales="free_y")+geom_smooth(method="lm")

#STATS check for divergence from mass
shape0 <- lmer(Femur_Adult ~ Mass_Adult* Temp * Light * sex * Site +
                (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)

shape0 <- lmer(log(Pronotum_Adult) ~ log(Mass_Adult)* Temp * Light * sex * Site +
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)

Anova(shape0, type = 3)
summary(shape0)

# 4 way interaction
plot_model(shape0, type="pred", terms=c("Site","Temp","Light","sex"), show.data=FALSE)

#----------
#femur length
ggplot(data = Instars, aes(x = Site, y = Femur_Adult, color=Temp)) +
  #rory_theme +
  facet_wrap(~sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Femur length (mm)") + xlab("Site")

#Pronotum_Adult length
ggplot(data = Instars, aes(x = Site, y = Pronotum_Adult, color=Temp)) +
  rory_theme +
  facet_wrap(~sex+Light, nrow = 1) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1.5, position = position_dodge(width = .6)) +
  stat_summary(aes(x = as.numeric(as.factor(Site))), fun.y = mean, geom = "line", size = 1.5, position = position_dodge(width = .6)) +
  #stat_summary(fun.data = Give.N, geom = "text", fun.y = median, position = position_dodge(width = .6), size = 7) +
  ylab("Femur length (mm)") + xlab("Site")

#------
#DROP PC ANALYSIS?
#PC analysis
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
PC1_se <- unlist(tapply(PC.scores1[, "PC1"], TS, std.error))
PC2_se <- unlist(tapply(PC.scores1[, "PC2"], TS, std.error))
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

#=================
#By instar

#Development
Instars_Adult <- Instars[is.na(Instars$Age_Adult) == FALSE,]
Instars_Adult_Long <- melt(Instars_Adult[,c(1:8,50:54)], id.vars = c("Individual", "CCode", "Female", "Female2", "Site", "Species", "sex", "Temp", "Light")) 

#Make instar numeric
Instars_Adult_Long$instar=3
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Age_4th"]=4
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Age_5th"]=5
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Age_Adult"]=6

Lmer0.1 <- lmer(as.numeric(value) ~ sex * Site * Temp * Light * instar + 
                  (1|Female2/Individual), na.action = 'na.omit', REML=FALSE, data = Instars_Adult_Long)

summary(Lmer0.1)
Anova(Lmer0.1)

plot_model(Lmer0.1, type="pred", terms=c("instar","Temp","Light","Site","sex"), show.data=FALSE)

#split male and female?
#split sex
Lmer0.1 <- lmer(as.numeric(value) ~ Site * Temp * Light * instar + 
                  (1|Female2/Individual), na.action = 'na.omit', REML=FALSE, data = Instars_Adult_Long[Instars_Adult_Long$sex=="F",])

plot_model(Lmer0.1, type="pred", terms=c("instar","Temp","Light","Site"), show.data=FALSE)

#---------------------
#Mass
Instars_Adult <- Instars[is.na(Instars$Age_Adult) == FALSE,]
Instars_Adult_Long <- melt(Instars_Adult[,c(1:8,18,21,24,30,54)], id.vars = c("Individual", "CCode", "Female", "Female2", "Site", "Species", "sex", "Temp", "Light")) 

#Make instar numeric
Instars_Adult_Long$instar=3
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Mass_4th"]=4
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Mass_5th"]=5
Instars_Adult_Long$instar[Instars_Adult_Long$variable=="Mass_Adult"]=6

Lmer0.1 <- lmer(as.numeric(value) ~ sex * Site * Temp * Light * instar + 
                  (1|Female2/Individual), na.action = 'na.omit', REML=FALSE, data = Instars_Adult_Long)

summary(Lmer0.1)
Anova(Lmer0.1)

plot_model(Lmer0.1, type="pred", terms=c("instar","Temp","Light","Site","sex"), show.data=FALSE)

#split male and female?
#split sex
Lmer0.1 <- lmer(as.numeric(value) ~ Site * Temp * Light * instar + 
                  (1|Female2/Individual), na.action = 'na.omit', REML=FALSE, data = Instars_Adult_Long[Instars_Adult_Long$sex=="M",])

plot_model(Lmer0.1, type="pred", terms=c("instar","Temp","Light","Site"), show.data=FALSE)



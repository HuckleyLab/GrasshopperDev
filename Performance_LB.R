#load Rory's Hopping.R for data

#hopping performance analyses
ggplot(data = Hop, aes(x = TestTemp, y = Mean, col = Site, linetype = Sex)) +
  rory_theme +
  facet_grid(Light~Temp) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1, position = position_dodge(width = 3)) +
  stat_summary(fun.y = mean, geom = "line", size = 2, position = position_dodge(width = 3)) +
  ylab("Mean hopping distance (cm)") + xlab("Test temperature") +
  scale_x_continuous(limits = c(5, 40), breaks = c(10, 17, 25, 35))

#plot based on individual medians
ggplot(data = Hop, aes(x = TestTemp, y = Median, col = Site, linetype = Sex)) +
  rory_theme +
  facet_grid(Light~Temp) +
  stat_summary(fun.y = mean, fun.ymin = MinSE, fun.ymax = MaxSE, na.rm = TRUE, size = 1) +
  stat_summary(fun.y = mean, geom = "line", size = 2) +
  ylab("Median hopping distance (cm)") + xlab("Test temperature") +
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

#==================================
#FEEDING
#Load Rory's Feeding.R

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




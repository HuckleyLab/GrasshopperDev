#DEVELOPMENT PLOT

#try ordered factor
#Instars$Site= factor(Instars$Site, ordered=FALSE, levels=c("A1","B1","C1") )

nLmer0 <- lmer(as.numeric(Age_Adult) ~ sex * Site * Temp * Light + 
                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars)

#split sexes
#nLmer0 <- lmer(as.numeric(Age_Adult) ~ Site * Temp * Light + 
#                 (1|Female2), na.action = 'na.omit', REML=FALSE, data = Instars[Instars$sex=="F",])

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



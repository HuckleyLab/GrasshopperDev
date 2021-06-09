#Abbreviate to focus on summer GDDs and quantiles

#load libraries
library(ggplot2)
library(plyr) 
library(dplyr)
library(colorRamps)
library(patchwork)
library(viridis)
library(tidyverse)

sites= c("CHA", "A1", "B1", "C1", "D1")  #Redfox: 1574
elevs= c(1752, 2195, 2591, 3048, 3739)

#source degree days function
setwd("/Users/lbuckley/HopperPhenology/")
source("degreedays.R")

#======================================================
#READ DATA

#fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
#clim= read.csv("AlexanderClimateAll_filled_May2018.csv")
clim2= read.csv("AlexanderClimateAll_filled_Oct2019.csv")

#drop NOAA
clim2= clim2[-which(clim2$Site=="NOAA"),]

#----

#calculate degree days
inds= which(!is.na(clim2$Min) & !is.na(clim2$Max))
clim2$dd=NA
clim2$dd[inds]= mapply(degree_days, T_min=clim2$Min[inds], T_max=clim2$Max[inds], LDT=0, UDT=100, method="single.sine")

#cummulative dd by year, starting March 1
clim2= clim2%>%
  group_by(Year)%>% arrange(Ordinal) %>%
  dplyr::mutate(cumsum=cumsum(replace_na(dd, 0)))

#---------------------
#cummulative degree days
#cumsum within groups
clim1 = clim2 %>% group_by(Year,Site) %>% arrange(Ordinal) %>% dplyr::mutate(cdd_sum = cumsum(replace_na(dd, 0)) ) 

#load hopper data
setwd( paste(fdir, "grasshoppers/SexCombined/", sep="") )
hop= read.csv("HopperData_May2018.csv")

#======================================================
#CALCULATE GDD METRICS

#subset years to survey
clim1= clim1[which(clim1$Year %in% c(1958, 1959, 1960, 2006:2015))  ,]

#set temp outside summer to NA
clim1$Mean= (clim1$Min + clim1$Max)/2
clim1$Mean_summer= clim1$Mean
clim1[clim1$Ordinal<60 | clim1$Ordinal>243,"Mean_summer"]=NA

#metrics across years
clim1= ddply(clim1, c("Site", "Year"), summarise,
             Mean = mean(Mean_summer, na.rm=TRUE), Cdd_seas = max(cdd_sum, na.rm=FALSE) )
#    Mean = mean(Mean_summer, na.rm=TRUE), Sd = sd(Mean_summer, na.rm=TRUE),Cdd_seas = max(cdd_sum, na.rm=TRUE),Cdd_june = max(cdd_june, na.rm=TRUE),Cdd_july = max(cdd_july, na.rm=TRUE),Cdd_aug = max(cdd_aug, na.rm=TRUE),Cdd_early = max(cdd_early, na.rm=TRUE),Cdd_mid = max(cdd_mid, na.rm=TRUE),Cdd_ac = max(cdd_ac, na.rm=TRUE),Cdd_mb = max(cdd_mb, na.rm=TRUE),Cdd_ms = max(cdd_ms, na.rm=TRUE), Sem = sd(Mean_summer, na.rm=TRUE)/sqrt(length(na.omit(Mean_summer))))

#---------
# Jadult and GDDadult ~year by sites

#subset to dates with adults
hop1= hop[which(hop$in6>0),]

#subset to focal species
#specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")
specs= c("Melanoplus boulderensis")

hop1= hop1[which(hop1$species %in% specs ),]

#trim columns
hop4= hop1[,c("ordinal","species","in6","in5","in4","in3","in2","in1","total","site","period","year","sjy","cdd_sum")]

#-----------------
# match phenology to temp and dd
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
hop4$siteyear= paste(hop4$site, hop4$year, sep="")

hop4$Tmean<-NA
match1= match(hop4$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop4$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop4$cdd_seas<-NA
hop4$cdd_seas[matched]<- clim1$Cdd_seas[match1[matched]]  

#---------------
#add elevation
hop4$elevation= as.factor(elevs[match(hop4$site, sites)])

#---------------
#set up clim var

#GDDS
#add elevation
clim1$elevation= elevs[match(clim1$Site, sites)]
#add period
clim1$period="resurvey"
clim1$period[which(clim1$Year<2000)]<-"initial"

#order varaibles
#period
clim1$period= factor(clim1$period, levels=c("resurvey", "initial") )

#Order by average across sites with complete gdd data
clim.ave= subset(clim1, clim1$Site %in% c("B1","C1"))
clim.ave= aggregate(clim.ave, list(clim.ave$Year),FUN=mean )
clim1$Cdd_siteave= clim.ave$Cdd_seas[match(clim1$Year, clim.ave$Year)]
clim1$Cdd_july_siteave= clim.ave$Cdd_july[match(clim1$Year, clim.ave$Year)]

#------------------------------------------------
#CALCULATE DEVELOMENT INDEX

dat=hop4[,c("ordinal","species","in6","in5","in4","in3","in2","in1","total","year","site","period","sjy","dd","dd_sum","cdd_sum","cdd_sumfall")]

#replace instar NAs with zero
dat[which(is.na(dat$in1)), "in1"]=0 
dat[which(is.na(dat$in2)), "in2"]=0
dat[which(is.na(dat$in3)), "in3"]=0
dat[which(is.na(dat$in4)), "in4"]=0
dat[which(is.na(dat$in5)), "in5"]=0
dat[which(is.na(dat$in6)), "in6"]=0

#Calculate development index
dat$DI=0
inds=which(dat$total>0)  
dat$DI[inds]= (dat$in1[inds] +dat$in2[inds]*2 +dat$in3[inds]*3 +dat$in4[inds]*4 +dat$in5[inds]*5 +dat$in6[inds]*6)/dat$total[inds]

#----------
#code period
dat$per=1
dat$per[dat$year>2000]=2

# #code early and late species
# dat$early_late=2
# dat$early_late[dat$species %in% specs[c(2,4)]]=1

#add elevation
dat$elev= as.factor(elevs[match(dat$site, sites)])

#add seasonal GDDs
dat$siteyear= paste(dat$site, dat$year, sep="")

dat$Tmean<-NA
match1= match(dat$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
dat$Tmean[matched]<- clim1$Mean[match1[matched]]  

dat$cdd_seas<-NA
dat$cdd_seas[matched]<- clim1$Cdd_seas[match1[matched]]  

dat$Cdd_siteave<-NA
dat$Cdd_siteave[matched]<- clim1$Cdd_siteave[match1[matched]]  

#clean up
dat$year= as.factor(dat$year)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#order variables
#period
dat$period= factor(dat$period, levels=c("resurvey", "initial") )
#elevation
dat$elev= factor(dat$elev, levels=c(3048,2591,2195,1752) )
#species
dat$species= factor(dat$species) #, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida", "Melanoplus sanguinipes", "Melanoplus dawsoni") )

#------------------------------------------------
#ESTIMATE ADULTHOOD BASED ON DI

## failed using broom
#library(broom)
#broom::augment(x=fm1, newdata = Data, type.predict = "response")

#Indexed calculation
dat$spsiteyear= paste(dat$siteyear, dat$species, sep="")
combs= unique(dat$spsiteyear)

#days to predict over
doys= 150:265

#make matrix to store output
dout= data.frame(spsiteyear=combs, doy_adult= rep(NA, length(combs)),gdd_adult= rep(NA, length(combs)) ) 

for(k in 1:length(combs)){
  dats= subset(dat, dat$spsiteyear==combs[k])
  
  #require at least 5 data points
  if(nrow(dats)>=5) { 
    #doy
    doys= seq(min(dats$ordinal), max(dats$ordinal+7),5)
    
    spl<- smooth.spline(x=dats$ordinal, y=dats$DI)
    pred.spl<- predict(spl, doys)
    #extract point where almost all adults DI>5.5
    dout[k,2]= doys[which.max(pred.spl$y>5.5)]
    
    #gdd
    #restrict to observed gdds
    gdds= seq(min(dats$cdd_sumfall), max(dats$cdd_sumfall+50),10)
    
    spl<- smooth.spline(x=dats$cdd_sumfall, y=dats$DI)
    pred.spl<- predict(spl, gdds)
    #extract point where almost all adults DI>5.5
    dout[k,3]= gdds[which.max(pred.spl$y>5.5)]
    
  } #end check length
  
} #end combs

#add estimate back to df
dat$doy_adult= dout[match(dat$spsiteyear, dout$spsiteyear),"doy_adult"]
dat$gdd_adult= dout[match(dat$spsiteyear, dout$spsiteyear),"gdd_adult"]

#====================================
## FIGURE

##read processed data
#setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/GrasshopperDev2015/data/")
#dat= read.csv("MB_phen_clim.csv")

#DEVELOPMENTAL INDEX
#Plot DI by ordinal date

#update elevation labels
dat$elev.lab= paste(dat$elev,"m",sep="")
dat$elev.lab= factor(dat$elev.lab, levels=c("3048m","2591m","2195m","1752m") )

di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=Cdd_siteave, group=siteyear, linetype=period))+facet_grid(elev.lab~.) +
  theme_bw()+xlim(120,200)+
  geom_point()+geom_line(aes(alpha=0.5))+ #+geom_smooth(se=FALSE, aes(alpha=0.5), span=2)+
  scale_colour_gradientn(colours =matlab.like(10))+ylab("development index")+xlab("day of year")+labs(color="mean season gdds")+
  theme(legend.position = "none") + guides(alpha=FALSE)+
  theme(axis.title=element_text(size=12))

#Plot DI by GDD
di.plot.gdd= ggplot(data=dat, aes(x=cdd_sum, y = DI, color=Cdd_siteave, group=siteyear, linetype=period))+facet_grid(elev.lab~.) +
  theme_bw()+xlim(0,200)+
  geom_point()+geom_line(aes(alpha=0.5))+ #+geom_smooth(se=FALSE, aes(alpha=0.5),span=2)+
  scale_colour_gradientn(colours =matlab.like(10))+ylab("development index")+xlab("cummulative growing degree days")+labs(color="mean season gdds")+
  theme(legend.position = "right") + guides(alpha=FALSE)+
  theme(axis.title=element_text(size=12))

#stats
mod1 <- lmer(DI~poly(ordinal,2)*Cdd_siteave*elev.lab +
                 (1|year), na.action = 'na.omit', REML=FALSE, data = dat)
mod1 <- lmer(DI~poly(cdd_sum,2)*Cdd_siteave*elev.lab +
               (1|year), na.action = 'na.omit', REML=FALSE, data = dat)

Anova(mod1, type=3)

#====================================
## M. boulderensis
# PLOT ADULT PHENOLOGY
# estimated by DI

#aggregate to spsiteyear
dat.ssy= dat[,c("elev","species","cdd_seas","doy_adult","gdd_adult","spsiteyear","elev.lab","period","year","Cdd_siteave")  ]
dups= duplicated(dat.ssy$spsiteyear)
dat.ssy=dat.ssy[which(dups==FALSE),]

dat.ssy$elevspec= paste(dat.ssy$elev, dat.ssy$species, sep="")
elevspec= matrix(unique(dat.ssy$elevspec))

#------

#DOY
plot.phen.doye=ggplot(data=dat.ssy, aes(x=Cdd_siteave, y = doy_adult, color=elev.lab))+
  geom_point(aes(shape=period, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, stroke=1), size=3)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+ylab("day of year")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))+theme(legend.position="none")+scale_colour_viridis_d()+
  theme(axis.title=element_text(size=12))
#GDD metrics: Cdd_siteave cdd_seas

#GDD
plot.phen.gdde=ggplot(data=dat.ssy, aes(x=Cdd_siteave, y = gdd_adult, color=elev.lab))+
  geom_point(aes(shape=period, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, stroke=1), size=3)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+ylab("cummulative gdds")+xlab("season growing degree days (C)")+
  labs(color="elevation (m)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))+scale_colour_viridis_d()+
  theme(axis.title=element_text(size=12))

#stats
mod1= lm(doy_adult~ Cdd_siteave*elev, data=dat.ssy)
mod1= lm(gdd_adult~ Cdd_siteave*elev, data=dat.ssy)
anova(mod1)

#single elevation
mod1= lm(gdd_adult~ Cdd_siteave, data=dat.ssy[dat.ssy$elev=="3048",])

#FIG?
setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/GrasshopperDev2015/figures/")

pdf("Fig4_phen_est.pdf", height = 8, width = 8)
(di.plot +di.plot.gdd)/(plot.phen.doye + plot.phen.gdde) +  plot_layout(heights = c(1, 0.4))+ 
  plot_annotation(tag_levels = 'A')
dev.off()

#=================================
#by elevation
ggplot(data=dat, aes(x=elev, y = doy_adult, color=factor(Cdd_siteave)))+
  geom_point() +geom_smooth(method="lm",se=F) #+geom_line() 

ggplot(data=dat, aes(x=elev, y = gdd_adult, color=factor(Cdd_siteave)))+
  geom_point() +geom_smooth(method="lm",se=F)

#Cdd_siteave





# Herbivore driven decline of invasive macroalgae in Kāne‘ohe Bay, Hawai‘i #
# Code & analysis by Morgan Winston & Mary Donovan #

# This code reads in four datasets:
# 1) Benthic transect level data with paired environmental and fish data - fish data in this dataset was summarized per month to faciliate merging with benthic data
# 2) Reef fish transect data 
# 3) Daily rain data for appendix
# 4) Daily sea surface temperature data for appendix
# The code formats the data for input into a GAM model and exports figures and model outputs. 

#### INITIALIZATION ####
## load required packages
library(plyr)
library(reshape2)
library(knitr)
library(kableExtra)
library(dplyr) # yes
library(mgcv) # yes
library(plotrix)
library(vegan) 
library(tidyr) # yes
library(ggplot2) #yes
library(scales) # yes

# set working directory
setwd()

# import & format data
raw <- read.csv("all_kbay_combined.csv") 
fishcombo_all <- read.csv("fishcombo_all.csv") 
sst_dat <- read.csv("sst_dat.csv")
rain_dat <- read.csv("rain_dat.csv")

# create an index for time
temp <- NA
for(i in 0:(length(unique(raw$Year))-1)){
  t <- rep(2011+i,12) # create 12 reps of each year
  temp <- c(temp,t)
}
temp <- temp[2:length(temp)] # remove the initial NA value
fulltime <- data.frame(Month=rep(NA,length(temp)),Year=rep(NA,length(temp))) # create a placeholder df with month and year columns
fulltime$Year <- temp # fill in reps of year into this column
fulltime$Month <- rep(seq(1:12),length(unique(raw$Year))) # add in months (1-12) for each year
fulltime <- fulltime[11:nrow(fulltime),] # remove months prior to 11/2011 from data frame (no surveys were done)
fulltime$time.ind <- seq(1:nrow(fulltime)) # add index for each month/year timepoint
raw <- left_join(raw,fulltime,by=c("Month","Year")) # add index to raw data
raw$dum <- 1 # dummy var

# Figure A1: survey effort over time across reefs ---------
surveys <- aggregate(list(n = raw$Transect), by = list(date = raw$SurveyDate, Reef = raw$Reef, Method = raw$Method, Treatment = raw$Treatment), function(x) length(unique(x)))
surveys$Reef <- as.factor(surveys$Reef)
surveys$Method <- as.factor(surveys$Method)
surveys$date <- as.Date(surveys$date)
surveys$Treatment <- as.factor(surveys$Treatment)

surveys$t_date <- NA
surveys$t_date <- as.Date(surveys$t_date)
for(i in c(1:nrow(surveys))){
  if(surveys$Reef[i] == 9){
    surveys$t_date[i] <- "2019-05-01"
  }
  if(surveys$Reef[i] == 23){
    surveys$t_date[i] <- "2019-03-01"
  }
  if(surveys$Reef[i] == 28){
    surveys$t_date[i] <- "2017-09-01"
  }
  if(surveys$Reef[i] == 10){
    surveys$t_date[i] <- "2014-10-01"
  }
  if(surveys$Reef[i] == 14){
    surveys$t_date[i] <- "2015-10-01"
  }
  if(surveys$Reef[i] == 16){
    surveys$t_date[i] <- "2015-07-01"
  }
  if(surveys$Reef[i] == 19){
    surveys$t_date[i] <- "2014-10-01"
  }
  if(surveys$Reef[i] == 26){
    surveys$t_date[i] <- "2011-01-01"
  }
  if(surveys$Reef[i] == 27){
    surveys$t_date[i] <- "2012-03-01"
  }
  if(surveys$Reef[i] == 29){
    surveys$t_date[i] <- "2012-08-01"
  }
}

surveys$date <- as.POSIXct(surveys$date, "%Y-%m-%d")
surveys$t_date <- as.POSIXct(surveys$t_date, "%Y-%m-%d")
surveys$Method <- as.character(surveys$Method)
surveys[which(surveys$Method == "photoquad"),]$Method <- "Photoquadrat"
surveys$Method <- as.factor(surveys$Method)

# setwd("") set working directory for output
tiff(file='Figure A1.tif',height=6000,width=4800,res=300, compression = "lzw")
ggplot(surveys, aes(x=date, y=n, color = Method)) + 
  geom_jitter(size = 3) +
  geom_vline(data = surveys, mapping = aes(xintercept = t_date), show.legend = F) +
  scale_linetype(name = " ") +
  facet_wrap(~ Reef, ncol = 2) +
  xlab("") +
  ylab("Transects Conducted (#) \n") +
  scale_x_datetime(limits = c(min(surveys[which(!is.na(surveys$t_date)),]$t_date), 
                              max(surveys$date)), 
                   date_breaks = "12 months", labels = date_format("%m/%Y")) +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "White"),
    legend.position = "bottom",
    text = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

# Figure 3: % cover of E/K over time ---------
temp <- raw %>% group_by(Habitat,Reef,Method,time.ind) %>%summarise('kapp'=mean(kapp,na.rm=T)) # summarize kapp per time.ind/method/reef/habitat
temp$dum <- 1
temp <- as.data.frame(temp)
temp$Habitat <- as.factor(temp$Habitat)
temp$Method <- as.factor(temp$Method)

tiff(file='Figure3.tif',height=3800,width=4800,res=300, compression = "lzw")
ggplot(temp, aes(x = time.ind, y = kapp)) +
  geom_point(size=2) +
  ylab("E/K Percent Cover (%)\n") +
  xlab("") + 
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "White"),
    text = element_text(size = 30)
  ) +
  stat_smooth(method = "loess") +
  scale_x_continuous(breaks=c(1,10,20,30,40,50,60,70,80,90), labels=format(timego$mo_yr[timego$time.ind %in% c(1,10,20,30,40,50,60,70,80,90)],'%m/%y'))
dev.off()


# Figure A2: % cover of E/K per reef over time ---------
temp <- raw
raw$MonthYr <- format(as.Date(raw$SurveyDate), "%Y-%m")
reef_sum <- data.frame(aggregate(list(kapp = raw$kapp), by = list(Reef = raw$Reef, date = raw$MonthYr), mean))
reef_sum2 <- data.frame(aggregate(list(kapp_se = raw$kapp), by = list(Reef = raw$Reef, date = raw$MonthYr), std.error))
reefs <- merge(reef_sum, reef_sum2)

reefs$Reef <- as.factor(reefs$Reef)
reefs$date <- format(as.Date(paste(reefs$date, "-01", sep = "")), "%Y-%m-%d")
reefs$date <- as.Date(reefs$date)

tiff(file='Figure A2.tif',height=7000,width=4800,res=300, compression = "lzw")
ggplot(reefs, aes(x = date, y = kapp, color = Reef)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=kapp-kapp_se, ymax=kapp+kapp_se), width=.2,
                position=position_dodge(0.05)) +
  geom_line(group = 1) + 
  facet_wrap(~ Reef, scales = 'free_x', ncol = 2) +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "White"),
    legend.position = "none",
    text = element_text(size = 22),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("E/K Percent Cover (%)") +
  scale_x_date(breaks = "6 months", labels = date_format(format = "%b %Y"))
dev.off()

#### ANALYSIS ####
fvar <- raw[,c('kapp','Habitat','Method','Reef','time.ind','SurveyDate',"Treatment", 
               "H_bio","Hbrow","Hscex","H_dens","totBio","sst","wind","rain","par",'dum')] # subset to what will go in the model
fvar <- fvar[!is.na(fvar$H_bio),] # remove NA vals for herb biomass -- there are quite a few from when fish surveys were not conducted at the same month as benthic surveys
fvar.s <- fvar # create copy
fvar.s$rain <- log(fvar.s$rain) # removes outlier
for(i in c("H_bio","Hbrow","Hscex","H_dens","totBio","sst","wind","rain","par")) fvar.s[i] <- scale(fvar.s[i])[,1] # scale each predictor
fvar.s$Treatment <- as.factor(fvar.s$Treatment) # make random effects factors for gam to run properly
fvar.s$Habitat <- as.factor(fvar.s$Habitat)
fvar.s$Method <- as.factor(fvar.s$Method)
fvar.s$Reef <- as.factor(fvar.s$Reef)

#### check for correlations #
# full_pred <- fvar.s[,c("H_bio",
#                        "sst",
#                        "wind",
#                        "rain",
#                        "par")]
# pred_cor <- cor(full_pred)
# head(round(pred_cor,2))
# corrplot::corrplot(pred_cor, method="circle") # looks pretty good (why? no coefficients > 0.6 or < -0.6)
# pairs(full_pred)

fvar.s$kapp_p <- fvar.s$kapp/100 # will run model on proportion data

# run GAMM: reef included in model
kapp.gam.all <- gam((sqrt(fvar.s$kapp_p)) ~  s(time.ind, k=4) + # note square-root transformation
                      s(Habitat, bs='re', by = dum) +
                      s(Method, bs='re', by = dum) +
                      s(Reef, bs='re', by = dum) +
                      s(Treatment, bs='re', by = dum) +
                      s(H_bio, k=4) +
                      s(sst, k=4) +
                      s(wind, k=4) +
                      s(rain, k=4) +
                      s(par, k=4),
                    family = 'gaussian',
                    data = fvar.s,
                    correlation=corAR1(~time.ind),
                    method = "REML",
                    na.action='na.fail')

# model checks to run:
# plot(kapp.gam.all) 
# par(mfrow = c(2, 2))
# gam.check(kapp.gam.all) 

# Figure 4: GAMM model results when reef included in model  ---------
tiff(file='Figure4.tif',height=1400,width=1100,res=300, compression = "lzw")
plot(kapp.gam.all,select=c(6),shade=T,xlab='Herbivore Biomass',ylab='s(Herbivore Biomass)',cex.lab=1.3,cex.axis=1.3)
dev.off()

# run GAMM: reef NOT included in model
kapp.gam.all <- gam((sqrt(fvar.s$kapp_p)) ~  s(time.ind, k=4) + # note square-root transformation
                      s(Habitat, bs='re', by = dum) +
                      s(Method, bs='re', by = dum) +
                      #s(Reef, bs='re', by = dum) +
                      s(Treatment, bs='re', by = dum) +
                      s(H_bio, k=4) +
                      s(sst, k=4) +
                      s(wind, k=4) +
                      s(rain, k=4) +
                      s(par, k=4),
                    family = 'gaussian',
                    data = fvar.s,
                    correlation=corAR1(~time.ind),
                    method = "REML",
                    na.action='na.fail')

summary(kapp.gam.all)

# Figure 5: GAMM model results when reef NOT included in model  ---------
tiff(file='Figure5.tif',height=1000,width=2400,res=300, compression = "lzw")
par(mfrow=c(1,3),mar=c(3,3,2,1),mgp=c(1.7,.6,0)) ## REEF NOT INCLUDED MODEL
plot(kapp.gam.all,select=c(5),shade=T,xlab='Herbivore Biomass',ylab='s(Herbivore Biomass)',cex.lab=1.3,cex.axis=1.3)
plot(kapp.gam.all,select=c(6),shade=T,xlab='SST',ylab='s(SST)',cex.lab=1.3,cex.axis=1.3)
plot(kapp.gam.all,select=c(8),shade=T,xlab='Rainfall',ylab='s(Rainfall)',cex.lab=1.3,cex.axis=1.3)
dev.off()


# investigate herbivore biomass over time - uses fishcombo data
temp <- left_join(fishcombo_all, fulltime, by=c('Month','Year')) 
temp <- temp[!is.na(temp$time.ind),]
temp$dum <- 1
str(temp)
summary(temp)
temp$Reef <- as.factor(temp$Reef)

# run herbivore biomass GAMM
H.gam <- gam(log(H_bio+1)~s(time.ind,k=5)+
                          s(Reef,bs='re',by=dum),
                          family='gaussian',
                          data=temp,
                          method = "REML",
                          correlation=corAR1(~time.ind))
summary(H.gam)
newdat <- data.frame(time.ind=temp$time.ind, Reef=temp$Reef, dum = 0)
pred.H <- predict.gam(H.gam,se.fit=T,newdata=newdat)
fit.H <- exp(pred.H$fit)-1
up.H <- exp(pred.H$fit+1.96*pred.H$se.fit)-1
down.H <- exp(pred.H$fit-1.96*pred.H$se.fit)-1
pred.H <- cbind(fit.H,up.H,down.H,temp[c('H_bio','time.ind')])
pred.H <- pred.H[with(pred.H,order(time.ind)),]

# Figure 6: herbivore biomass over time  ---------
tiff(file='Figure6.tif',height=2000,width=2900,res=300, compression = "lzw")
plot(pred.H$H_bio~pred.H$time.ind,ylab="",xlab="",pch=21,col='grey45',bg='grey85',xaxt='n',type='n',yaxt="n")
axis(1, at=c(1,10,20,30,40,50,60,70,80,90),labels=format(timego$mo_yr[timego$time.ind %in% c(1,10,20,30,40,50,60,70,80,90)],'%m/%y'), cex.axis=1.3)
polygon(c(rev(pred.H$time.ind),pred.H$time.ind),c(rev(pred.H$up),pred.H$down),border=NA,col=rgb(0,0,255,100,max=255))
points(pred.H$fit~pred.H$time.ind,type="l",lwd=2,col=rgb(0,0,255,200,max=255))
points(pred.H$H_bio~pred.H$time.ind,pch=21,col='black',bg=rgb(0,0,255,60,max=255))
axis(2, cex.axis=1.3)
mtext(expression("Herbivore Biomass"~~bgroup("(",'g '*m^{-2},")")),side=2, line = 2, cex = 1.3)
dev.off()

# Figure 7: herbivore community composition  ---------
herbs <- raw %>% group_by(Month,Year,Reef) %>% summarise(Hgd=median(Hgd,na.rm=T),Hbrow=median(Hbrow,na.rm=T),Hscex=median(Hscex,na.rm=T),kapp=median(kapp,na.rm=T)) %>% ungroup() ## do something better for the kapp numbers?
herbs <- herbs[!is.na(herbs$Hgd),]
herbs$sum <- rowSums(herbs[4:6])
herbs <- herbs[!herbs$sum==0,]
herbs.f <- log(herbs[4:6]+1)
herbs.mds <- suppressMessages(metaMDS(herbs.f,trace=0))

tiff(file='Figure7.tif',height=1800,width=2700,res=300, compression = "lzw")
plot(herbs.mds,xlim=c(-3,2),xaxt='n',yaxt='n',ylab='',xlab='')
points(herbs.mds, cex=herbs$kapp/10, pch=19, col='darkgreen')
reg.arrow <- scores(herbs.mds,choices=1:3,display="sp",scaling=2)
arrows(0,0,reg.arrow[,1],reg.arrow[,2],length=0,lty=1,cex=3,lwd=1.5,col="black")
text(-0.2207784,0.2108962+0.1,'grazers',cex=1.5)
text(-2.2918944-.3,0.6242135+.1,'browsers',cex=1.5)
text(-0.1774290,-0.2522467-.1,'scrapers',cex=1.5)
dev.off()

# Figure A3: SST over time  ---------
tiff(file='FigureA3.tif',height=1800,width=2900,res=300, compression = "lzw")
ggplot(sst_dat[ which(sst_dat$date >= "2011-11-03" & sst_dat$date <= "2019-10-01"),], aes(x=date, y=sst)) +
  geom_line() +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "White"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  xlab("") +
  ylab("Sea Surface Temperature (?C)") +
  scale_x_datetime(date_breaks = "12 months", labels = date_format("%m/%Y"))
dev.off()

# Figure A4: rain over time  ---------
tiff(file='FigureA4.tif',height=1800,width=2900,res=300, compression = "lzw")
ggplot(rain_dat[ which(rain_dat$date >= "2011-11-03" & rain_dat$date <= "2019-10-01"),], aes(x=date, y=rain)) +
  geom_line() +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "White"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  xlab("") +
  ylab("Rainfall (mm)") +
  scale_x_datetime(date_breaks = "12 months", labels = date_format("%m/%Y"))
dev.off()

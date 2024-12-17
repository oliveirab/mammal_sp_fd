################################################################
#### Ecological Limits on Diversity 
#### Brunno Oliveira, 2014
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
################################################################

# Load packages
library(maptools)
library(rgdal)
library(ggplot2)
library(devtools)
library(ppcor)
library(sem)
library(raster)
library(ape)
library(mgcv)
library(vegan)
library(sp)
require(reshape2)
library(picante)
library(nlme)
library(lavaan)
library(rgeos)
library(piecewiseSEM)#install_github("jslefche/piecewiseSEM")
library(AICcmodavg)
library(gridExtra)
library(geepack)
library(spdep)             #http://cran.r-project.org/src/contrib/PACKAGES.html
library(ncf)                # http://onb.ent.psu.edu/onb1/
library(plspm)
library(sesem)
library(lavaan)
library(semPlot)
library(nlme)
library(foreach);library(parallel);library(doParallel)
library(car) # outlier analysis


gc()

#Set up stuff
rm(list=ls())

### PREPARING DATA FOR THE FIRST TIME 

#setwd('/home/brunno/Dropbox/Doutorado Brunno/Manuscritos/Chap1 Age and FD/Script')
setwd('C:/Users/Brunno/Documents/brunno')
load("PD&FD.RData")

rm(list=setdiff(ls(),c("metadata","Envivars")))

metadata<-cbind(metadata,Envivars)
rm(Envivars)

metadata1 <- metadata

############################################################################
#### WORKING (ALL SET)
#setwd('/home/brunno/Dropbox/Doutorado Brunno/Manuscritos/Chap1 Age and FD')
setwd('C:/Users/Brunno/Documents/brunno')
load("Analysis_FDvsPDvsRS.RData")
#save.image("Analysis_FDvsPDvsRS.RData")

############################################################################
### Load functions:

# Remove outliers function
source('remove_outliers_brunno.R')
# Multiplot
source("multiplot.R")

#### Organizing stuff

# Rename Oceanian Realm
Realms <- metadata$Realm
Realms <- gsub('Oceanina', 'Australian', Realms) # Combine Oceanian and Australian
metadata <- subset(metadata,select=-c(Realm))
metadata <- cbind(metadata,Realm=Realms)

# vars names
prednames<- c("Functional Dispersion","Richness", "AET","DivRates","Evolutionary Time","ses.MPD","ses.FD")

# vars cod
predvar<- c("FDis","SPD","aet","DR","MPD",'ses.MPD','ses.FD')

# remove outliers in FD
FDis1 <- metadata$FDis # save backup
metadata$FDis <- reout(FDis1,na.q=1,probs.u = c(.05,.95))

# dataset with vars
data<-data.frame(metadata[,predvar])
names(data)<-predvar # name columns

# Subset the Australian and Oceanian realms (They are outliers in MPD)
AUS<-subset(metadata,Realm=="Australian")
AUSxy <- cbind(x=AUS$x,y=AUS$y)
AUSxy <- as.data.frame(AUSxy)
AUSRealm<-as.vector(AUS$Realm)

AUS<-data.frame(AUS[,predvar])
names(AUS)<-c(predvar)

# The world without Australian and Oceanian realms
AUS0 <- subset(metadata, Realm=="Palearctic" | Realm=="Afrotropical" |
                 Realm=="Saharo-Arabian" | Realm=="Oriental" | Realm=="Sino-Japanese" | Realm=="Madagascan"  |
                 Realm=="Nearctic" |Realm=="Neotropical" | Realm=="Panamanian")

AUS0xy <- cbind(x=AUS0$x,y=AUS0$y)
AUS0xy <- as.data.frame(AUS0xy)
AUS0Realm<-as.vector(AUS0$Realm)
AUS0<-data.frame(AUS0[,predvar])
names(AUS0)<-predvar

save.image("Analysis_FDvsPDvsRS.RData")

############################################################################
#See histograns
par(mfrow=c(3,3))
for(i in 1:length(predvar)){  hist(metadata[predvar[i]][,1],main=prednames[i])}
dev.off()

par(mfrow=c(3,3))
for(i in 1:length(predvar)){  hist(data[predvar[i]][,1],main=prednames[i])}
dev.off()

par(mfrow=c(3,3))
for(i in 1:length(predvar)){  hist(AUS0[predvar[i]][,1],main=prednames[i])}
dev.off()

par(mfrow=c(3,3))
for(i in 1:length(predvar)){  hist(AUS[predvar[i]][,1],main=prednames[i])}
dev.off()

par(mfrow=c(1,1))
############################################################################

# Map of mundinents
mund <- readShapePoly("F:/GIS/Shp files/Mundi.shp")
lps <- getSpPPolygonsLabptSlots(mund)
IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
mund   <- unionSpatialPolygons(mund,IDOneBin)
crs(mund) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
mund <- spTransform(mund, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

# Subset New and Old World
Old<-subset(metadata,mundi=="Africa" | mundi=="Europe" | mundi=="Asia-Temperate" | mundi=="Asia-Tropical" | mundi=="Australasia")
Old$mundi <- "Old World"
New<-subset(metadata,mundi=="Southern America" | mundi=="Northern America")
New$mundi <- "New World"
NO<-rbind(Old,New)

############################################################################
### Write maps
for (i in seq_along(prednames)){
  png(paste("Figures/#Maps/",predvar[i],".png",sep=""),width=600,height=400)
  x <- cbind(metadata[,c('x','y')],data[predvar[i]])
  coordinates(x)<-~x+y
  gridded(x) <- TRUE
  x <- raster(x)
  plot(mund)
  plot(x,main=paste(prednames[i]),xaxt='n',yaxt='n',axes=F,box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
  dev.off()
}

############################################################################
### Plot maps evolutionary rates

varsplot <- metadata[,27:32]
length(varsplot)

pdf(paste("teste.pdf",sep=""),width=12,height=6,useDingbats=F)
par(mfrow = c(3, 2),oma = c(0,0,0,1),mar=c(0,0,1.5,0))
for (i in seq_along(varsplot)){
  x <- cbind(metadata[,c('x','y')],varsplot[,i])
  coordinates(x)<-~x+y
  gridded(x) <- TRUE
  x <- raster(x)
  plot(x,main=names(varsplot[i]),xaxt='n',yaxt='n',axes=F,box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
  plot(add=T,mund)
}
dev.off()

## which model is better?

models <- varsplot[,c(2,4,6)]

res <- NA
for(i in 1:nrow(models)){
  res[i] <- as.numeric(which.min(apply(varsplot[i,],MARGIN = 2,min)))
}

## Create maps manually

SR <- cbind(metadata[,c('x','y')],Var=metadata$ses.FD)
coordinates(SR)<-~x+y
gridded(SR) <- TRUE
SR <- raster(SR)
plot(SR,main='',axes=F,box=F, col=colorRampPalette(c('blue','yellow','red'))(40))
#plot(mund,add=TRUE,lwd=.5,axes=F,box=F)
dev.off()
hist(sqrt(SR),main="Mammal richness",xlab="")
boxplot(SR,main="Mammal richness",xlab="")
#Where are the outliers?
nor <- metadata[,'AGM']
out <- remove_outliers1(metadata[,'AGM'])
higher <- rep(NA,length(nor))
lower <- rep(NA,length(nor))
for (i in 1:length(nor)){higher[i]<-ifelse(nor[i]>out[i],2,0)}
for (i in 1:length(nor)){lower[i]<-ifelse(nor[i]<out[i],1,0)}
hilo <- higher+lower
temp <- cbind(cbind(metadata[,c('x','y')],Var=hilo))
coordinates(temp)<-~x+y
gridded(temp) <- TRUE
temp <- raster(temp)
plot(temp,main="Median community age outliers",xaxt='n',yaxt='n',box=F,axes=F,col=colorRampPalette(c("white", "blue", "red"))(255))
plot(mund,add=T)

### crop mundi (remove Antarctica)
SR <- cbind(metadata[,c('x','y')],Var=metadata$FDis); coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
e <- extent(SR)
mundi.crop <- crop(mund, e, snap="out", filename="mundi.crop.tif")

### Figure 1
pdf(paste("Fig1_.pdf",sep=""),width=12,height=6,useDingbats=F)
par(mfrow = c(2, 2),oma = c(0,0,0,1),mar=c(0,0,1.5,0))
#1) Rich
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$SPD)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR)),2),round(maxValue(na.omit(SR)),2));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Species richness'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#2) FDis
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$FDis)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR)),2),round(maxValue(na.omit(SR)),2));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Functional diversity'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#3.1) MPD AUS
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(AUSxy[,c('x','y')],Var=AUS$MPD)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
image(SR, add=T,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
#3.2) MPD global
SR <- cbind(AUS0xy[,c('x','y')],Var=AUS0$MPD); coordinates(SR)<-~x+y; gridded(SR) <- TRUE
SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR))+.1,1),round(maxValue(na.omit(SR))-.1,1));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Evolutionary time'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#4) DR
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$DR)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR))+.1,1),round(maxValue(na.omit(SR))-.1,1));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Diversification rate'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
dev.off()

par(mfrow=c(1,1))

### Figure 1 SES
pdf(paste("Fig1_SES.pdf",sep=""),width=12,height=6,useDingbats=F)
par(mfrow = c(2, 2),oma = c(0,0,0,1),mar=c(0,0,1.5,0))
#1) Rich
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$SPD)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR)),2),round(maxValue(na.omit(SR)),2));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Species richness'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#2) FDis
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$ses.FD)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR)),2),round(maxValue(na.omit(SR)),2));ax <- list(at=zlim, labels=zlim)
plot(SR, title('S.E.S. Functional diversity'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#3.1) MPD AUS
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(AUSxy[,c('x','y')],Var=AUS$ses.MPD)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
image(SR, add=T,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
#3.2) MPD global
SR <- cbind(AUS0xy[,c('x','y')],Var=AUS0$ses.MPD); coordinates(SR)<-~x+y; gridded(SR) <- TRUE
SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR))+.1,1),round(maxValue(na.omit(SR))-.1,1));ax <- list(at=zlim, labels=zlim)
plot(SR, title('S.E.S. Evolutionary time'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
#4) DR
SR<-raster(extent(e)); values(SR)<-NA; plot(SR,main='',axes=F,box=F)
SR <- cbind(metadata[,c('x','y')],Var=metadata$DR)
coordinates(SR)<-~x+y; gridded(SR) <- TRUE; SR <- raster(SR)
zlim=c(round(minValue(na.omit(SR))+.1,1),round(maxValue(na.omit(SR))-.1,1));ax <- list(at=zlim, labels=zlim)
plot(SR, title('Diversification rate'),add=T,axis.args=ax,legend.shrink=.4,legend.width=1.5,main='',axes=F,col=colorRampPalette(c("#3d46a9","#3d46a9","#3d46a9","#49bff0","#49bff0","#58c53f","#f2f327", "#f2f327", "#ec2236", "#ec2236"))(40))
plot(mundi.crop,add=T,lwd=.5)
dev.off()

par(mfrow=c(1,1))

############################################################################

#Gam formula
stat_smooth(method="gam",formula=y~ s(x,k=3),level=0.95,se=T,size=1.5,
            colour='black',show_guide=F) # Add regression line

#LOESS formula
stat_smooth(method="loess",formula=y~ x,level=0.95,se=T,size=1.5,
            colour='black',show_guide=F) # Add regression line

#LM formula
stat_smooth(method="lm",level=0.95,se=T,size=1.5,
            colour='black',show_guide=F) # Add regression line

#Polynomial formuma
stat_smooth(method="lm",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
            colour='black',show_guide=F) # Add regression line

# Figue2a
pdf(paste("Fig 2a1",".pdf",sep=""),width=8.333,height=5.55,useDingbats=F)
f2a<-ggplot(na.omit(metadata), aes(SPD,FDis))+
  geom_point(aes(colour=factor(Realm), fill = factor(Realm), shape = factor(Realm)),
             size=2,name="") + 
  scale_colour_manual(values=c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_fill_manual(values=c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_shape_manual(values=c(21,22,21,22,21,22,21,22,21,22,21))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F) # Add regression line
f2a
dev.off()


# Figue2a SES
pdf(paste("Fig 2a1 SES",".pdf",sep=""),width=8.333,height=5.55,useDingbats=F)
f2a.ses<-ggplot(na.omit(metadata), aes(SPD,ses.FD))+
  geom_point(aes(colour=factor(Realm), fill = factor(Realm), shape = factor(Realm)),
             size=2,name="") + 
  scale_colour_manual(values=c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_fill_manual(values=c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_shape_manual(values=c(21,22,21,22,21,22,21,22,21,22,21))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity (SES)")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  #scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F) # Add regression line
f2a.ses
dev.off()


############################################################################
### Correlations FD vs SR

summary(lm(FDis~SPD,data=metadata))
#linear model explains 20%
AICc(lm(FDis~SPD,data=metadata))
#AIC: -52557.93
summary(lm(FDis ~ poly(SPD, 2, raw=TRUE),data=metadata))
#polynomial model explains 47%
AICc(lm(FDis ~ poly(SPD, 2, raw=TRUE),data=metadata))
#AIC: -56580.18

# Does MPD explains the residuals of FD~SPD?
#Make the graph
AUS$x<-"AUS"
AUS<-cbind(AUS,Realm=AUSRealm)
res=residuals(lm(FDis~log(SPD),AUS))
AUS<-cbind(AUS,res)

AUS0<-cbind(AUS0,Realm=AUS0Realm)
AUS0$x<-"AUS0"
res=residuals(lm(FDis~log(SPD),AUS0))
AUS0<-cbind(AUS0,res=res)

#xx<-rbind(na.omit(AUS0),na.omit(AUS))

#xx$Realm<-factor(xx$Realm,levels=sort(unique(levels(xx$Realm))))#Reorder Realms

#Figure 2b
# Lets do independent plots for AUS0 and AUS(it goes to the Supp Docs) . If one wnats to plot both AUS0 and AUS together use xx as dataset and use - aes(MPD,res,group=factor(AUS)) - in stat_smooth
# for AUS0
pdf(paste("Figure 2b_Global",".pdf",sep=""),width=8.33,height=5.55,useDingbats=F)
f2b1<-ggplot(na.omit(AUS0), aes(MPD,res))+
  geom_point(aes(colour=factor(Realm), fill = factor(Realm), shape = factor(Realm)),
             size=2,name="") + 
  scale_colour_manual(values=c("#F8766D","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_fill_manual(values=c("#F8766D","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3", "#FF62BC"))+
  scale_shape_manual(values=c(21,21,22,21,22,21,22,21,22,21))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  #geom_abline(intercept = -0.6180711, slope = 0.003807)+
  #geom_abline(intercept = -0.25327850, slope = 0.0011039921)
  stat_smooth(method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='black',show_guide=F) # Add regression line
f2b1
dev.off()

# for AUS
pdf(paste("Figure 2b_Australian",".pdf",sep=""),width=8.33,height=5.55,useDingbats=F)
f2b2<-ggplot(na.omit(AUS), aes(MPD,res))+
  geom_point(aes(colour=factor(Realm), fill = factor(Realm), shape = factor(Realm)),
             size=2,name="") + 
  scale_shape_manual(values=c(21,22,21,22,21,22,21,22,21,22,21))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  #geom_abline(intercept = -0.6180711, slope = 0.003807)+
  #geom_abline(intercept = -0.25327850, slope = 0.0011039921)
  stat_smooth(method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='black',show_guide=F) # Add regression line
f2b2
dev.off()

### Plot Figure 2 (3 figures in 1 plot)
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(f2a)

pdf(paste("Figure 2",".pdf",sep=""),width=5.28*2,height=5.55*3,useDingbats=F)
f2 <- grid.arrange(arrangeGrob(f2a + theme(legend.position="none"),
                               f2b1 + theme(legend.position="none"),
                               f2b2 + theme(legend.position="none"),
                               ncol =1),
                   mylegend, ncol=2,heights=c(10, 1))
dev.off()

## teste the correlation
modelAUS0 <- lm(res~MPD,data=AUS0)
summary(modelAUS0) # R2 global = 0.22
modelAUS <- lm(res~MPD,data=AUS)
summary(modelAUS) # R2 Australia = 0.34

## Do the residuals correlate with divrates and aet?
modelAUS0 <- lm(res~DR,data=AUS0)
summary(modelAUS0) # R2 global = 0.01
modelAUS <- lm(res~DR,data=AUS)
summary(modelAUS) # R2 Australia = 0.01

modelAUS0 <- lm(res~aet,data=AUS0)
summary(modelAUS0) # R2 global = 0.02
modelAUS <- lm(res~aet,data=AUS)
summary(modelAUS) # R2 Australia = 0.02



##########################################################################

#######################
# Spatial models

#######################
#Make a matrix of coordinates (X and Y coordinates)
dataset <- na.omit(cbind(AUS0xy, data.frame(AUS0[,predvar]),Realm = AUS0Realm))
coords <- as.matrix(dataset[,c("x","y")])

#######################
# First run OLS to use as control analysis
# OLS models

# Predicting species richness
SPD.ols<-lm(SPD ~ aet + MPD + DR ,data=dataset)
summary(SPD.ols)

# Predicting species richness
FD.ols<-lm(FDis ~ aet + MPD + DR ,data=dataset)
summary(FD.ols)

#######################
# OLS models SES

# Predicting species richness
SPD.ols.ses<-lm(SPD ~ aet + ses.MPD + DR ,data=dataset)
summary(SPD.ols)

# Predicting species richness
FD.ols.ses<-lm(ses.FD ~ aet + ses.MPD + DR ,data=dataset)
summary(FD.ols)

#######################
# GLS models 

model.runs <- list(as.formula(SPD ~ aet + DR + MPD),
                   as.formula(FDis ~ aet + DR + MPD),
                   as.formula(MPD ~ aet),
                   as.formula(DR ~ aet + MPD))
# set the cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

gls.models <- foreach(i=1:length(model.runs), .packages='nlme') %dopar% {
  gls(model.runs[[i]], data = dataset, corr = corSpatial(form = ~ x + y))
}

# stop the cluster
stopCluster(cl)

save.image("Analysis_FDvsPDvsRS.RData")

# Function to calculate standardized beta coefficient (base on the code from: http://www.r-bloggers.com/example-8-14-generating-standardized-regression-coefficients/)
std.beta <- 
  function (MOD,DATA,Y) 
  {
    b <- summary(MOD)$tTable[-1, 1]
    sx <- apply(DATA[,names(b)],2,function(x) sd(x))
    sy <- sd(DATA[,Y])
    beta <- b * sx/sy
    return(beta)
  }

std.beta(gls.models[[1]],dataset,'SPD')
std.beta(gls.models[[2]],dataset,'FDis')



#######################
# GLS models SES

model.runs <- list(as.formula(SPD ~ aet + DR + ses.MPD),
                   as.formula(ses.FD ~ aet + DR + ses.MPD),
                   as.formula(ses.MPD ~ aet),
                   as.formula(DR ~ aet + ses.MPD))
# set the cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

gls.models.ses <- foreach(i=1:length(model.runs), .packages='nlme') %dopar% {
  gls(model.runs[[i]], data = dataset, corr = corSpatial(form = ~ x + y))
}

# stop the cluster
stopCluster(cl)

save.image("Analysis_FDvsPDvsRS.RData")

# Standardized coefficients
std.beta(gls.models.ses[[1]],dataset,'SPD')
std.beta(gls.models.ses[[2]],dataset,'FDis')

#######################
# LME models
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme = lme(SPD ~ aet + MPD + DR, random = ~ 1 | Realm, data = dataset)
SPD.lme2 = lme(SPD ~ FDis + aet + MPD + DR, random = ~ 1 | Realm, data = dataset)
# Predicting functional diversity
FD.lme = lme(FDis ~ aet + MPD + DR, random = ~ 1 | Realm, data = dataset)
# Predicting diversification rates
DR.lme = lme(DR ~ aet + MPD, random = ~ 1 | Realm, data = dataset)
# MPD
MPD.lme = lme(MPD ~ aet, random = ~ 1 | Realm, data = dataset)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme, FD.lme, DR.lme , MPD.lme)

# Run goodness-of-fit tests
sem.fit(SpSEM1, dataset)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, dataset, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


#######################
# LME models SES
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme.ses = lme(SPD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset)
SPD.lme2.ses = lme(SPD ~ ses.FD + aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset)
# Predicting functional diversity
FD.lme.ses = lme(ses.FD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset)
# Predicting diversification rates
DR.lme.ses = lme(DR ~ aet + ses.MPD, random = ~ 1 | Realm, data = dataset)
# ses.MPD
MPD.lme.ses = lme(ses.MPD ~ aet, random = ~ 1 | Realm, data = dataset)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme.ses, FD.lme.ses, DR.lme.ses , MPD.lme.ses)

# Run goodness-of-fit tests
sem.fit(SpSEM1, dataset)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, dataset, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


#######################

### What is the best model?

#######################
#Define neighbourhood
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
nb1.5<-dnearneigh(coords,0,all.linked)

#Spatial weights, illustrated with coding style "W" (row standardized)
nb1.5.w<-nb2listw(nb1.5, glist=NULL, style="W", zero.policy=FALSE)


#######################
# Extract residuals and calculate Moran's I

# Predicting species richness
SPD.res.ols <- residuals(SPD.ols)
SPD.res.lme <- residuals(SPD.lme)
SPD.res.gls <- residuals(gls.models[[1]])

# Predicting species richness
FD.res.ols <- residuals(FD.ols)
FD.res.lme <- residuals(FD.lme)
FD.res.gls <- residuals(gls.models[[2]])

hist(SPD.res.ols) # residuals are normally distributed. 
hist(FD.res.ols) # residuals are normally distributed. 
hist(SPD.res.lme) # residuals are normally distributed.
hist(FD.res.lme) # residuals are normally distributed. 
hist(SPD.res.gls) # residuals are normally distributed.
hist(FD.res.gls) # residuals are normally distributed. 

#calculate Moran's I values explicitly for a certain distance,
# and to test for its significance:
require(spdep)
SPD.moran <- moran.test(SPD.res.ols, listw=nb1.5.w)
FD.moran <- moran.test(FD.res.ols, listw=nb1.5.w)
SPD.moran.lme <- moran.test(SPD.res.lme, listw=nb1.5.w)
FD.moran.lme <- moran.test(FD.res.lme, listw=nb1.5.w)
SPD.moran.gls <- moran.test(SPD.res.gls, listw=nb1.5.w)
FD.moran.gls <- moran.test(FD.res.gls, listw=nb1.5.w)


# Moran I is greater than expected, thus there is spatial autocorrelation.


#######################
#Correlograms

SPD.cor.ols <- correlog(dataset$x, dataset$y, SPD.res.ols, na.rm = T, increment = 500, resamp = 0)
SPD.cor.lme <- correlog(dataset$x, dataset$y, SPD.res.lme, na.rm = T, increment = 500, resamp = 0)
SPD.cor.gls <- correlog(dataset$x, dataset$y, SPD.res.gls, na.rm = T, increment = 500, resamp = 0)

FD.cor.ols <- correlog(dataset$x, dataset$y, FD.res.ols, na.rm = T, increment = 500, resamp = 0)
FD.cor.lme <- correlog(dataset$x, dataset$y, FD.res.lme, na.rm = T, increment = 500, resamp = 0)
FD.cor.gls <- correlog(dataset$x, dataset$y, FD.res.gls, na.rm = T, increment = 500, resamp = 0)


#Set plotting options to plot correlogram
pdf(paste("correlogram.pdf",sep=""),width=6,height=10,useDingbats=F)

#x11(width = 10, height = 6)

par(mar=c(5,5,2,0.1), mfrow=c(1,2))

distances <- seq(1,1000,50)

ylim_ = c(min(c(SPD.cor.ols$correlation[distances],SPD.cor.gls$correlation[distances],SPD.cor.lme$correlation[distances])), 
          max(c(SPD.cor.ols$correlation[distances],SPD.cor.gls$correlation[distances],SPD.cor.lme$correlation[distances])))

plot(SPD.cor.ols$correlation[distances], type="b", pch=1, cex=1.5, lwd=1.5, ylim = ylim_,
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, cex.axis=1.5)
abline(h=0)
# then lme model residuals
points(SPD.cor.lme$correlation[distances], pch=4, cex=1.2)
lines(SPD.cor.lme$correlation[distances], lwd=1.5)
# then gls model residuals
points(SPD.cor.gls$correlation[distances], pch=2, cex=1.2)
lines(SPD.cor.gls$correlation[distances], lwd=1.5)

# annotate
legend(x=15, y=0.8, legend=c("OLS"), pch=c(1), bty="n", cex=1.2)
legend(x=15, y=0.7, legend=c("LME"), pch=c(4), bty="n", cex=1.2)
legend(x=15, y=0.6, legend=c("GLS"), pch=c(2), bty="n", cex=1.2)

# FD
ylim_ = c(min(c(FD.cor.ols$correlation[distances],FD.cor.gls$correlation[distances],FD.cor.lme$correlation[distances])), 
          max(c(FD.cor.ols$correlation[distances],FD.cor.gls$correlation[distances],FD.cor.lme$correlation[distances])))

plot(FD.cor.ols$correlation[distances], type="b", pch=1, cex=1.5, lwd=1.5, ylim = ylim_,
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, cex.axis=1.5)
abline(h=0)
# then lme model residuals
points(FD.cor.lme$correlation[distances], pch=4, cex=1.2)
lines(FD.cor.lme$correlation[distances], lwd=1.5)
# then gls model residuals
points(FD.cor.gls$correlation[distances], pch=2, cex=1.2)
lines(FD.cor.gls$correlation[distances], lwd=1.5)

# annotate
legend(x=15, y=0.8, legend=c("OLS"), pch=c(1), bty="n", cex=1.2)
legend(x=15, y=0.7, legend=c("LME"), pch=c(4), bty="n", cex=1.2)
legend(x=15, y=0.6, legend=c("GLS"), pch=c(2), bty="n", cex=1.2)

dev.off()

# save matrix residuals

res.matrix <- data.frame(x=dataset$x, y=dataset$y, SPDols = SPD.res.ols, SPDlme = SPD.res.lme, SPDgls = SPD.res.gls,
                         FDols = FD.res.ols, FDlme = FD.res.lme, FDgls = FD.res.gls)
write.csv(res.matrix, "resid.matrix.csv")

################################
# What is the best model?
library(AICcmodavg)

AICc(SPD.ols)
AICc(SPD.lme)
AICc(gls.models[[1]])

AICc(FD.ols)
AICc(FD.lme)
AICc(gls.models[[2]])


################################
# Run LME again for Australia

#######################
#Make a matrix of coordinates (X and Y coordinates)
dataset2 <- na.omit(cbind(AUSxy, data.frame(AUS[,predvar]),Realm = AUSRealm))

#######################
# LME models
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme = lme(SPD ~ aet + MPD + DR, random = ~ 1 | Realm, data = dataset2)
SPD.lme2 = lme(SPD ~ FDis + aet + MPD + DR, random = ~ 1 | Realm, data = dataset2)
# Predicting functional diversity
FD.lme = lme(FDis ~ aet + MPD + DR, random = ~ 1 | Realm, data = dataset2)
# Predicting diversification rates
DR.lme = lme(DR ~ aet + MPD, random = ~ 1 | Realm, data = dataset2)
# MPD
MPD.lme = lme(MPD ~ aet, random = ~ 1 | Realm, data = dataset2)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme, FD.lme, DR.lme , MPD.lme)

# Run goodness-of-fit tests
sem.fit(SpSEM1, dataset2)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, dataset2, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


#######################
# LME models SES Australia
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme.ses = lme(SPD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset2)
SPD.lme2.ses = lme(SPD ~ ses.FD + aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset2)
# Predicting functional diversity
FD.lme.ses = lme(ses.FD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = dataset2)
# Predicting diversification rates
DR.lme.ses = lme(DR ~ aet + ses.MPD, random = ~ 1 | Realm, data = dataset2)
# ses.MPD
MPD.lme.ses = lme(ses.MPD ~ aet, random = ~ 1 | Realm, data = dataset2)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme.ses, FD.lme.ses, DR.lme.ses , MPD.lme.ses)

# Run goodness-of-fit tests
sem.fit(SpSEM1, dataset2)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, dataset2, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)



################################
# Run LME again for young

#######################
#Make a matrix of coordinates (X and Y coordinates)
young <- na.omit(metadata[metadata$MPD<quantile(metadata$MPD)[2],])

#######################
# LME models
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme = lme(SPD ~ aet + MPD + DR, random = ~ 1 | Realm, data = young)
# Predicting functional diversity
FD.lme = lme(FDis ~ aet + MPD + DR, random = ~ 1 | Realm, data = young)
# Predicting diversification rates
DR.lme = lme(DR ~ aet + MPD, random = ~ 1 | Realm, data = young)
# MPD
MPD.lme = lme(MPD ~ aet, random = ~ 1 | Realm, data = young)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme, FD.lme, DR.lme , MPD.lme)

# Run goodness-of-fit tests
sem.fit(SpSEM1, young)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, young, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


#######################
# LME models SES young
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme.ses = lme(SPD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = young)
# Predicting functional diversity
FD.lme.ses = lme(ses.FD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = young)
# Predicting diversification rates
DR.lme.ses = lme(DR ~ aet + ses.MPD, random = ~ 1 | Realm, data = young)
# ses.MPD
MPD.lme.ses = lme(ses.MPD ~ aet, random = ~ 1 | Realm, data = young)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme.ses, FD.lme.ses, DR.lme.ses , MPD.lme.ses)

# Run goodness-of-fit tests
sem.fit(SpSEM1, young)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, young, standardize = "scale")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)

################################
# Run LME again for old

#######################
#Make a matrix of coordinates (X and Y coordinates)
old <- na.omit(metadata[metadata$MPD>quantile(metadata$MPD)[4],])

#######################
# LME models
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme = lme(SPD ~ aet + MPD + DR, random = ~ 1 | Realm, data = old)
# Predicting functional diversity
FD.lme = lme(FDis ~ aet + MPD + DR, random = ~ 1 | Realm, data = old)
# Predicting diversification rates
DR.lme = lme(DR ~ aet + MPD, random = ~ 1 | Realm, data = old)
# MPD
MPD.lme = lme(MPD ~ aet, random = ~ 1 | Realm, data = old)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme, FD.lme, DR.lme , MPD.lme)

# Run goodness-of-fit tests
sem.fit(SpSEM1, old)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, old, standardize = "range")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


#######################
# LME models SES
# Lets fit the same SEM but know with only Realm as a random effect (just like Jetz and Fine 2012 did)

# Predicting species richness
SPD.lme.ses = lme(SPD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = old)
# Predicting functional diversity
FD.lme.ses = lme(ses.FD ~ aet + ses.MPD + DR, random = ~ 1 | Realm, data = old)
# Predicting diversification rates
DR.lme.ses = lme(DR ~ aet + ses.MPD, random = ~ 1 | Realm, data = old)
# ses.MPD
MPD.lme.ses = lme(ses.MPD ~ aet, random = ~ 1 | Realm, data = old)

#save.image("Script/Analysis_FDvsPDvsRS.RData")

SpSEM1 = list(  SPD.lme.ses, FD.lme.ses, DR.lme.ses , MPD.lme.ses)

# Run goodness-of-fit tests
sem.fit(SpSEM1, old)

# Obtain standardized regression coefficients
SpSEM.coefs1 = sem.coefs(SpSEM1, old, standardize = "range")

SpSEM.coefs1[, 3:5] = round(SpSEM.coefs1[, 3:5], 3)

SpSEM.coefs1[order(SpSEM.coefs1$response), ]

# Explore individual model fits
sem.model.fits(SpSEM1)


###########################################################################
# rerun the analysis without outliers

library(car)


fit <- lm(FDis ~ poly(log(SPD), 2, raw=TRUE),data=na.omit(metadata))
summary(fit)


# Assessing Outliers
### remove the botton left points, as suggested by the reviwer.

blpoints <- which(metadata$FDis<.26)

# plot relationship

newdata2<-metadata
newdata2$highlight <- NA
newdata2$highlight[blpoints] <- "outlier"
newdata2$highlight[-blpoints] <- "normal"

fit2 <- lm(FDis ~ poly(SPD, 2, raw=TRUE),data=newdata2[-blpoints,])
summary(fit2)
#polynomial model explains 42.6%

# plot line of fit without outilers
P3A1 <- ggplot(na.omit(newdata2), aes(SPD,FDis))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata2[-blpoints,],method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='red',show_guide=F)+ # Add regression line without outliers
  geom_text(aes(x = 100, y = .27, label = paste("R^2 == ", round(summary(fit2)$r.squared, 3))), parse = TRUE, size = 8) +
  annotate("text", label = paste(length(blpoints),"outliers"), size = 8, x = 90, y = .25)


#### removing extreme residuals - 5% most extreme residuals
resfit <- sort(abs(resid(fit)), decreasing = T) # absolut value of residuals in crescent order
threa <- round((length(resfit)*5)/100, 0) # how much observations corresponde to 5% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation


newdata<-na.omit(metadata)
newdata$res <- abs(resid(fit)) # absolut value of residuals
highlight <- rep(NA, nrow(newdata))
highlight[newdata$res>=threa] <- "outlier"
highlight[newdata$res<threa] <- "normal"
newdata$highlight <- highlight

fit3 <- lm(FDis ~ poly(log(SPD), 2, raw=TRUE),data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P3A2 <- ggplot(na.omit(newdata), aes(SPD,FDis))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='red',show_guide=F)+ # Add regression line without outliers
  geom_text(aes(x = 100, y = .27, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),"outliers"), size = 8, x = 90, y = .25)



#### removing extreme residuals - 10% most extreme residuals
resfit <- sort(abs(resid(fit)), decreasing = T) # absolut value of residuals in crescent order
threa <- round((length(resfit)*10)/100, 0) # how much observations corresponde to 10% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation


newdata<-na.omit(metadata)
newdata$res <- abs(resid(fit)) # absolut value of residuals
highlight <- rep(NA, nrow(newdata))
highlight[newdata$res>=threa] <- "outlier"
highlight[newdata$res<threa] <- "normal"
newdata$highlight <- highlight

fit4 <- lm(FDis ~ poly(log(SPD), 2, raw=TRUE),data=newdata[highlight=="normal",])
summary(fit4)
#polynomial model explains 53.8%

# plot line of fit without outilers
P3A3 <- ggplot(na.omit(newdata), aes(SPD,FDis))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 100, y = .27, label = paste("R^2 == ", round(summary(fit4)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),"outliers"), size = 8, x = 90, y = .25)





pdf(paste("Figure S5",".pdf",sep=""),width=5.28,height=5.55*2,useDingbats=F)
f2 <- multiplot(f2ax + theme(legend.position="none"),
                f2ax2 + theme(legend.position="none"),
                f2ax3 + theme(legend.position="none"),
                cols =1)
dev.off()



########################

# Do the same for AUS


# extract the res of SR-FD and fit the model res ~ MPD
AUS<-na.omit(AUS)

fit <- lm(res ~ MPD, AUS)
resfit <- sort(abs(resid(fit)), decreasing = T)

#### removing extreme residuals - 5% most extreme residuals
threa <- round((length(resfit)*5)/100, 0) # how much observations corresponde to 5% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation


newdata<-AUS
highlight <- rep(NA, nrow(newdata))
highlight[abs(resid(fit))>=threa] <- "outlier"
highlight[resid(fit)<threa] <- "normal"
newdata$highlight <- highlight

fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P3C1 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F)+ # Add regression line without outliers
  geom_text(aes(x = 150, y = .05, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 150, y = .03)



#### removing extreme residuals - 10% most extreme residuals
threa <- round((length(resfit)*10)/100, 0) # how much observations corresponde to 5% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation


newdata<-AUS
highlight <- rep(NA, nrow(newdata))
highlight[abs(resid(fit))>=threa] <- "outlier"
highlight[resid(fit)<threa] <- "normal"
newdata$highlight <- highlight


fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P3C2 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F)+ # Add regression line without outliers
  geom_text(aes(x = 150, y = .05, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 150, y = .03)


########################

# Do the same for AUS0

# extract the res of SR-FD and fit the model res ~ MPD
AUS0<-na.omit(AUS0)

fit <- lm(res ~ MPD, AUS0)
resfit <- sort(abs(resid(fit)), decreasing = T)


#### removing extreme residuals - 5% most extreme residuals
threa <- round((length(resfit)*5)/100, 0) # how much observations corresponde to 5% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation

newdata<-AUS0
highlight <- rep(NA, nrow(newdata))
highlight[abs(resid(fit))>=threa] <- "outlier"
highlight[resid(fit)<threa] <- "normal"
newdata$highlight <- highlight

fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P3B1 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 120, y = .1, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 120, y = .08)





#### removing extreme residuals - 10% most extreme residuals
threa <- round((length(resfit)*10)/100, 0) # how much observations corresponde to 5% of the total number of observations?
threa <- resfit[threa] # residual value at the 5% most extreme observation

newdata<-AUS0
highlight <- rep(NA, nrow(newdata))
highlight[abs(resid(fit))>=threa] <- "outlier"
highlight[resid(fit)<threa] <- "normal"
newdata$highlight <- highlight

fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P3B2 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 120, y = .1, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 120, y = .08)



pdf(paste("Figure S6",".pdf",sep=""),width=5.28*3,height=5.55*3,useDingbats=F)
f2 <- multiplot(P3A2 + theme(legend.position="none"),
                P3B1 + theme(legend.position="none"),
                P3C1 + theme(legend.position="none"),
                P3A3 + theme(legend.position="none"),
                P3B2 + theme(legend.position="none"),
                P3C2 + theme(legend.position="none"),
                P3A1 + theme(legend.position="none"),
                cols =3)
dev.off()




png(paste("Figure S6",".png",sep=""),width=528*3,height=555*3)
f2 <- multiplot(P3A2 + theme(legend.position="none"),
                P3B1 + theme(legend.position="none"),
                P3C1 + theme(legend.position="none"),
                P3A3 + theme(legend.position="none"),
                P3B2 + theme(legend.position="none"),
                P3C2 + theme(legend.position="none"),
                P3A1 + theme(legend.position="none"),
                cols =3)
dev.off()


########################
###### Alternative figure S6 - remove Studentized residuals
########################

fit <- lm(FDis ~ poly(log(SPD), 2, raw=TRUE),data=na.omit(metadata))
summary(fit)

resfit <- names(outlierTest(fit)$rstudent)


newdata<-metadata
highlight <- rep("normal", nrow(newdata))
names(highlight) <- rownames(newdata)
highlight[resfit] <- "outlier"
newdata$highlight <- highlight

fit3 <- lm(FDis ~ poly(log(SPD), 2, raw=TRUE),data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P4A <- ggplot(na.omit(newdata), aes(SPD,FDis))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  theme(text = element_text(size = 26),legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) + #encrease size of points in legend
  labs(x = "Species richness", y = "Functional diversity")+ # XY Labels
  scale_x_log10(breaks=c(10,50,200))+
  scale_y_continuous(breaks=c(0.25,0.3,0.35,0.4))+
  stat_smooth(method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="gam",formula = y ~ poly(x, 2),level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 100, y = .27, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),"outliers"), size = 8, x = 90, y = .25)



# Do the same for AUS0

# extract the res of SR-FD and fit the model res ~ MPD
AUS0<-na.omit(AUS0)

fit <- lm(res ~ MPD, AUS0)

resfit <- names(outlierTest(fit)$rstudent)


newdata<-AUS0
highlight <- rep("normal", nrow(newdata))
names(highlight) <- rownames(newdata)
highlight[resfit] <- "outlier"
newdata$highlight <- highlight

fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P4B1 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 180, y = -.07, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 180, y = -.09)


# Do the same for AUS

# extract the res of SR-FD and fit the model res ~ MPD
AUS<-na.omit(AUS)

fit <- lm(res ~ MPD, AUS)

resfit <- names(outlierTest(fit)$rstudent)


newdata<-AUS
highlight <- rep("normal", nrow(newdata))
names(highlight) <- rownames(newdata)
highlight[resfit] <- "outlier"
newdata$highlight <- highlight

fit3 <- lm(res ~ MPD, data=newdata[highlight=="normal",])
summary(fit3)
#polynomial model explains 52%

# plot line of fit without outilers
P4B2 <- ggplot(newdata, aes(MPD,res))+
  geom_point(size=2,name="", aes(colour=highlight)) + 
  scale_colour_manual(values=c("blue","red"))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))+ #encrease size of points in legend
  theme(text = element_text(size = 26),legend.title = element_blank())+
  labs(x = "Evolutionary time", y = "residuals")+ # XY Labels
  stat_smooth(method="lm",formula = y ~ x,se=T,size=1.5,
              colour='black',show_guide=F)+ # Add regression line
  stat_smooth(data = newdata[highlight=="normal",],method="lm",formula = y ~ x,level=0.95,se=T,size=1.5,
              colour='red',show_guide=F) + # Add regression line without outliers
  geom_text(aes(x = 230, y = -.08, label = paste("R^2 == ", round(summary(fit3)$r.squared, 3))), parse = TRUE, size = 8)+
  annotate("text", label = paste(length(which(newdata$highlight=="outlier")),
                                 "outliers"), size = 8, x = 230, y = -.09)




pdf(paste("Figure S2",".pdf",sep=""),width=5.28,height=5.55*3,useDingbats=F)
f2 <- multiplot(P4A + theme(legend.position="none"),
                P4B1 + theme(legend.position="none"),
                P4B2 + theme(legend.position="none"),
                cols =1)
dev.off()


png(paste("Figure S2",".png",sep=""),width=528,height=555*3)
f2 <- multiplot(P4A + theme(legend.position="none"),
                P4B1 + theme(legend.position="none"),
                P4B2 + theme(legend.position="none"),
                cols =1)
dev.off()

pdf(paste("Figure S3",".pdf",sep=""),width=5.28,height=5.55,useDingbats=F)
P3A1+ theme(legend.position="none")
dev.off()

png(paste("Figure S3",".png",sep=""),width=528,height=555)
P3A1+ theme(legend.position="none")
dev.off()

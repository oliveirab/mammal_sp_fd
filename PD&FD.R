
################################################################
#### Ecological Limits on Diversity 
#### Brunno Oliveira, 2014
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
################################################################

#### 1) CREATE A LIST OF SPECIES BELONGUING TO ALL DATA SETS (TRAIT, TREE AND TAXA)
#### 2) CALCULATE PHYLOMETRICS
#### 3) ADD ENVIRONMENTAL DATA


library(geiger)
library(picante)
library(laser)
library(apTreeshape)
library(FD)
library(ape)
library(raster)
library(reshape)
library(sp)
library(maptools)
library(rgdal)
library(plyr)
library(phytools)
library(maps)
library(moments)
library(phylocom)
library(BioGeoBEARS)
library(phylobase)
library(foreach);library(parallel);library(doParallel)

rm(list=ls())

#setwd('/home/brunno/Dropbox/Doutorado Brunno/Manuscritos/R_scripts')
setwd('C:/Users/Brunno/Documents/brunno')
load("PD&FD.RData")
save.image("PD&FD.RData")

### 1 - DISTRIBUTIONAL DATA
dataraw <- read.csv("Working data/siteXsppXenv_1dgr.csv")
spdata <- names(dataraw)
spdata <- data.frame(spdata)
#dim(spdata) #5141 species

### 2 - PHYLOGENY 
tre <- read.tree("Working data/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
tredata <- data.frame(tre$tip.label)
#dim(tredata) #5364 tips

### 3 - TRAIT DATA
trait<-data.frame(read.table("Working data/imputedmammals_with_categorical_sub.csv",header=TRUE,sep=","))
traitdata <- trait[,1]
traitdata <- data.frame(traitdata)
#dim(traitdata) #4866 sps

# CREATING A LIST OF SPECIES WHICH BELONGS TO THE THREE DATA SETS
names(spdata)<-names(tredata)<-names(traitdata) <- "sps"

lista1 <- merge(spdata,tredata,by=intersect(names(spdata), names(tredata)))
lista2 <- merge(lista1,traitdata,by=intersect(names(lista1), names(traitdata))) # This is the list of species belonging to all three datasets
lista <-as.vector(lista2[,1])

rm(lista1,lista2)

### CREATING NEW OCCR, TREE AND TRAITS DATA WITH SPS BELONGING TO ALL DATASETS

# 1 - DISTRIBUTIONAL DATA
data <- dataraw[lista] #just species
data$Rich <- rowSums(data) #add rich column
bios <- dataraw[,1:21] #21th is the last column before sps
XY<-dataraw[,1:2]
metadata <- cbind(bios,data)
metadata <- subset(metadata, Rich > 4) # Exclude cells with poor SR. See "IS THERE ANY PROBLEM IN DELETE SPECIES POOR SITES?"
metadata <- subset(metadata, select=-c(Rich)) # Exclude Richness column
XY<-metadata[,1:2]
bios <- metadata[,3:21]
occr <-  metadata[,-1:-21]
occr <- occr[,colSums(occr) != 0] # Exclude species occurring nowhere
occr <- data.frame(t(occr))
names(occr) <- paste("Cell",1:ncol(occr),sep="") # This is the occurence table
lista <- as.vector(rownames(occr)) # New species list

rm(data)

# 2 - PHYLOGENY
phylo_nqro <- drop.tip(tre,lista) #Remove as espécies que quero e ficam as que não quero
tip_nqro<-as.vector(phylo_nqro$tip.label) #Lista das sps que nao quero
tre<-drop.tip(tre, tip_nqro) #Eliminar as sp que nao quero e ficam as que quero

rm(tip_nqro,phylo_nqro)

# 3 - TRAIT DATA
trait <- trait[trait$IUCN.binomial%in%lista,]
rownames(trait) <- as.vector(trait[,1])
trait <- subset(trait, select=-c(IUCN.binomial))#delete the column IUCN.binomial

# 4 - REARRANGE DATA
# The species data in the data frames containing the distribution and trait measurements 
# must be in the same row order as the names of the species in the tree object (mytree$tip.label). 
occr <- occr[match(tre$tip.label,rownames(occr)),]
trait <- trait[match(tre$tip.label,rownames(trait)),]
# check whether our community, tree and traits are in the same order.
all.equal(rownames(occr), tre$tip.label)
all.equal(rownames(occr), rownames(trait))

rm(traitdata,tredata,spdata)

#### SEEING SPECIES RICHNESS MAP
test <- cbind(metadata[,c('x','y')],Rich=rowSums(t(occr)))
coordinates(test)<-~x+y
gridded(test) <- TRUE
rasterTes <- raster(test)
plot(rasterTes,col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327","#e69500", "#ec2236"))(40))

#### CHECK FOR DIFFERENCES BETWEEN RAW DATA AND MATCHED DATA
test2 <- cbind(dataraw[,c('x','y')],Rich=rowSums(dataraw[,-1:-21]))
test2 <- test2[which(test2$x %in% test$x),]
test2 <- test2[which(test2$y %in% test$y),]
coordinates(test2)<-~x+y
gridded(test2) <- TRUE
rasterTes2 <- raster(test2)
plot(rasterTes2,col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327","#e69500", "#ec2236"))(40))

#### IS THIS SUBSET A BIASED GEOGRAPHY OF SPECIES RICHNESS?
#### WHERE DID WE LOOSE MORE DATA ON SPECIES?
raster3 <- rasterTes2/rasterTes
plot(raster3)
rast12<-na.omit(cbind(rasterTes2[],rasterTes[]))
cor.test(rast12[,1],rast12[,2])
#We have a subset of the original sp dist data but it's not a problem because our subset actualy 
#reflect the original data with a righ confidential level --> R2 = 0.99 ; p < 2.2e-16

#### IS THERE ANY PROBLEM IN DELETE SPECIES POOR SITES?
spnum <- 5
plot(rasterTes <= spnum)
hist(rowSums(dataraw[,-1:-21]), breaks=1000)
# Species poor sites (<5) are located in islands and high latitudinal zones.
# However, few points in deserts...
# Avoiding these site will not compromise our analysis.

rm(dataraw)

#### RUN PHYLOGENETIC DIVERSITY INDEXES

### Species Ages
# Get the age at each node
trc <- read.tree("Working data/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

# Get evolutionary distinctiveness
spp.ED <- evol.distinct(trc, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                               scale = FALSE, use.branch.lengths = TRUE)


### Phylometrics

Ncell <- ncol(occr)

SPD<-PHD<-GAM<-DIV<-MRD.t<-AGX<-AGM<-AGML<-TEH<-ED<-ED.a<-DR<-RBL<-
  SAG<-mSAG<-clades<-rep(NA,Ncell)

for (i in 1:Ncell){
  cat("\r",i,"of", Ncell)

  tipnames <- which(!t(occr[,i]))
  trx <- drop.tip(tre, tipnames) ### construct tree of all species within a cell

  #sub<-subset(spp.ages,Species%in%trx$tip.label)
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  
  SPD[i] <- length(trx$tip.label) ### cell richness
  #PHD[i] <- sum(trx$edge.length)### phylogenetic diversity within cell
  #GAM[i] <- gamStat(branching.times(trx),return.list=FALSE) ### gamma statistic (stemmy versus tippy tree)
  #IMY[i] <- colless(as.treeshape(trx),norm="yule")### tree imbalance (radiations and old lineages)
  #IMP[i] <- colless(as.treeshape(trx),norm="pda") ### tree imbalance (radiations and old lineages)
  #DIV[i] <- log(length(trx$tip.label))/(max(branching.times(trx)))  ### diversification rate (species accumulation over time)
  #MRD.t[i] <- mean(max(branching.times(trx))-sub$age) ### median root distance - how far from the base species arise
  #AGX[i] <- max(branching.times(trx)) ### maximum lineage age (new and old diversity)
  #AGM[i] <- mean(branching.times(trx)) ### median lineage age (new and old diversity)
  #AGML[i] <- mean(log(branching.times(trx))) ### median lineage age (new and old diversity)
  #TEH[i] <- sum(trx$edge.length) ### total evolutionary history -> total branching lengths of in the phylogeny of cell's species Davies et al. 2007 (PNAS)
  #ED[i] <- mean(subED[,2]) ### Evolutionary distinctiveness
  #ED.a[i] <- mean(evol.distinct(trx, type = c("fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)[,2]) ### Evolutionary distinctiveness assemblage
  DR[i] <- mean(1/subED[,2]) # Diversification rates in Jetz et al.2012
  
  #RBL[i] <- median((max(branching.times(trx))-branching.times(trx))/max(branching.times(trx)))
  
  #SAG[i] <- mean(sub$age)
  #mSAG[i] <- max(sub$age)
  
  #clades[i] <- length(getCladesofSize(trx, clade.size=2)) ### clade rich
  }

#MRD <- MRD_eliot(matrix2sample(t(occr)),tre)

#AGE <- residuals(loess(TEH~log(SPD))) # 'old' vs. 'young' diversity in Davies et al. 2008

#metadata <- cbind(data.frame(XY),SPD,PHD,GAM,DIV,MRD.t,AGX,AGM,AGML,TEH,ED,ED.a,DR,RBL,SAG,mSAG,clades,MRD=MRD$MRD,AGE)
metadata <- cbind(data.frame(XY),SPD,DR)

save.image("PD&FD.RData")

### COMMUNITY PHYLOGENETICS

COM1 <- ses.mpd (t(occr), cophenetic (tre), null.model = "taxa.labels", runs=99, iterations=100)
NRI  <- -1 * COM1[,"mpd.obs.z"]
MPD  <- COM1[,"mpd.obs"]
ses.MPD <- COM1[,"mpd.obs.z"]

#COM2 <- ses.mntd(t(occr), cophenetic (tre), null.model = "taxa.labels", runs=99, iterations=100)
#NTI  <- -1 * COM2[,"mntd.obs.z"]
#MNTD <- COM2[,"mntd.obs"]

#COM3 <-psv(t(occr), cophenetic (tre), compute.var=TRUE)
#PSV <- COM3[,"PSVs"]

#metadata <- cbind(metadata,NRI,MPD,NTI,MNTD,PSV)
metadata <- cbind(metadata,MPD,NRI,ses.MPD)

save.image("PD&FD.RData")


### FUNCTIONAL DIVERSITY INDEXES WITH SUBSET OF TRAITS (include categorical)
### including functional dispersion, evenesss, rao's entropy, etc. 

complet <- read.table('Working data/All_Mammal_trait_Data-9-6-13.txt',header=T)
complet <- complet[,c("BodyMass.g","HabitatBreadth","DietBreadth","LitSz","TrophicD","HabMode","ActivCycle")]
# % completeness for each trait
sapply(complet, function(x) 1-(sum(length(which(is.na(x))))/nrow(complet)))
#mean completeness
mean(sapply(complet, function(x) 1-(sum(length(which(is.na(x))))/nrow(complet))))
#sd completeness
sd(sapply(complet, function(x) 1-(sum(length(which(is.na(x))))/nrow(complet))))

# traits are:
# logBodyMass.g / TrophicD / HabMode / logHabitatBreadth / logDietBreadth
# ActivCycle / logLitSz

# Use the variation of Gower's Distance publicated by Pavoine in Oikos (2009)

# Make a file from the Quantitative Variable
tabQ <- trait[c("logBodyMass.g","logHabitatBreadth","logDietBreadth","logLitSz")]
# Now with the Nominal Variables
tabN <- trait[c("TrophicD","HabMode","ActivCycle")]

ktab1  <- ktab.list.df(list(data.frame(tabN),data.frame(tabQ)))
TraitDis <- dist.ktab(ktab1 , c("N", "Q"), c("scaledBYrange"))

#FD_ <-dbFD(trait, t(occr)) # dont work...

FD <- fdisp(TraitDis, t(occr))

str(FD)

metadata <- cbind(metadata,FDis=FD$FDis)

## Calculate ses.FD
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

obs.null.output <- foreach(i=1:100, .combine='cbind', .packages=c('FD','picante')) %dopar% { 
  fdisp(TraitDis, randomizeMatrix(t(occr), null.model = "richness"))$FDis
}
# stop the cluster
stopCluster(cl)

obs.null.output <- data.frame(obs.null.output)

ses.FD <- NA
for(i in 1:dim(occr)[2]){ cat("\r", i, 'from', dim(occr)[2])
  ses.FD[i] <- (FD$FDis[i] - mean(as.numeric(obs.null.output[i,])))/sd(as.numeric(obs.null.output[i,]))
}

metadata <- cbind(metadata, ses.FD)

######################################################
### Calculate PCA traits

# transform categorical in numeric
# in TrophicD (her omni car - 1 2 3 .:. represent trophic level)
# in HabMode (fos aqu ter arb vol  - 1 2 3 4 5 .:. represent vertical level)
# in ActivCycle (D B N - 1 2 3.:. represent diel activity)

tabN2<- tabN

tabN2$TrophicD <- gsub('her',1,tabN2$TrophicD)
tabN2$TrophicD <- gsub('omni',2,tabN2$TrophicD)
tabN2$TrophicD <- gsub('car',3,tabN2$TrophicD)
unique(tabN2$TrophicD)
tabN2$TrophicD <- as.numeric(tabN2$TrophicD)

tabN2$HabMode <- gsub('fos',1,tabN2$HabMode)
tabN2$HabMode <- gsub('aqu',2,tabN2$HabMode)
tabN2$HabMode <- gsub('ter',3,tabN2$HabMode)
tabN2$HabMode <- gsub('arb',4,tabN2$HabMode)
tabN2$HabMode <- gsub('vol',5,tabN2$HabMode)
unique(tabN2$HabMode)
tabN2$HabMode <- as.numeric(tabN2$HabMode)

tabN2$ActivCycle <- gsub('D',1,tabN2$ActivCycle)
tabN2$ActivCycle <- gsub('B',2,tabN2$ActivCycle)
tabN2$ActivCycle <- gsub('N',3,tabN2$ActivCycle)
unique(tabN2$ActivCycle)
tabN2$ActivCycle <- as.numeric(tabN2$ActivCycle)

trait3 <- cbind(tabQ,tabN2) # bind continuous and categorical traits

pcatraits <- princomp(trait3, scale.unit=TRUE)

pcatraits$loadings

summary(pcatraits)
# % of variance
#PC1 = 0.71
#PC2 = 0.12
# cumulative % of the first two pca = 0.8346726
pc1Traits <- data.frame(pc1Traits=pcatraits$scores[,1])

phevoBM <- NA
phevoBMaic <- NA
phevoOU <- NA
phevoOUaic <- NA
phevoEB <- NA
phevoEBaic <- NA

rownames(pc1Traits) <- lista

for (i in 1:Ncell){
  cat("\r",i,"of", Ncell)
  
  tipnames <- which(!t(occr[,i]))
  trx <- drop.tip(tre, tipnames) ### construct tree of all species within a cell
  
  fitraits <- pc1Traits[trx$tip.label,]
  names(fitraits) <- trx$tip.label
  
  phevoBM[i] <- fitContinuous(trx,dat=fitraits,model = c("BM"))$opt$sigsq
  phevoBMaic[i] <- fitContinuous(trx,dat=fitraits,model = c("BM"))$opt$aic
  phevoOU[i] <- fitContinuous(trx,dat=fitraits,model = c("OU"))$opt$sigsq
  phevoOUaic[i] <- fitContinuous(trx,dat=fitraits,model = c("OU"))$opt$aic
  phevoEB[i] <- fitContinuous(trx,dat=fitraits,model = c("EB"))$opt$sigsq
  phevoEBaic[i] <- fitContinuous(trx,dat=fitraits,model = c("EB"))$opt$aic
}

metadata <- cbind(metadata, phevoBM, phevoBMaic, phevoEB, phevoEBaic, phevoOU, phevoOUaic)

  
### ADD ENVIRONMENTAL DATA

coordinates(XY)=~x+y
mollP<-"+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
projection(XY)<-mollP

# REFERENCE PROJECTION MAP
mapa <- raster("F:/GIS/Mollweide/mundo_1Dg.grd") 
#AET
aet <- raster("F:/GIS/#Environment/AET/ataet00.tif")
crs(aet) <- "+proj=longlat +datum=WGS84" # Set projection to WGS84.
extent(aet) <- extent(c(-180, 180, -90, 90)) # Set extention as WGS84. The original extention doesnt make sense
aet <- projectRaster(aet, mapa)
aet <- extract(aet,metadata[,1:2],method='simple')
aet[which(aet==0)] = NA # NA in original map is zero. Replace zero by NA.

#Continents
mundi <- readShapePoly("F:/GIS/Shp files/Mundi.shp")
crs(mundi) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
mundi <- spTransform(mundi, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
extent(mapa) <- extent(mundi)
mundinames<-sort(mundi$CONTINENT_) # put alphabetic order because rasterize will put number in acordance with alphabetic order
mundiraster <- rasterize(mundi,mapa,'CONTINENT_')
mundi <- over(SpatialPoints(XY),mundi)
mundi <- as.vector(extract(mundiraster,metadata[1:2],method='simple'))
mundi <- as.character(mundinames)[ match(mundi, c(1:9))]

### Realms
#1) Get names
Reanames <- readShapePoly("F:/GIS/CMEC regions & realms/newRealms.shp")
crs(Reanames) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
extent(mapa) <- extent(Reanames)
Reanames<-data.frame(Reanames$Realm)
names(Reanames)<-'Realm'
Reanames<-cbind(Realms=rownames(Reanames),Reanames)
#2) Extract values
Realm <- readShapePoly("F:/GIS/CMEC regions & realms/Realms.shp")
crs(Realm) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
realmnames<-sort(Realm$fullupgmar) # put alphabetic order because rasterize will put number in acordance with alphabetic order
extent(mapa) <- extent(Realm)
Realm <- rasterize(Realm,mapa,'fullupgmar')
Realm <- projectRaster(Realm, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
Realm <- as.vector(extract(Realm,XY,method='simple'))

#Some cells comprise more than one Realm and the extract function give a medium value for that cell. 
#This value have no meaning as we need the value of the Realm for each cell.
check.integer <- function(N){ # function for check if number is integer
  !length(grep("[^[:digit:]]", as.character(N)))
}

for (i in seq_along(Realm)){ # put NA in non-integer numbers
  if (check.integer(Realm[i]))
  {Realm[i]<-Realm[i]
  }else{Realm[i]<-NA}
}

#Dont know why...but some codes come with whitespaces -To see use: unique(Realm)
trim <- function( x ) { #Function to remove whitespaces
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}
Realm<-trim(Realm)

#Give names for Realm code...
Realm <- as.character(Reanames$Realm)[ match(Realm, c(1:11))]

### Solve problem with NAs between realms
Realm2 <- readShapePoly("F:/GIS/CMEC regions & realms/newRealms.shp")
crs(Realm2) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
Realm2 <- spTransform(Realm2, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
coordinates(XY)=~X+Y
crs(XY) <-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
Realm2<-over(XY,Realm2)
Realm2<-as.vector(Realm2$Realm)

#Give names for NAs based on the previous 
Realm3<-NA
for(i in 1:length(Realm2)){
  if(is.na(Realm[i])){
  Realm3[i]<-Realm2[i]
  }else{
    Realm3[i]<-Realm[i]
  } 
}

Realm<-Realm3
rm(Realm2,Realm3)


metadata <- cbind(metadata, aet, Realm)

save.image("PD&FD.RData")

####

rm(list=ls())
# install.packages("rgdal")
# install.packages("sp")
library(sp)
library(car)
library(rgdal)
library(spdep)                               # Spatial Analysis; also opens libraries "maptools" and "sp"
library(RColorBrewer)                        # see "http://colorbrewer2.org/"
library(classInt)                            # functions of data classification
library(maptools)

# list <- c("sp", "care", "rgdal", "spdep", "RColorBrewer", "classInt", "maptools")
# install.packages(list)
##### Import the given functions to plot
plotColorQual <- function(var.name,shape,my.title="",
                          my.legend=deparse(substitute(var.name)),
                          addToMap=F) {
  ##
  ## Plot a qualitative colors for factor "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  if (!is.factor(var.name)) stop("plotColorQual: Not a factor.")
  
  qualVal <- as.numeric(unclass(var.name))
  qualName <- levels(var.name)
  pal.Qual <- brewer.pal(12,"Set3")
  map.col <- pal.Qual[qualVal]
  
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend, legend=qualName,
         fill=pal.Qual[1:length(qualName)],bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotColorQual

plotColorRamp <- function(var.name,shape,n.breaks=8,my.title="",
                          my.legend=deparse(substitute(var.name)),
                          addToMap=F) {
  ##
  ## Plot a color ramp variable "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  
  ## define breaks and color assignment
  q.breaks <- classIntervals(var.name, n=n.breaks, style="quantile")
  pal.YlOrRd <- brewer.pal(n.breaks, "Oranges")
  #pal.YlOrRd <- brewer.pal(n.breaks, "YlOrRd")
  map.col <- pal.YlOrRd[findInterval(var.name,q.breaks$brks,rightmost.closed=T)]
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend,legend=leglabs(round(q.breaks$brks,digits=3)),
         fill=pal.YlOrRd,bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotColorRamp 

plotBiPolar <- function(var.name,shape,
                        neg.breaks=4,pos.breaks=neg.breaks,break.value=0,
                        my.title="",my.legend=deparse(substitute(var.name)),
                        addToMap=F) {
  ##
  ## Plot bipolar map theme for variable "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  
  ## define quantile breaks and color assignment
  q.neg.breaks <- classIntervals((var.name[var.name < break.value]), n=neg.breaks, style="quantile")
  q.pos.breaks <- classIntervals((var.name[var.name > break.value]), n=pos.breaks, style="quantile")
  q.breaks <- c(q.neg.breaks$brks[-(neg.breaks+1)],break.value,q.pos.breaks$brks[-1])     # combine neg and pos over zero
  
  pal.neg <- brewer.pal(neg.breaks, "Blues")
  pal.pos <- brewer.pal(pos.breaks, "Reds")
  pal <- c(rev(pal.neg),pal.pos)                                                # combine palettes
  
  map.col <- pal[findInterval(var.name,q.breaks,rightmost.closed=T)]
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend,legend=leglabs(round(q.breaks,digits=3)),
         fill=pal,bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotBiPolar

####

directory <- "C:\\Users\\hxw170008\\Downloads\\TXCnty2018"
setwd(directory)

### Read and transformation
TX <- rgdal::readOGR(dsn=getwd(), layer="TXCnty", integer64="warn.loss")
NEIGHBORS <- rgdal::readOGR(dsn=getwd(), layer="TXNeighbors", integer64="warn.loss")
HWY <- rgdal::readOGR(dsn=getwd(), layer="InterStateHwy", integer64="warn.loss")


TX <- spTransform(TX, CRS("+proj=longlat"))
NEIGHBORS <- spTransform(NEIGHBORS, CRS("+proj=longlat"))
HWY <- spTransform(HWY, CRS("+proj=longlat"))


### Calculate the votion population
### Vote/total people voted

TX$VotingPOP <-  TX$TURNOUT16 * TX$REGVOT16
TX$clinton_rate <- TX$CLINTONVOT/TX$VotingPOP
TX$trump_rate <- TX$TRUMPVOT16/TX$VotingPOP


# Drop the smallist county
plot(TX$VotingPOP, ylim = c(0, 200))
#drop the row with min Voting pop
data <- as.data.frame(TX)
toDROP <- which(data[,which(colnames(data)=="VotingPOP")] == min(data$VotingPOP))
TX <- TX[-c(toDROP),]
## Identified the small counties and removed them

#Try this part later

# ## Why is Bolzano-Bozen an extreme observation? Shall we delete it?
# ##
# idx.max <- which.max(abs(fertResid))        # Get index of a record with "outlying" observation 
# ## Map potential outlier
# extremeObs <- rep(0, length(fertResid))
# extremeObs[idx.max] <- 1
# extremeObs <- factor(extremeObs, labels=c("inPop","extreme"))
# table(extremeObs)
# plot(neig.shp,axes=T,col=grey(0.9),border="white",                 # background: neighboring countries
#      xlim=prov.bbox[1,],ylim=prov.bbox[2,])                     
# plotColorQual(extremeObs, prov.shp, my.title="Potential Outlier",
#               my.legend="Outliers", addToMap=T)        
# box()  
# 
# ## Inspect and delete outlier
# prov.shp@data[idx.max,]                     # List info of record centre is outlier
# prov.shp <- prov.shp[-idx.max, ]            # Delete extreme observation from shapefile 
# 


# distribution (how to improve it)
hist(TX$trump_rate)
hist(TX$clinton_rate)
# Skewness test show trump is better
e1071::skewness(TX$trump_rate)
e1071::skewness(TX$clinton_rate)

#Box cox transformation of rate
lambda1 <- powerTransform(trump_rate~1, data=TX)$lambda
lambda2 <- powerTransform(clinton_rate~1, data=TX)$lambda

# Test if use log would be better!?
TX$bc.trump <- car::bcPower(TX$trump_rate, lambda=lambda1)
TX$bc.clinton <- car::bcPower(TX$clinton_rate, lambda=lambda2)


####Things to add use anova to test if log is sufficient


#Use trump because of lower skewness 
hist(data$bc.trump)
hist(data$bc.clinton)

TX$ID<-as.factor(TX$ID)

### Try to map
prov.bbox <- bbox(TX)                                        # province bounding box for map region
plot(NEIGHBORS,axes=T,col=grey(0.9),border="white",                 # background: neighboring countries
     xlim=prov.bbox[1,],ylim=prov.bbox[2,]) 
# plotColorQual(TX$REGION,TX,my.title="Texas Counties", my.legend = "TX COUNTIES"
#               , addToMap=T) 

plotColorRamp(TX$trump_rate,TX,n.breaks=8,    # second add map
              my.title="Spatial Pattern of Trump Vote Rate",
              my.legend="Trump Rate",addToMap=T)     # addToMap=T over-plots provinces over neighbors

#Plot the link
TX.centroid <- coordinates(TX) 
TX.link <- poly2nb(TX, queen=F)

plot(NEIGHBORS,axes=T,col=grey(0.9),border="white",                 # background: neighboring countries
     xlim=prov.bbox[1,],ylim=prov.bbox[2,]) 
plot(TX,col="palegreen3" ,border=grey(0.9), axes=T, add=T) # Second plot areas
plot(TX.link,coords=TX.centroid, pch=19, cex=0.1,            # Third plot links focused at centroids
     col="blue", add=T)
title("Augmented Spatial Links among counties")                 # Forth add title
box()  
## Spatial Residual Analysis


#Task 1 
# Drop the smallist county
plot(data$VotingPOP, ylim = c(0, 200))
#drop the row with min Voting pop
data <- data[-c(which(data[,which(colnames(data)=="VotingPOP")] == min(data$VotingPOP))),]

#Task 2
# get potential regressors


# Ploting is not really useful
# setwd("D:\\GIS Statitics\\Lab works\\Labwork 4\\plots")
# list <- colnames(data)
# attach(data)
# for(i in 1:length(list)){
#   jpeg(paste(list[i],".png"))
#   plot(data$bc.trump, data[,i], xlab = "trump", ylab = list[i])
#   dev.off() 
# }
# setwd(directory)

# Run f test for all parameters vs. Trump
# Functions Get f-statistics
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

list<-colnames(data)
listofnotable <- c()
for(i in 2:length(list)){
  mod<-lm(formula(paste("bc.trump","~",list[i])),data = data)
  x <- lmp(mod)
  if(is.na(x)){
    x <- 1
  }
  if(x < 0.05){
    listofnotable <- c(listofnotable,list[i])
  }
}

reduced_list<- 
  
  ## to do (use weighted.residuals for heteroskedesticity)
  # resfun	
  # default: weighted.residuals; the function to be used to extract residuals from the lm object, may be residuals, weighted.residuals, rstandard, or rstudent

##The principal coordinates of neighbor matrices (PCNM) analysis was performed as this method can decompose the total spatial variation into a finite set of explanatory spatial variables, each of which corresponds to a specific spatial structure or scale
##Before perform the PCNM analysis, please try to library these packages as following:
library(ape)
library(spdep)
library(ade4)
library(spacemakeR)
library(AEM)
library(PCNM)
library(SoDA)

##change to the directory on your computer that contains Latitude and longitude of each sample
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/PCNM")

##Read in the longitude and Latitude of each sample
coordinate <- read.csv("coordinate.csv", row.names=1)

##Convert the coordinates in latitude and longitude form to corresponding coordinates in X (east-west) and Y (north-south) distances along the surface of the earth
coor_XY <- geoXY(coordinate[,1], coordinate[,2], unit = 1000)#X=Latitude, Y=longitude, unit=1000 causes the coordinates to be in kilometers

##Constructing the Euclidean distance matrix among samples
distance_euc <- dist(coor_XY)

##perform PCNM analysis
dis.PCNM.auto <- PCNM(distance_euc)
summary(dis.PCNM.auto)

##Expected value of Moran's I, represent no spatial correlation
dis.PCNM.auto$expected_Moran
dis.PCNM.auto$Moran_I

##Eigenfunctions with positive spatial correlation
dis.select <- which(dis.PCNM.auto$Moran_I$Positive==TRUE)
length(dis.select)#Number of PCNM with I > E(I)
dis.PCNM.pos <- as.data.frame(dis.PCNM.auto$vectors)[,dis.select]

##write the spatial variables generated from PCNM analysis
write.csv(dis.PCNM.pos,"spatial_variables.csv")

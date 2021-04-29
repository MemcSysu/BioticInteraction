##change to the directory on your computer that contains the data of normalized bacterial and pico-nano eukaryotic OTU table
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/Mantel")

##Load Package
library(vegan)
library(basicTrendline)

##Read in the normalized OTU table in csv format
Spe_16S <- read.csv("OTU_table_16S_norm.csv", row.names = 1)
bacteria <- t(Spe_16S)#sample in row, OTU in column
Spe_18S <- read.csv("OTU_table_18S_norm.csv", row.names = 1)
eukaryotes <- t(Spe_18S)#sample in row, OTU in column

##Compute dissimilarity indices
bacteria_bray <- vegdist(bacteria, method='bray')#'bray' represent bray-curtis dissimilarity index
eukaryotes_bray <- vegdist(eukaryotes, method='bray')

#mantel test
mantel(bacteria_bray, eukaryotes_bray)
pdf("Mantel_test_plot.pdf")
plot(eukaryotes_bray,bacteria_bray)
trendline(eukaryotes_bray, bacteria_bray, model="line2P", linecolor = "black", ePos.x = "topleft", summary=TRUE, eDigit=5)
dev.off()

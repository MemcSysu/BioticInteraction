##change to the directory on your computer that contains the data of raw bacterial OTU table and biotic traits of pico-nano eukaryotic community
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/network")

##Load Package
library(psych)

##Filter out the bacterial OTUs with relative abundance less than 0.01%
bac_raw_otu <- read.csv("OTU_table_16S.csv", row.names = 1)
bac_raw_otu_relative <- apply(bac_raw_otu, 2, function(x){x/sum(x)*100})
bac_otu_relative_remove0.01 <- bac_raw_otu_relative[apply(bac_raw_otu_relative, 1, mean) > 0.01,]
Bac_remove0.01 <- t(bac_otu_relative_remove0.01)#samlpes in row, OTU in column

##Find the correlations between bacteria and pico-nano eukaryotes
Bio <- read.csv("biotic_traits.csv", row.names = 1)#read in biotic traits
occor_none <- corr.test(Bac_remove0.01, Bio, use="pairwise", method="spearman", adjust="none", alpha=0.05)
occor.r <- occor_none$r
occor.p <- occor_none$p
occor.r[occor.p>0.01] <- 0
write.csv(occor.r[which(rowSums(occor.r) != 0),],file="correlation_result.csv")#export the r values of correlation
occor.p[occor.p>0.01]<-0
write.csv(occor.p[which(rowSums(occor.p) != 0),],file="pvalueresult.csv")#export the p values of correlation

##Note that further analysis and visualization are finished in Cytoscape software 

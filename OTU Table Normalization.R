##Normalize the OTU table using the 'edgeR' and 'limma' packages
##To install core packages, type the following in an R command window:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##Install 'edgeR' packages, with
BiocManager::install("edgeR")
install.packages('limma')

library(edgeR)
library(limma)

##change to the directory on your computer that contains the OTU tables of bacteria and pico-nanoeukaryotes
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/normalization")

##read in the OTU table of bacteria
OTU_table_16S <- read.csv("OTU_table_16S.csv", row.names = 1)#make your OTU table in .csv or table format with the OTU name in row and sample name in column

##Creates a DGEList object from a table of counts (rows=OTUs, columns=samples)
OTU_16S_list <- DGEList(counts=OTU_table_16S)

##Calculate normalization factors to scale the raw library sizes
OTU_16S_Normfactors <- calcNormFactors(OTU_16S_list)
OTU_16S_Normfactors$samples#look the normalization factors of each sample

##Normalize each library size using the normalization factors
normalization_16S <- OTU_16S_Normfactors$samples$lib.size*OTU_16S_Normfactors$samples$norm.factors
for(i in 1:20){
  OTU_table_16S[,i] <- OTU_table_16S[,i]/normalization_16S[i]
}

##write the normalized bacterial OTU table in csv format
write.csv(OTU_table_16S, "OTU_table_16S_norm.csv")


##read in the OTU table of pico-nano eukaryotes
OTU_table_18S <- read.csv("OTU_table_18S.csv", row.names = 1)#make your OTU table in .csv or table format with the OTU name in row and sample name in column

##Creates a DGEList object from a table of counts (rows=OTUs, columns=samples)
OTU_18S_list <- DGEList(counts=OTU_table_18S)

##Calculate normalization factors to scale the raw library sizes
OTU_18S_Normfactors <- calcNormFactors(OTU_18S_list)
OTU_18S_Normfactors$samples#look the normalization factors of each sample

##Normalize each library size using the normalization factors
normalization_18S <- OTU_18S_Normfactors$samples$lib.size*OTU_18S_Normfactors$samples$norm.factors
for(i in 1:20){
  OTU_table_18S[,i] <- OTU_table_18S[,i]/normalization_18S[i]
}

##write the normalized pico-nano eukaryotic OTU table in csv format
write.csv(OTU_table_18S, "OTU_table_18S_norm.csv")


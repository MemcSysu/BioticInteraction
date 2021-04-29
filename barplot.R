##change to the directory on your computer that contains the data of relative abundance of major taxa
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/barplot")

##load this package
##if not already installed, use install.packages('RColorBrewer')
library(RColorBrewer)

##barplot for major taxa of bacterial community
##Read in relative abundance data of mjor bacterial taxa
barplot_16S <- read.csv("barplot_16S.csv", row.names = 1)
barplot_16S_matrix <- as.matrix(barplot_16S)

##plot and save the barplot for bacterial 16S as pdf format
pdf("barplot_16S.pdf")
opar <- par(no.readonly = TRUE)
par(mai = c(0.5,0.5,0.1,2.0))#Resize the graphic margin to make enough space for the legend
xy_16S <- par("usr")
barplot(barplot_16S_matrix, col = c(brewer.pal(8,"Dark2")[1:6], brewer.pal(12,"Set3")[1:12],
                                    brewer.pal(8,"Dark2")[8]), legend=rownames(barplot_16S_matrix),
                                    axis.lty = 1, args.legend = list(x=xy_16S[2]+xinch(5.0),
                                                                     y=xy_16S[4]+yinch(5.0), ncol=1, xpd=T,
                                                                     x.intersp=0.4, y.intersp=1.2))
par(opar)
dev.off()

##barplot for major taxa of pico-nano eukaryotic community
##Read in relative abundance data of major eukaryotic taxa
barplot_18S <- read.csv("barplot_18S.csv", row.names = 1)
barplot_18S_matrix <- as.matrix(barplot_18S)

##plot and save the barplot for eukaryotic 18S as pdf format
pdf("barplot_18S.pdf")
par(mai=c(0.5,0.2,1.6,1.6))
xy_18S <- par("usr")
barplot(barplot_18S_matrix, col = c(brewer.pal(9,"Set1")[1:8], brewer.pal(8,"Set2")[1:2],
                                   brewer.pal(9,"Set1")[9]), legend=rownames(barplot_18S_matrix),
                                   axis.lty = 1, args.legend = list(x=xy_18S[2]+xinch(4.9),
                                                                    y=xy_18S[3]+yinch(3.5), ncol=1, xpd=T,
                                                                    x.intersp=0.4, y.intersp=1))
par(opar)
dev.off()

##barplot for functional composition of pico-nano eukaryotic community
##Read in relative abundance data of functional composition of pico-nano eukaryotic community
barplot_functional <- read.csv("barplot_functional.csv", row.names = 1)
barplot_functional_matrix <- as.matrix(barplot_functional)

##plot and save the barplot for functional composition as pdf format
pdf("barplot_functional.pdf")
par(mai=c(0.5,0.2,1.6,1.6))
xy_functional <- par("usr")
barplot(barplot_functional_matrix, col = c(brewer.pal(9,"Set1")[1:4],brewer.pal(9,"Set1")[9]),
                                          legend=rownames(barplot_functional_matrix),
                                          axis.lty = 1, args.legend = list(x=xy_functional[2]+xinch(4.9),
                                                                           y=xy_functional[3]+yinch(3.5), xpd=T,
                                                                           x.intersp=1, y.intersp=1))
par(opar)
dev.off()
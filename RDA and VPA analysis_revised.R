##change to the directory on your computer that contains the data of normalized bacterial OTU table, spatial variables generated from PCNM analysis, environmental factors and biotic traits of pico-nano eukaryotic community. The biotic traits are shown in Table S2 and Table S3 of the paper
##note that the 'slash' needs to be changed to a forward slash like this /
setwd("E:/Interaction/quantify")

##Load Package
library(vegan)

##Detrended correspondence analysis (DCA) is performed first to determine linear or unimodal model
##Read in the normalized bacterial OTU table in csv format
spe16 <- read.csv("OTU_table_16S_norm.csv", row.names = 1)
bacteria <- t(spe16)#sample in row, OTU in column
decorana(bacteria)#According to the length of the first ordination axis,RDA was selected

##identify the environmental factors that significantly co-vary with the bacterial community structure
Env <- read.csv("Env.csv", row.names = 1)#read in environmental factors data
scaleenv <- as.data.frame(scale(Env))#environmental variables are z score©\transformed
rda_env <- rda(bacteria ~ ., scaleenv)#RDA of the normalized 16S OTU table constrained by all the environmental variables contained in scaleenv

##To search for parsimony and reduce the explanatory variables that possible strong linear correlations in the RDA model, forward selection, a method often applied in RDA, is carried out using the function 'ordistep', see reference 52 in paper
step.forward_env <- ordistep(rda(bacteria ~ 1, data=scaleenv),scope = formula(rda_env), direction="forward", permutations = 1000)

##Construct parsimonious model
step.forward_env$anova
attr(step.forward_env$terms,"term.labels")
Env_signif <- as.data.frame(scaleenv[,c(attr(step.forward_env$terms,"term.labels"))])
new.rda_env <- rda(bacteria ~ ., data=Env_signif)
RsquareAdj(new.rda_env)$r.squared#r.squared of new RDA model
RsquareAdj(new.rda_env)$adj.r.squared#adjusted r.squared of new RDA model

#Calculate the value of AICc for Environmental model
Env_K <- length(Env_signif)
AIC_Env <- extractAIC(new.rda_env)[2]#extract the AIC value of new RDA model
AICc_Env <- sum(AIC_Env,2*Env_K*(Env_K+1)/(19-Env_K))#AICc = AIC+[2k(k+1)/(n-k-1)],K is the number of selected variables, n is the number of samples (also called observations)

permutest(new.rda_env, permu=999)#calculate the significance of new RDA model
vif.cca(new.rda_env)#calculate the variance inflation factors (VIF) of all environmental variables in new RDA model


##identify the spatial variables that significantly co-vary with the bacterial community structure
Spa <- read.csv("spatial_variables.csv", row.names = 1)#read in spatial vaiables data from PCNM analysis
rda_spa <- rda(bacteria ~ ., Spa)#RDA of the normalized 16S OTU table constrained by all the spatial variables contained in Spa

##search for parsimonious RDA model
step.forward_spa <- ordistep(rda(bacteria ~ 1, data=Spa),scope = formula(rda_spa), direction="forward", permutations = 1000)

##Construct parsimonious model
step.forward_spa$anova
attr(step.forward_spa$terms,"term.labels")
Spa_signif <- as.data.frame(dbmem_16S_frame[,c(attr(step.forward_spa$terms,"term.labels"))])
new.rda_spa <- rda(bacteria ~ ., data=Spa_signif)
RsquareAdj(new.rda_spa)$r.squared
RsquareAdj(new.rda_spa)$adj.r.squared

#Calculate the value of AICc for Spatial model
Spa_K <- length(Spa_signif)
AIC_Spa <- extractAIC(new.rda_spa)[2]#extract the AIC value of new RDA model
AICc_Spa <- sum(AIC_Spa,2*Spa_K*(Spa_K+1)/(19-Spa_K))#AICc = AIC+[2k(k+1)/(n-k-1)],K is the number of selected variables, n is the number of samples (also called observations)

permutest(new.rda_spa, permu=999)#calculate the significance of new RDA model
vif.cca(new.rda_spa)#calculate the variance inflation factors (VIF) of all spatial variables in new RDA model


##identify the biotic traits that significantly co-vary with the bacterial community structure
Bio <- read.csv("biotic_traits.csv", row.names = 1)#read in biotic traits
rda_bio <- rda(bacteria ~ ., Bio)#RDA of the normalized 16S OTU table constrained by all the biotic traits contained in Bio

##search for parsimonious RDA model
step.forward_bio <- ordistep(rda(bacteria ~ 1, data=Bio),scope = formula(rda_bio), direction="forward", permutations = 1000)

##Construct parsimonious model
step.forward_bio$anova
attr(step.forward_bio$terms,"term.labels")
Euk_signif <- Bio[,c(attr(step.forward_bio$terms,"term.labels"))]
new.rda_bio <- rda(bacteria ~ ., data=Euk_signif)
RsquareAdj(new.rda_bio)$r.squared
RsquareAdj(new.rda_bio)$adj.r.squared

#Calculate the value of AICc for Biotic model
Euk_K <- length(Euk_signif)
AIC_Euk <- extractAIC(new.rda_bio)[2]#extract the AIC value of new RDA model
AICc_Euk <- sum(AIC_Euk,2*Euk_K*(Euk_K+1)/(19-Euk_K))#AICc = AIC+[2k(k+1)/(n-k-1)],K is the number of selected variables, n is the number of samples (also called observations)

permutest(new.rda_bio, permu=999)
vif.cca(new.rda_bio)#calculate the variance inflation factors (VIF) of all Biotic variables in new RDA model, Note that if a variable has a VIF value greater than 10, please identify the biotic traits that significantly co-vary with the bacterial community structure using the functions 'rda', 'add1' and 'update' according to the P value and AIC value for manual selection and updating

##calculate the adjusted r.squared of each biotic trait
RsquareAdj(rda(bacteria ~ Dino_Group_II, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota + Basidiomycota, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae + Photo_Parasite, data=Euk_signif))
RsquareAdj(rda(bacteria ~ Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae + Photo_Parasite + Telonemia, data=Euk_signif))


##To quantify the relative contributions of environmental variables and spatial factors to the variation in bacterial community structure, We firstly constructed the RDA model with the significantly environmental factors and spatial variables selected using forward selection, and then used the 'varpart' function to quantify the contributions of these two parts
Env_Spa <- data.frame(Env_signif,Spa_signif)#Create a data.frame containing Temperature and PCNM_V1
names(Env_Spa) <- c("Temperature","PCNM_V1")#Modify the column name
row.names(Env_Spa) <- row.names(Spa)#Modify the row name
new.rda_Env_Spa <- rda(bacteria ~ ., data=Env_Spa)#construct new rda model with Temperature and PCNM_V1, note that proportion explained for each RDA axis can be obtained from the function 'summary' 
RsquareAdj(new.rda_Env_Spa)$r.squared
RsquareAdj(new.rda_Env_Spa)$adj.r.squared

Env_Spa_K <- length(Env_Spa)
AIC_Env_Spa <- extractAIC(new.rda_Env_Spa)[2]
AICc_Env_Spa <- sum(AIC_Env_Spa,2*Env_Spa_K*(Env_Spa_K+1)/(19-Env_Spa_K))

permutest(new.rda_Env_Spa, permu=999)
vif.cca(new.rda_Env_Spa)

##calculate the adjusted r.squared of Temperature and PCNM_V1
RsquareAdj(rda(bacteria~Temperature, data=Env_Spa))
RsquareAdj(rda(bacteria~Temperature + PCNM_V1,data=Env_Spa))


##Draw RDA plot including Temperature and PCNM_V1
site_scores_Env_Spa <- scores(new.rda_Env_Spa, display="lc")#Obtain the scores of samples(sites), "lc" for fitted site scores (linear combinations of explanatory variables)
group_Env_Spa <- c(1,3,6,6,1,4,3,7,5,3,2,2,7,4,6,8,5,2,7,7)#Group the samples
grouptotal <- cbind(site_scores_Env_Spa, group_Env_Spa)
factor(group_Env_Spa)#encode the values of group as factors
pdf("RDA_Env_Spa.pdf")
plot(new.rda_Env_Spa, display=c("lc","cn"),type="n")#choose the elements to be plotted, using the argument display=c(), "lc" for fitted site scores (linear combinations of explanatory variables), "cn" for constraints (i.e. the explanatory variables)
for(i in 1:4){
  points(grouptotal[group_Env_Spa==i,],pch=(14+i),cex=1.5)
}
for(i in 5:7){
  points(grouptotal[group_Env_Spa==i,],pch=(i-5),cex=1.5)
}
points(grouptotal[group_Env_Spa==8,1],grouptotal[group_Env_Spa==8,2],pch=5,cex=1.5)
text(new.rda_Env_Spa, display=c("cn"),col="black")
dev.off()


##Variation partitioning with environmental factors (Temperature) and Spatial variables (PCNM_V1)
varpart_Env_Spa <- varpart(bacteria, Env_signif, Spa_signif)
pdf("VPA_Env_Spa.pdf")
plot(varpart_Env_Spa,digits=3)
dev.off()

##Test of fractions of environmental factors and spatial variables
anova.cca(rda(bacteria,Env_Spa[,1]), step = 1000)
anova.cca(rda(bacteria,Env_Spa[,2]), step = 1000)


##To quantify the relative contributions of environmental variables, spatial factors and biotic traits to the variation in bacterial community structure, We constructed the RDA model with the significantly environmental factors, spatial variables and biotic traits selected using forward selection, and then used the 'varpart' function to quantify the contributions of these three parts
Env_Spa_Bio <- data.frame(Env_signif, Spa_signif, Euk_signif)#Create a data.fram containing Temperature, PCNM_V1 and seven biotic traits
new.rda_Env_Spa_Bio <- rda(bacteria ~ ., data=Env_Spa_Bio)
RsquareAdj(new.rda_Env_Spa_Bio)$r.squared
RsquareAdj(new.rda_Env_Spa_Bio)$adj.r.squared

#Calculate the value of AICc for Env_Spa_Bio model
Env_Spa_Bio_K <- length(Env_Spa_Bio)
AIC_Env_Spa_Bio <- extractAIC(new.rda_Env_Spa_Bio)[2]
AICc_Env_Spa_Bio <- sum(AIC_Env_Spa_Bio,2*Env_Spa_Bio_K*(Env_Spa_Bio_K+1)/(19-Env_Spa_Bio_K))

permutest(new.rda_Env_Spa_Bio, permu=999)
vif.cca(new.rda_Env_Spa_Bio)

##calculate the adjusted r.squared of Temperature, PCNM_V1 and seven biotic traits
RsquareAdj(rda(bacteria ~ Temp, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota + Basidiomycota, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae + Photo_Parasite, data=Env_Spa_Bio))
RsquareAdj(rda(bacteria ~ Temp + PCNM_V1 + Dino_Group_II + Ascomycota + Basidiomycota + Photo_Mixo + Phaeophyceae + Photo_Parasite + Telonemia, data=Env_Spa_Bio))


##Draw RDA plot including Temperature, PCNM_V1 and seven biotic traits
site_scores_Env_Spa_Bio <- scores(new.rda_Env_Spa_Bio, display="lc")
group_Env_Spa_Bio <- c(1,3,6,6,1,4,3,7,5,3,2,2,7,4,6,8,5,2,7,7)
grouptotal_Env_Spa_Bio <- cbind(site_scores_Env_Spa_Bio, group_Env_Spa_Bio)
factor(group_Env_Spa_Bio)
pdf("RDA_Env_Spa_Bio.pdf")
plot(new.rda_Env_Spa_Bio, display=c("lc","cn"),type="n")
for(i in 1:4){
  points(grouptotal_Env_Spa_Bio[group_Env_Spa_Bio==i,],pch=(14+i),cex=1.5)
}
for(i in 5:7){
  points(grouptotal_Env_Spa_Bio[group_Env_Spa_Bio==i,],pch=(i-5),cex=1.5)
}
points(grouptotal_Env_Spa_Bio[group_Env_Spa_Bio==8,1], grouptotal_Env_Spa_Bio[group_Env_Spa_Bio==8,2], pch=5, cex=1.5)
text(new.rda_Env_Spa_Bio, display=c("cn"), col="black")
dev.off()


##Variation partitioning with environmental factors(Temperature), Spatial variables (PCNM_V1), and seven biotic traits
varpart_Env_Spa_Bio <- varpart(bacteria, Env_signif, Spa_signif, Euk_signif)
pdf("VPA_Env_Spa_Bio.pdf")
plot(varpart_Env_Spa_Bio,digits=3)
dev.off()

#Test of fractions of environmental factors, spatial variables and biotic traits
anova.cca(rda(bacteria,Env_Spa_Bio[,1]), step = 1000)
anova.cca(rda(bacteria,Env_Spa_Bio[,3:9]), step = 1000)
anova.cca(rda(bacteria,Env_Spa_Bio[,2]), step = 1000)


#add Temporal component in varpart analysis
Temporal <- c("S","S","A","A","S","S","S","A","A","S","S","S","A","S","A","A","A","S","A","A")#Sping and Autumn
Temporal <- factor(Temporal,levels = c("S", "A"))
table(Temporal)
varpart_Env_Spa_Bio_Temporal <- varpart(bacteria, Env_signif,  Spa_signif, Temporal, Euk_signif)
plot(varpart_Env_Spa_Bio_Temporal,digits=3)


#add Temporal component in RDA plot
Env_Spa_Bio_Temporal <- data.frame(Env_signif, Spa_signif, Euk_signif, Temporal)#Create a data.fram containing Temperature, PCNM_V1, seven biotic traits and Temporal component
new.rda_Env_Spa_Bio_Temporal <- rda(bacteria ~ ., data=Env_Spa_Bio_Temporal)#construct new rda model
RsquareAdj(new.rda_Env_Spa_Bio_Temporal)$r.squared
RsquareAdj(new.rda_Env_Spa_Bio_Temporal)$adj.r.squared

#Calculate the value of AICc for Env_Spa_Bio_Temporal model
Env_Spa_Bio_Temporal_K <- length(new.rda_Env_Spa_Bio_Temporal)
AIC_Env_Spa_Bio_Temporal <- extractAIC(new.rda_Env_Spa_Bio_Temporal)[2]
AICc_Env_Spa_Bio_Temporal <- sum(AIC_Env_Spa_Bio_Temporal,2*Env_Spa_Bio_Temporal_K*(Env_Spa_Bio_Temporal_K+1)/(19-Env_Spa_Bio_Temporal_K))

permutest(new.rda_Env_Spa_Bio_Temporal, permu=999)
vif.cca(new.rda_Env_Spa_Bio_Temporal)

#Draw RDA plot
site_scores_Env_Spa_Bio_Temporal <- scores(new.rda_Env_Spa_Bio_Temporal, display="lc")
group_Env_Spa_Bio_Temporal <- c(1,3,6,6,1,4,3,7,5,3,2,2,7,4,6,8,5,2,7,7)
grouptotal_Env_Spa_Bio_Temporal <- cbind(site_scores_Env_Spa_Bio_Temporal, group_Env_Spa_Bio_Temporal)
factor(grouptotal_Env_Spa_Bio_Temporal)
plot(new.rda_Env_Spa_Bio_Temporal, display=c("lc","cn"),type="n")
for(i in 1:4){
  points(grouptotal_Env_Spa_Bio_Temporal[group_Env_Spa_Bio_Temporal==i,],pch=(14+i),cex=1.5)
}
for(i in 5:7){
  points(grouptotal_Env_Spa_Bio_Temporal[group_Env_Spa_Bio_Temporal==i,],pch=(i-5),cex=1.5)
}
points(grouptotal_Env_Spa_Bio_Temporal[group_Env_Spa_Bio_Temporal==8,1], grouptotal_Env_Spa_Bio_Temporal[group_Env_Spa_Bio_Temporal==8,2], pch=5, cex=1.5)
text(new.rda_Env_Spa_Bio_Temporal, display=c("cn"), col="black")

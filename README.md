# BioticInteraction

This repository has scripts and processed data for "Quantifying relative contributions of biotic interactions to bacterial diversity and community assembly using a trait-based approach"

Script explanations:
  
Mantel Test.R             - Find the Mantel statistic as a matrix correlation between two dissimilarity matrices from bacteria and pico-nano eukaryotes  
OTU Table Normalization.R - Normalize the raw OTU table of bacteria and pico-nano eukaryotes using the edgeR method  
PCNM analysis.R           - Perform the principal coordinates of neighbor matrices (PCNM) analysis to decompose the total spatial variation into a finite set of explanatory spatial variables  
RDA and VPA analysis_revised.R    - Identify the significant environmental variables, spatial factors and biotic traits to the variation in bacterial community structure, and quantify the relative contributions of  these three parts  
barplot.R                 - Barplot for relative abundance of major bacterial, pico-nano eukaryotic taxa and functional composition


data/

spatial_variables.csv 	- spatial variables generated from PCNM analysis  


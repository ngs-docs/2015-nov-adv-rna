#!/usr/bin/Rscript

library(ggplot2)
library("reshape2")
library(scales)

#Get the data into an R dataframe:

#set working directory --> path to your salmon quants folder
setwd("~/dib-training-materials/salmon_quants")

## list all quant.sf files:
files <- (Sys.glob("*/quant.sf"))

#read all quant.sf files in; save *only* the TPM column from each file.
salmonTPM <- do.call("cbind",lapply(files,read.csv,sep = "\t", skip=9, colClasses=c("NULL","NULL","numeric", "NULL"))) 

#name the TPM columns by the filename they came from: 
colnames(salmonTPM) <- sapply(strsplit(files, '.q'), '[', 1)

# All salmon output files are in the same order, so we can grab the gene names from one of the files; store them as a column:
rownames(salmonTPM)  <- read.csv(files[1],sep = "\t", skip=9, colClasses=c("character", "NULL","NULL","NULL"))[[1]]

#sum TPM
sumTPMPerContig <- rowSums(salmonTPM) #get summed TPM

#Plotting:
# -->add code for simple plot (coverage per contig: x-axis, count/frequency: y-axis) 

# Reshape the dataframe from wide --> long form for easier plotting with ggplot2
salmonTPM$Contig  <- read.csv(files[1],sep = "\t", skip=9, colClasses=c("character", "NULL","NULL","NULL"))[[1]]
meltedTPM <- melt(salmonTPM, id.vars="Contig",value.name="TPM", variable.name="Treatment") 

#plot contigs on x axis, coverage in each treatment on the y axis. WAY TOO NOISY/complicated with this much data... change into something simpler.
ggplot(data=meltedTPM,aes(x=Contig, y =TPM, color=Treatment)) + geom_point(size=.5) + scale_y_continuous(trans=log2_trans()) + theme(axis.ticks = element_blank(), axis.text.x  = element_blank()) + theme_bw()






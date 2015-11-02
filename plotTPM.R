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
salmonTPM$Contig  <- read.csv(files[1],sep = "\t", skip=9, colClasses=c("character", "NULL","NULL","NULL"))[[1]]

# Reshape the dataframe from wide --> long form for easier plotting with ggplot2
meltedTPM <- melt(salmonTPM, id.vars="Contig",value.name="TPM", variable.name="Treatment") 

#Plot a histogram of the TPM per Contig
ggplot(meltedTPM, aes(x=TPM)) + geom_histogram(binwidth=1) + xlab("Coverage Per Contig") + scale_x_continuous(limits = c(0, 200)) + scale_y_continuous(limits = c(0, 10000)) + theme_bw()



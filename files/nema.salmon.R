library("edgeR")

files <- c("0Hour_ATCACG_L002001.quant.counts",
           "0Hour_ATCACG_L002002.quant.counts",
           "0Hour_ATCACG_L002003.quant.counts",
           "0Hour_ATCACG_L002004.quant.counts",
           "0Hour_ATCACG_L002005.quant.counts",
           "12Hour_TTAGGC_L002001.quant.counts",
           "12Hour_TTAGGC_L002002.quant.counts",
           "12Hour_TTAGGC_L002003.quant.counts",
           "12Hour_TTAGGC_L002004.quant.counts",
           "18Hour_TGACCA_L002001.quant.counts",
           "18Hour_TGACCA_L002002.quant.counts",
           "18Hour_TGACCA_L002003.quant.counts",
           "18Hour_TGACCA_L002004.quant.counts",
           "18Hour_TGACCA_L002005.quant.counts",
           "18Hour_TGACCA_L002006.quant.counts",
           "18Hour_TGACCA_L002007.quant.counts",
           "18Hour_TGACCA_L002008.quant.counts",
           "24HourA_ACAGTG_L002001.quant.counts",
           "24HourA_ACAGTG_L002002.quant.counts",
           "24HourA_ACAGTG_L002003.quant.counts",
           "24HourA_ACAGTG_L002004.quant.counts",
           "24HourA_ACAGTG_L002005.quant.counts",
           "24HourA_ACAGTG_L002006.quant.counts",
           "24HourA_ACAGTG_L002007.quant.counts",
           "24HourA_ACAGTG_L002008.quant.counts",
           "24HourA_ACAGTG_L002009.quant.counts",
           "24HourA_ACAGTG_L002010.quant.counts",
           "24HourB_GCCAAT_L002001.quant.counts",
           "6Hour_CGATGT_L002001.quant.counts",
           "6Hour_CGATGT_L002002.quant.counts",
           "6Hour_CGATGT_L002003.quant.counts",
           "6Hour_CGATGT_L002004.quant.counts",
           "6Hour_CGATGT_L002005.quant.counts"
)

labels=c("0Hour_1", "0Hour_2", "0Hour_3", "0Hour_4", "0Hour_5",
"12Hour_1", "12Hour_2", "12Hour_3", "12Hour_4",
"18Hour_1", "18Hour_2", "18Hour_3", "18Hour_4",
"18Hour_5", "18Hour_6", "18Hour_7", "18Hour_8",
"24HourA_1", "24HourA_2", "24HourA_3", "24HourA_4",
"24HourA_5", "24HourA_6", "24HourA_7", "24HourA_8",
"24HourA_9", "24HourA_10", "24HourB_11",
"6Hour_1", "6Hour_2", "6Hour_3",
"6Hour_4", "6Hour_5")

data <- readDGE(files)

print(data)
head(data$counts)

###

group <- c(rep("0Hour",5), rep("12Hour",4), rep("18Hour", 8),
          rep("24Hour", 11), rep("6Hour", 5))

dge = DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# make an MA-plot of 0 vs 6 hour

et <- exactTest(dge, pair=c("0Hour", "6Hour"))
etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC
pdf("nema-edgeR-MA-plot.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "red", "black" ) )
dev.off()

# plot MDS
pdf("nema-edgeR-MDS.pdf")
plotMDS(dge, labels=labels)
dev.off()

# output CSV for 0-6 hr
write.csv(etp$table, "nema-edgeR.csv")


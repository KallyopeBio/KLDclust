library(ggplot2)
library(GGally)
library(reshape2)
library(hexbin)
library(gplots)
library(LDheatmap)
library(dplyr)

# Set input variables
rangeStart <- 22000000
rangeStop <-  24000000

# Read in correlation data from PLINK file
cc <- rep('NULL', 3)
cc[c(1, 2)] <- 'integer'
cc[c(3)] <- 'numeric'
corrTable <- read.table("/home/sarah/Data/1000genomes_phase3/Fisher_Penalty/MAF05/chr21/AFR/Correlations/all.txt", sep=" ", header=FALSE, colClasses=cc, nrows=2000000000)
colnames(corrTable) <- c("BP_A", "BP_B", "R")
print("Loading complete")

# Clean variables in table
corrTable$R <- abs(as.numeric(corrTable$R))
##corrTable$R <- replace(corrTable[["R"]], corrTable[["R"]]>.0001, .0001)
##corrTable$R <- corrTable$R/.0001
print("Data clean")

# Subset correlation data (NOTE RANDOM SAMPLING)
corrTable <- subset(corrTable, BP_A >= rangeStart & BP_A <= rangeStop & BP_B >= rangeStart & BP_B <= rangeStop & R>0.0001)
##corrTable <- sample_n(corrTable, 10000)
print("Sampling complete")

# Convert data to matrix

# get number of rows in data
ndata = dim(corrTable)[1]

# extract variant names
names1 = unique(corrTable$BP_A)
names2 = unique(corrTable$BP_B)
rownames=sort(unique(c(names1,names2)))
nrows = length(rownames)
print(nrows)

# initialize matrix
  out_matrix = matrix(NA, 
                      nrow = nrows, ncol = nrows,
                      dimnames=list(rownames, rownames))

  # iterate over input data rows
  for (i1 in 1:ndata) {
    # get matrix-row and matrix-column indices for the current data-row
    iR = which(rownames==corrTable[["BP_A"]][i1])
    iC = which(rownames==corrTable[["BP_B"]][i1])
    out_matrix[iR, iC] = corrTable[["R"]][i1]
  }

# Read in list of Pickrell blocks
#EUR_blockPTable <- read.table("/home/sarah/Data/testing/pickrell/output/pickrell/EUR_fourier_ls-chr21.txt")
#colnames(EUR_blockPTable) <- c("chrx", "blockStart", "blockStop")
#EUR_blockPTable <- subset(EUR_blockPTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
# Find number of Pickrell blocks 
#EUR_nPBlocks = dim(EUR_blockPTable)[1]

# Read in list of our blocks
AFR_blockKTable <- read.table("/home/sarah/Data/1000genomes_phase3/Fisher_Penalty/MAF05/chr21/AFR/clustering/FINAL_optBreaks.txt", sep=" ")
colnames(AFR_blockKTable) <- c("blockStart", "blockStop")
AFR_blockKTable <- subset(AFR_blockKTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
AFR_blockKTable <- subset(AFR_blockKTable, blockStart != blockStop) 
AFR_nKBlocks = dim(AFR_blockKTable)[1]

# Read in list of our blocks
SAS_blockKTable <- read.table("/home/sarah/Data/1000genomes_phase3/Mod_Penalty/MAF05/chr21/AFR/clustering/FINAL_optBreaks.txt", sep=" ")
colnames(SAS_blockKTable) <- c("blockStart", "blockStop")
SAS_blockKTable <- subset(SAS_blockKTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
SAS_blockKTable <- subset(SAS_blockKTable, blockStart != blockStop)
SAS_nKBlocks = dim(SAS_blockKTable)[1]

# Read in list of our blocks
#EAS_blockKTable <- read.table("/home/sarah/Data/1000genomes_phase3/MAF10/chr21/EAS/clustering/FINAL_optBreaks.txt", sep=" ")
#colnames(EAS_blockKTable) <- c("blockStart", "blockStop")
#EAS_blockKTable <- subset(EAS_blockKTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
#EAS_blockKTable <- subset(EAS_blockKTable, blockStart != blockStop)
#EAS_nKBlocks = dim(EAS_blockKTable)[1]

# Read in list of our blocks
#EUR_blockKTable <- read.table("/home/sarah/Data/1000genomes_phase3/MAF10/chr21/EUR/clustering/FINAL_optBreaks.txt", sep=" ")
#colnames(EUR_blockKTable) <- c("blockStart", "blockStop")
#EUR_blockKTable <- subset(EUR_blockKTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
#EUR_blockKTable <- subset(EUR_blockKTable, blockStart != blockStop)
#EUR_nKBlocks = dim(EUR_blockKTable)[1]

# Read in list of our blocks
#ALL_blockKTable <- read.table("/home/sarah/Data/1000genomes_phase3/MAF10/chr21/MAX/clustering/FINAL_optBreaks.txt", sep=" ")
#colnames(ALL_blockKTable) <- c("blockStart", "blockStop")
#ALL_blockKTable <- subset(ALL_blockKTable, blockStart >= rangeStart & blockStart <= rangeStop & blockStop >= rangeStart & blockStop <= rangeStop)
#ALL_blockKTable <- subset(ALL_blockKTable, blockStart != blockStop)
#ALL_nKBlocks = dim(ALL_blockKTable)[1]

# Plot heatmap of matrix
rgb.palette <- terrain.colors(256)
myHeatmap <- LDheatmap(out_matrix, genetic.distances=rownames, distances="physical", LDmeasure="r", color=rgb.palette, add.map=FALSE, title="Correlation Heatmap for Variants on Chromosome 21 with BP Locations between 22 and 25 million")
for (i in 1:AFR_nKBlocks) { 
    LDheatmap.highlight(myHeatmap, which.min(abs(AFR_blockKTable[["blockStart"]][i]-rownames)), which.min(abs(AFR_blockKTable[["blockStop"]][i]-rownames)), col="red", lwd=1)}
for (i in 1:SAS_nKBlocks) { 
    LDheatmap.highlight(myHeatmap, which.min(abs(SAS_blockKTable[["blockStart"]][i]-rownames)), which.min(abs(SAS_blockKTable[["blockStop"]][i]-rownames)), col="purple", lwd=1)}
#for (i in 1:EAS_nKBlocks) { 
#    LDheatmap.highlight(myHeatmap, which.min(abs(EAS_blockKTable[["blockStart"]][i]-rownames)), which.min(abs(EAS_blockKTable[["blockStop"]][i]-rownames)), col="pink", lwd=1)}
#for (i in 1:EUR_nKBlocks) { 
#    LDheatmap.highlight(myHeatmap, which.min(abs(EUR_blockKTable[["blockStart"]][i]-rownames)), which.min(abs(EUR_blockKTable[["blockStop"]][i]-rownames)), col="green", lwd=1)}
#for (i in 1:ALL_nKBlocks) { 
#    LDheatmap.highlight(myHeatmap, which.min(abs(ALL_blockKTable[["blockStart"]][i]-rownames)), which.min(abs(ALL_blockKTable[["blockStop"]][i]-rownames)), col="blue", lwd=1)}

#for (i in 1:EUR_nPBlocks) { 
#    LDheatmap.highlight(myHeatmap, which.min(abs(EUR_blockPTable[["blockStart"]][i]-rownames)), which.min(abs(EUR_blockPTable[["blockStop"]][i]-rownames)), col="red", lwd=1)}

# Plot heatmap
#d <- ggplot(corrTable, aes(BP_A, BP_B, z = RT))
#d + stat_summary_2d(fun="mean", bins=c(10000,10000))
#geom_hline(yintercept=19480375, color = "red")+geom_vline(xintercept=19480375, color = "red")+
#geom_hline(yintercept=20962339, color = "red")+geom_vline(xintercept=20962339, color = "red")+
#geom_hline(yintercept=22072764, color = "red")+geom_vline(xintercept=22072764, color = "red")+
#geom_hline(yintercept=23362875, color = "red")+geom_vline(xintercept=23362875, color = "red")+
#geom_hline(yintercept=22184972, color = "green")+geom_vline(xintercept=22184972, color = "green")





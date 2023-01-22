##################### Differential Expression with DESeq2 ################################
### Date: Tuesday 27 OCTOBER 2020 ############################################################ 
### By: Octavio Palacios ######################################################################
### Tags: Differential Expression/Chrm_A Zebra finch ######################################

# Input file: ChrA_counts.only.LTR_ERV.txt

#1.Setup:
#First, install DESeq2:
source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')

#Then load the libraries we will need into R:
library('DESeq2')
library('RColorBrewer')

#2.Read gene counts into a data frame
#Read sample gene counts a tab-delimited file. The rows of the data frame are genes while
#the columns are samples.
sampleNames <- c('Female.ovary_1Aligned.sortedByCoord.out.bam', 'Female.ovary_2Aligned.sortedByCoord.out.bam', 'Female.ovary_3Aligned.sortedByCoord.out.bam', 'Female.ovary_4Aligned.sortedByCoord.out.bam', 'Female.ovary_5Aligned.sortedByCoord.out.bam',
                 'Female.brain_1Aligned.sortedByCoord.out.bam', 'Female.brain_2Aligned.sortedByCoord.out.bam', 'Female.brain_3Aligned.sortedByCoord.out.bam', 'Female.brain_4Aligned.sortedByCoord.out.bam',
                 'Female_embryonic_brain_1Aligned.sortedByCoord.out.bam', 'Female_embryonic_brain_2Aligned.sortedByCoord.out.bam',
                 'Male.testis_1Aligned.sortedByCoord.out.bam', 'Male.testis_2Aligned.sortedByCoord.out.bam', 'Male.testis_3Aligned.sortedByCoord.out.bam')
filePath = '~/Desktop/Plotting/TPM/Emu/ChrA_counts.only.LTR_ERV.txt'
countData = read.table(file = filePath, header = TRUE, row.names = 1, sep = '\t')

# conver countData column to integer
#countData <- transform(countData, FO1=as.integer(FO1), FO2=as.integer(FO2), FO3=as.integer(FO3), FO4=as.integer(FO4), FO5=as.integer(FO5),
                       #FO6=as.integer(FO6), FO7=as.integer(FO7), FO8=as.integer(FO8), FO9=as.integer(FO9), FO10=as.integer(FO10),
                       #MT1=as.integer(MT1), MT2=as.integer(MT2), MT3=as.integer(MT3), MT4=as.integer(MT4), MT5=as.integer(MT5),
                       #MT6=as.integer(MT6), MT7=as.integer(MT7), MT8=as.integer(MT8))

dim(countData)  #view number of rows and columns
View(countData) #view countData
new_countData <- countData[,c(6:10,17:19)] #selecting the colums 1 and from 7 to 10
View(new_countData) #view new_countData

#Now create a second data frame for sample information, such as the experimental
#condition that each sample belongs to:
condition <- c('Female', 'Female', 'Female','Female', 'Female', 'Male', 'Male', 'Male')  #vector of column names for the data frame
colData <- data.frame(row.names=colnames(new_countData), 
                      condition=factor(condition, levels=c('Female','Male')))
colData

#3.Run DESeq2
#First, create a DESeqDataSet by specifying the gene counts data frame, the sample
#information data frame and a design model:
dataset <- DESeqDataSetFromMatrix(countData = new_countData,
                                  colData = colData,
                                  design = ~condition)
dataset

# Pre-filtering (by removing rows in which there are no reads or nearly no read)
dataset <- dataset[ rowSums(counts(dataset)) > 1, ]

#Then run the DESeq2 algorithm and extract results for our two-class comparison:
dds <- DESeq(dataset)

result <- results(dds, contrast=c('condition','Male','Female'), alpha = 0.05)
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)
View(result)
​
#Export results to Desktop
filePath = '~/Desktop/Emu/result.txt'
write.table(result, file = filePath, row.names = TRUE, col.names = TRUE)
​
#4.View results
#View a summary of DESeq2 results:
summary(result)
​
###### Convert read counts to transcripts per million 
d <- counts(dds, normalized = TRUE)
View(d)
​
# Usage: counts_to_tpm(countMat, geneLengths)
# Define counts_to_tpm
counts_to_tpm = function(countMat, geneLengths) {
  
  rpk = countMat / (geneLengths/1000)          # Reads per kilobase
  scalingFactors = colSums(rpk, na.rm=TRUE) / 10^6   # "Per million" scaling factor
  tpm = t( t(rpk) / scalingFactors)         # Transcripts per million
  
  return(tpm)
  
}
# Run counts_to_tpm
TPM <- counts_to_tpm(d, countData$Length)
View(TPM)
​
filePath = '~/Desktop/Plotting/TPM/Emu/ChrA_TPM_normalized.txt'
write.table(TPM, file = filePath, row.names = TRUE, col.names = TRUE, 
            sep = "\t")
​
#The top 150 up-regulated and down-regulated genes by p-value:
n = 20
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
View(topResults)

#Export topResults to Desktop
filePath = '~/Desktop/Emu/topResults_20.txt'
write.table(topResults, file = filePath, row.names = TRUE, col.names = TRUE)

#Plot log fold change vs. mean expression for all genes, with genes where p < 0.1 
#colored red:
plotMA(result, main='DESeq2: Male testis vs. Female ovary', ylim=c(-15,15))

#Male vs Female
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld)

#Plot counts for a single gene. Below is the plot for the gene with the lowest p-value:
plotCounts(dds, gene=which.min(result$padj), intgroup='condition', pch = 19)

#Display top genes normalized counts in a heatmap:
hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, 
        mar = c(13,5))

library(pheatmap)
#heatmap(as.matrix(nCounts),
nCounts <- counts(dds, normalized=TRUE)
Dat <- as.data.frame(nCounts[ row.names(topResults), ])
pheatmap(log(Dat+0.5), cluster_rows = TRUE, scale = "none")

# CytoReason - Exercise for bioinformatics candidates
# Submitted by Nomi Hadar
# Januar 2020 

library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(rstatix) #for wilcox_effsize
library(coin) #for wilcox_effsize



#######################################
# Definitions  
#######################################

#define input paths
COUNT_PATH = "data/counts.txt"
GENE_ANNOT_PATH = "data/gene-annotation.txt"
SAMPLE_ANNOT_PATH = "data/sample-annotation.txt"


#define output paths
OUT_FILE1 = "norm_log_CPM_counts.rda" #output binary file
OUT_FILE2 = "wilcox_test.txt" #output tab-separated file


#define constants
CUT_OFF = 0.75 #cut off - percent of samples 
MIN_CPM = 1 #minimu value of CPM

#CPM calculation: 
# divide by lib-size, and multiply by 10^6. 
# If log- apply log2.

#######################################
# Load Data 
#######################################

#load data
counts = read.table(COUNT_PATH, header= TRUE, row.names= 1)
sample_annot = read.delim(SAMPLE_ANNOT_PATH, header= TRUE)
gene_annot = read.delim(GENE_ANNOT_PATH, sep="\t", header= TRUE) 
 
#show dimentions 
dim(counts) #[57992, 178]
dim(sample_annot) #[178, 2]
dim(gene_annot) #[25503, 4]

#get number of genes 
n_genes = dim(counts)[1]

#get groups sizes
n_normal = sum(sample_annot$type=='normal') #83 normal samples 
n_lesional = sum(sample_annot$type=='lesional') #95 lesional samples 
n_normal+n_lesional == dim(sample_annot)[1]
  
#make sure the count and annotation data are consistent with each other
#by ordering counts and annotation the samples ids 
counts = counts[,order(names(counts))]
sample_annot = sample_annot[order(sample_annot$sample_id),]

#######################################
# Filter Data 
#######################################

#Filter the count data for lowly-expressed genes
#only keep genes  with a CPM >= 1 in at least 75% samples, 
#in at least one of the groups.

filterCounts <- function(row, sample_annot){
  
  #get normal and lesional samples
  normal = row[sample_annot[sample_annot$type=='normal',]$sample_id]
  lesional = row[sample_annot[sample_annot$type=='lesional',]$sample_id]
  
  #filter by CPM value and a cut off
  (mean(normal >= MIN_CPM) >= CUT_OFF) | (mean(lesional >= MIN_CPM) >= CUT_OFF)
  
}

#compute CPM: Count-Per-Million 
counts_cpm = cpm(counts)

#genes to keep for the normal and lesional groups
keep = apply(counts_cpm, 1, filterCounts, sample_annot=sample_annot)

#filter counts 
counts_filtered = counts[keep,]
dim(counts_filtered) #[18286, 178]

#39,706 genes were removed
n_genes - dim(counts_filtered)[1] 


#######################################
# Normalize and compute log-CPM
#######################################

#noramlize by total number of counts per sample
#and compute log-CPM
norm_log_cpm = cpm(counts_filtered, log=TRUE, normalized.lib.sizes=TRUE)

#save to a binary file 
save(norm_log_cpm ,file=OUT_FILE1)


#######################################
# Investigate main properties
#######################################

#-------------------------
# Histograms of lib sizes
#-------------------------

#get lib sizes 
lib_sizes = apply(norm_log_cpm, 2, sum)

#plot histograms of lib sizes: for all samples, and per group
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(lib_sizes, main="Hist of all samples")
hist(lib_sizes[normal_samples], main="Hist of normal samples", xlab="Normal")
hist(lib_sizes[lesional_samples], main="Hist of lesional samples", xlab="Lesional")

#-------------------------
# PCA
#-------------------------

#compute pca object
pca.obj = prcomp(t(norm_log_cpm), center = TRUE,scale. = TRUE)

summary(pca.obj);

#create a dataframe of first and second components + groups
dtp = data.frame('Group' = sample_annot$type, pca.obj$x[,1:2])

#plot pca
ggplot(data = dtp) + 
  geom_point(aes(x= PC1, y= PC2, col = Group)) 
  #+geom_text(aes(x= PC1, y= PC2, label = 1:nrow(dtp))) #uncomment to get points labeld 

#-------------------------
# Remove samples
#-------------------------

#remove outlier/mis-labeled samples from the downstream analysis

outlier_idx = 3 #SPR1146078, probably an outlier
mis_idx = 141 #SRR1146216, probably a mis-labeled 
col_to_drop = noquote(rownames(dtp[c(outlier_idx, mis_idx), ]))

#remove from both normalized and not-normalized matrices

counts_filtered = counts_filtered[,!(colnames(counts_filtered) %in% col_to_drop)]
norm_log_cpm = norm_log_cpm[,!(colnames(norm_log_cpm) %in% col_to_drop)]

sample_annot = sample_annot[!(sample_annot$sample_id %in% col_to_drop),]

dim(counts_filtered)
dim(norm_log_cpm)
dim(sample_annot)

#######################################
# Differential expression analysis
#######################################

#-------------------------
# Wilcoxon Rank Sum
#-------------------------

myWilcoxTest <- function(row, sample_annot){
  
  "
  Wilcoxon Rank Sum test between normal and lesional
  "
  
  #get normal and lesional samples
  x = row[sample_annot[sample_annot$type=='normal',]$sample_id]
  y = row[sample_annot[sample_annot$type=='lesional',]$sample_id]
  
  #compute Wilcoxon Rank Sum test 
  t = wilcox.test(x, y, alternative = "two.sided")
  
  #return statistic value and and p-value  
  c(t$statistic, t$p.value)
  
}

#compute Wilcox-Rank-Sum Test for each gene
WilcoxTestResult = apply(counts_filtered, 1, myWilcoxTest, sample_annot=sample_annot)

#bind Wilcox test's statistics and p-values into one dataframe
Wilcox.statistic = WilcoxTestResult[1,]
pvalue = WilcoxTestResult[2,]
Wilcox_test = as.data.frame(cbind(Wilcox.statistic, pvalue))

#-------------------------
# Adjust P-values for Multiple Comparisons
#-------------------------

#get p-values adjusted, using FDR method
Wilcox_test$FDR = p.adjust(Wilcox_test$pvalue, method = "fdr")

head(Wilcox_test)

#-------------------------
# Save Wilcoxon result to file
#-------------------------

#filter genes annotation 
gene_annot_filtered = gene_annot[gene_annot$ENSEMBL, ]

#set rows names as a column to be able to merge
Wilcox_test$ENSEMBL = row.names(Wilcox_test)

#merge statistics and gene annotation
Wilcox_merged = merge(x= Wilcox_test, y= gene_annot , by= "ENSEMBL", all.x= TRUE)

dim(Wilcox_merged) #[18286, 6]
head(Wilcox_merged)

#sort by p-value 
Wilcox_merged = Wilcox_merged[with(Wilcox_merged, order(pvalue)), ]

#save to file 
write.table(Wilcox_merged, OUT_FILE2, sep="\t", row.names= FALSE, quote= FALSE)

#-------------------------
# get significant genes
#-------------------------

N_TOP = 100

#Select the top 100 most significant annotated genes
top_genes = Wilcox_merged[!is.na(Wilcox_merged$ENTREZID ),][1:N_TOP,]$ENSEMBL

#select counts of top genes 
top_lcpm = norm_log_cpm[top_genes,]
dim(top_lcpm)

#-------------------------
# Plot Heatmap
#-------------------------

#set sample IDs as rows names (needed for heatmap)
rownames(sample_annot) = sample_annot$sample_id

#plot heatmap 
pheatmap(top_lcpm, annotation_col = sample_annot["type"], 
         fontsize = 6, main= "Top 100 significant genes")


#-------------------------
# Volcano plot
#-------------------------

par(mfrow=c(1,1))
#with(Wilcox_test, plot(Wilcox_test, -log10(FDR),  main="Volcano plot"))
#, colour=threshold


myWilcoxEffectSize <- function(x, sample_annot){
  "
  Effect size of Wilcoxon Rank Sum test 
  between normal and lesional
  "
  #get samples annotations  
  annot = sample_annot$type
 
  #create a df of count values and sample annotation  
  data = data.frame(cbind(x, annot))
  
  #compute Wilcoxon Rank effect size
  effsize = wilcox_effsize(data, x ~ annot, alternative = "two.sided")
  
  #extract value of effect size 
  effsize$effsize
  
}


#compute Wilcox-Rank-Sum Test for each gene
#takes serveral minutes =(
Wilcox_effectsize = apply(counts_filtered, 1, myWilcoxEffectSize, 
                          sample_annot=sample_annot) 

#to df
Wilcox_effectsize_df = data.frame(Wilcox_effectsize)

#apply log 2 
log2effectSize = log2(Wilcox_effectsize)

#set genes as a column
Wilcox_effectsize_df["ENSEMBL"] = row.names(Wilcox_effectsize_df)

#merge adjsted Wilcox test p-values and effect sizes
merged = merge(x=Wilcox_test , y= Wilcox_effectsize_df, by="ENSEMBL")

#mark true top genes
rownames(merged) <- merged$ENSEMBL
merged["is_top"] = FALSE
merged$is_top[which(merged$ENSEMBL %in% top_genes)] = TRUE 


#plot volcano
ggplot(merged) +
  geom_point(aes(x=log2(Wilcox_effectsize), y=-log10(FDR),  colour = factor(is_top)   )) +

  geom_point(data = subset(merged, is_top == TRUE),
          aes(x=log2(Wilcox_effectsize), y=-log10(FDR),  colour = factor(is_top)   )) +
  
  
  ggtitle("Volcano plot") +
  xlab("log2 Wilcox effect size") + 
  ylab("-log10 adjusted p-value (FDR)") +
  
  
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
  







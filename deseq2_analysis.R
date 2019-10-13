library("DESeq2")
library(RColorBrewer)
try(library(gplots))
library("lattice")
library("readr")
library("plyr")

arguments = commandArgs(trailingOnly = T)
print(arguments)

if (length(arguments)!=2) {
  stop("You should provide sample description file and count matrix.")
}

#sample description is the first argument
samples <- read.table(arguments[1], sep = "\t", h = T)

#count matrix is the second argument

countdf = read.table(arguments[2], sep = "\t", h = T, stringsAsFactors = F)

#set working directory to the directory of count matrix
setwd(dirname(arguments[2]))

#annotation for GRCz11

ann <- read.table("/ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/gene_names_locations_mart_export_nozfin.txt", 
                  header = T, sep = "\t", stringsAsFactors = F, strip.white = T, fill = TRUE, quote="")

colnames(ann)[1] = "GeneID"


#prepare matrix

rownames(countdf) = countdf$GeneID

countdf = countdf[,-1]


#analyze all samples
dds_counts <- DESeqDataSetFromMatrix(countData = countdf,
                                     colData = samples,
                                     design = ~ condition)
dds <- DESeq(dds_counts)

dir.create("whole-matrix-output")

#save normalized counts for all samples with annotations
counts_df <- as.data.frame(counts(dds, normalized=TRUE))
norm_counts <- cbind(GeneID = row.names(counts_df), counts_df, stringsAsFactors = F)
norm_counts_ann <-  merge(norm_counts, ann, by = "GeneID", all.norm_counts = T)
write.table(norm_counts_ann, "whole-matrix-output/normalized_counts_deseq.tsv", sep = "\t", row.names = F, quote = F, dec = ".")

#produce sample pairs, factors to strings
sample_vect = unique(as.character(samples$condition))
combinations <-combn(sample_vect, 2, simplify = F)
#write combinations into file
for (i in c(1:length(combinations))){
  write(paste(combinations[[i]][[1]], combinations[[i]][[2]], sep= "-"), "whole-matrix-output/combinations.txt", append = T)
}





#pairwise DE analysis, with all the samples previously normalized together
deseq_pair_result <- function(sample1, sample2, dds){
  fname = paste0("whole-matrix-output/", sample1, "-", sample2, "-whole-matrix.tsv")
  fname_sign = paste0("whole-matrix-output/", sample1, "-", sample2, "-whole-matrix-sign.tsv")
  res <- results(dds, contrast=c("condition", sample2, sample1))
  res_df <- cbind(as.data.frame(res), GeneID = row.names(res), stringsAsFactors = F)
  res_df_ann <- merge(res_df, ann, by = "GeneID", all.res_df = T)
  res_df_ann_pval <- res_df_ann[res_df_ann$padj < 0.05 & !is.na(res_df_ann$padj),]
  write.table(res_df_ann, fname, sep = "\t", row.names = F, quote = F, dec = ".")
  write.table(res_df_ann_pval, fname_sign, sep = "\t", row.names = F, quote = F, dec = ".")
  return(res_df_ann_pval)
}

for (i in c(1:length(combinations))){
  sample1 = combinations[[i]][[1]]
  sample2 = combinations[[i]][[2]]
  deseq_pair_result(sample1, sample2, dds)
}

#visualize the whole dds matrix
rld <- rlogTransformation(dds)

#pca
png(paste0("whole-matrix-output/all_samples_", "lab_pca_1x2.png"), w=800, h=800, pointsize=8)
print(plotPCA(rld, intgroup="condition"), ntop = 500, returnData = T)
dev.off()

#visualizations

condition <- samples$condition
sampname = "whole-matrix-output/all_samples_matrix"
png(paste0(sampname, "_qc-dispersions.png"), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

head(assay(rld))
#hist(assay(rld))
mycols <- c(brewer.pal(12, "Set3"), brewer.pal(4, "Dark2"))[1:length(unique(samples$condition))]
sampleDists <- as.matrix(dist(t(assay(rld))))


png(paste0(sampname, "_heatmap.png"), w=1200, h=1000, pointsize=20)
try(heatmap.2(as.matrix(sampleDists), key=F, trace="none",
              col=colorpanel(100, "black", "white"),
              ColSideColors=mycols[condition], RowSideColors=mycols[condition],
              margin=c(10, 10), main="Sample Distance Matrix"))
#legend("right", title = "conditions",legend=condition, 
#       fill=mycols, cex=0.5, box.lty=0)
dev.off()

print("Differential gene expression analysis complete")

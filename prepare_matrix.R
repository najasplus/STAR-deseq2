
#merge individual count files into one matrix

#find all matrices 
all_files <- list.files(path = "gene_counts", pattern = ".count")
file1 <- read.table(paste0("gene_counts/", all_files[1]), sep = "\t", h = F)
gnames <- data.frame(matrix(nrow = nrow(file1), ncol = 0), stringsAsFactors = F)


#check if all count files have the same length and same gene order
for (f in all_files){
  f_table <- read.table(paste0("gene_counts/", f), sep = "\t", h = F)
  gnames <- data.frame(gnames, f_table[,1])
  ifelse(gnames[,1] == f_table[,1], TRUE, stop("The order of genes is not the same"))
  }

#set up a matrices for each library type 
# col2 unstranded
# col3, stranded forward
# col4 (=htseq -s rev), stranded, reverse

col4_df <- data.frame(matrix(nrow = nrow(gnames),ncol = 0))
col3_df <- data.frame(matrix(nrow = nrow(gnames),ncol = 0))
col2_df <- data.frame(matrix(nrow = nrow(gnames),ncol = 0))

sampnames <- c()

for (f in all_files){
  sample_name <- unlist(strsplit(f, split = ".count"))
  f_table <- read.table(paste0("gene_counts/", f), sep = "\t", h = F)
  sampnames = c(sampnames, sample_name)
  col2_df <- data.frame(col2_df, f_table[,2])
  col3_df <- data.frame(col3_df, f_table[,3])
  col4_df <- data.frame(col4_df, f_table[,4])
}

#attach gene IDs and sample names as column names, write matrixes into file
df_names = c("col2_df", "col3_df", "col4_df")

for (df_name in df_names) {
  df = get(df_name)
  df = data.frame(gnames[,1], df)
  colnames(df) = c("GeneID", sampnames)
  write.table(df, paste0(df_name, "_raw_counts.tsv"), sep = "\t", row.names = F, col.names = T)
}


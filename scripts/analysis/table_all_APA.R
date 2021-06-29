# Script for creating binary table of significant APA events 
# for all comparisons

library(biomaRt)
library(dplyr)

# Load tables

files <- c('results/luad_matched/lambrechts/auc_comparisons/analysis_out/cpm2_TEs_shorter_IQR_cancer_vs_alveolar_normal.tsv',
           'results/luad_matched/laughney/auc_comparisons/analysis_out/cpm2_TEs_shorter_IQR_cancer_vs_alveolar_normal.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_alveolar.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_B_cell.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_endothelial.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_fibroblast.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_mast_cell.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_myeloid.tsv',
           'results/luad_matched/lambrechts/auc/analysis_out/cpm2_TEs_shorter_IQR_T_cell.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_alveolar.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_B_cell.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_endothelial.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_fibroblast.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_mast_cell.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_myeloid.tsv',
           'results/luad_matched/laughney/auc/analysis_out/cpm2_TEs_shorter_IQR_T_cell.tsv',
           'results/luad_matched/patient3/auc/analysis_out/cpm2_TEs_shorter_IQR_endothelial.tsv',
           'results/luad_matched/patient3/auc/analysis_out/cpm2_TEs_shorter_IQR_fibroblast.tsv',
           'results/luad_matched/patient3/auc/analysis_out/cpm2_TEs_shorter_IQR_myeloid.tsv',
           'results/luad_matched/patient3/auc/analysis_out/cpm2_TEs_shorter_IQR_T_cell.tsv',
           'results/luad_matched/patient4/auc/analysis_out/cpm2_TEs_shorter_IQR_endothelial.tsv',
           'results/luad_matched/patient4/auc/analysis_out/cpm2_TEs_shorter_IQR_fibroblast.tsv',
           'results/luad_matched/patient4/auc/analysis_out/cpm2_TEs_shorter_IQR_myeloid.tsv',
           'results/luad_matched/patient4/auc/analysis_out/cpm2_TEs_shorter_IQR_T_cell.tsv',
           'results/luad_matched/patient6/auc/analysis_out/cpm2_TEs_shorter_IQR_endothelial.tsv',
           'results/luad_matched/patient6/auc/analysis_out/cpm2_TEs_shorter_IQR_fibroblast.tsv',
           'results/luad_matched/patient6/auc/analysis_out/cpm2_TEs_shorter_IQR_myeloid.tsv',
           'results/luad_matched/patient6/auc/analysis_out/cpm2_TEs_shorter_IQR_T_cell.tsv')

names(files) <- c('cancer_LB', 'cancer_LN', 
                  'alveolar_LB', 'B_cell_LB', 'endothelial_LB', 'fibroblast_LB', 'mast_cell_LB', 'myeloid_LB', 'T_cell_LB',
                  'alveolar_LN', 'B_cell_LN', 'endothelial_LN', 'fibroblast_LN', 'mast_cell_LN', 'myeloid_LN', 'T_cell_LN',
                  'endothelial_p3', 'fibroblast_p3', 'myeloid_p3', 'T_cell_p3',
                  'endothelial_p4', 'fibroblast_p4', 'myeloid_p4', 'T_cell_p4',
                  'endothelial_p6', 'fibroblast_p6', 'myeloid_p6', 'T_cell_p6')

tbls <- sapply(files, function(x) read.table(x, sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,1])

# Obtain all TE ids and create table for conversion
new_config <- httr::config(ssl_verifypeer = FALSE)
httr::set_config(new_config, override = FALSE)

conv_ttog <- function(myvector, mytable){
  res <- sapply(myvector, function(x) mytable[mytable[,1] == x, 2][1])
  return(unname(res))
}

# select proper organism!
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
transcripts <- unique(unlist(tbls))
transcript.tbl <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),
                           values = transcripts, mart= mart)

tbls.genes <- lapply(tbls, function(x) conv_ttog(x, transcript.tbl))

# omit NA and only unique genes
tbls.genes <- lapply(tbls.genes, function(x) unique(na.omit(x)))

# Construct complete table
df <- data.frame(genes = unique(unlist(tbls.genes)))

for(i in 1:length(tbls.genes)){
  genes <- tbls.genes[[i]]
  binary_vec <- vector(mode = "logical", length = nrow(df))
  for(gene in genes){
    for(k in 1:nrow(df)){
      if(df$genes[k] == gene){
        binary_vec[k] = TRUE
      }
    }
  }
  df[names(tbls.genes)[i]] <- binary_vec
}

# Display some statistics
## Distribution of how many TRUEs per gene
df %>% mutate()
df$n_TRUE <- apply(df[,-1], 1, function(x) length(which(unlist(x))))

hist(df$n_TRUE, breaks = 30)

# Convert logical to numeric
df.num <- df %>% mutate_if(is.logical, as.numeric)

# re-order table
df.num <- arrange(df.num, -n_TRUE)

# Save table
write.table(df.num,
            file = 'results/luad_matched/comparison/table_all_APA.tsv',
            sep = "\t", quote = FALSE,
            row.names = FALSE)


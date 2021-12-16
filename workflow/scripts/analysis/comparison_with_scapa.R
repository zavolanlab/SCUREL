# Comparison of number of genes found significant in
# SCUREL and
# scAPA
# for activated T cells, spermatocytes and LUAD datasets

setwd('~/polyasite')

indir.auc <- 'results/spermatocytes/20210421_A/auc_comparisons'
prefix <- 'cpm2_'
ct <- 'es_vs_sc'

indir.scapa <- 'scAPA/spermatocytes/scAPA/outs'

outdir <- 'comparison_scAPA/spermatocytes_cpm2'

# load list of genes
df.auc.all <- read.table(file.path(indir.auc, paste0('coverage_3p/auc_', ct, '.tsv')),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df.auc <- list()
df.auc[1] <- read.table(file.path(indir.auc, paste0('analysis_out/', prefix, 'TEs_shorter_IQR_', ct, '.tsv')),
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
df.auc[2] <- read.table(file.path(indir.auc, paste0('analysis_out/', prefix, 'TEs_longer_IQR_', ct, '.tsv')),
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(df.auc) <- c('shortening', 'lengthening')

df.scapa.all <- read.table(file.path(indir.scapa, 'genes_sig_two.tsv'),
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df.scapa <- list()
df.scapa[1] <- read.table(file.path(indir.scapa, 'genes_shortening_ct1.txt'),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
df.scapa[2] <- read.table(file.path(indir.scapa, 'genes_lengthening_ct1.txt'),
                          sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(df.scapa) <- c('shortening', 'lengthening')

# transform gene names to external gene name
library(biomaRt)
new_config <- httr::config(ssl_verifypeer = FALSE)
httr::set_config(new_config, override = FALSE)

conv_ttog <- function(myvector, mytable){
  res <- sapply(myvector, function(x) mytable[mytable[,1] == x, 2][1])
  return(unname(res))
}

# for refseq transcripts in AUC
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
transcripts <- df.auc.all$TE_id
mouse_transcripts <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),
               values = transcripts, mart= mart)

df.auc.all$gene_name <- conv_ttog(df.auc.all$TE_id, mouse_transcripts)
df.auc <- lapply(df.auc, function(x) conv_ttog(x, mouse_transcripts))

# for ensembl genes in scAPA
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- union(df.scapa$shortening, df.scapa$lengthening)
mouse_genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
                                 values = genes, mart= mart)


df.scapa.all$gene_name <- conv_ttog(df.scapa.all$gene, mouse_genes)
df.scapa <- lapply(df.scapa, function(x) conv_ttog(x, mouse_genes))

# compare the two sets
library(ggplot2)
library(ggpolypath)
venn::venn(c(df.auc, df.scapa),
           snames = c('auc.shortening', 'auc.lengthening', 'scapa.shortening', 'scapa.lengthening'),
           ggplot = TRUE)
ggsave(paste0(outdir, '_venn_all.png'),
       device = "png", units = "cm", width = 8, height = 8)

venn::venn(list(df.auc$shortening, df.scapa$shortening),
           snames = c('AUC', 'scAPA'),
           ggplot = TRUE)
ggsave(paste0(outdir, '_venn_shortening.png'),
       device = "png", units = "cm", width = 5, height = 5)

venn::venn(list(df.auc$lengthening, df.scapa$lengthening),
           snames = c('AUC', 'scAPA'),
           ggplot = TRUE)
ggsave(paste0(outdir, '_venn_lengthening.png'),
       device = "png", units = "cm", width = 5, height = 5)



# check PPUI usage
library(dplyr)
df.all <- full_join(df.auc.all, df.scapa.all,
          by = "gene_name")
df.all <- mutate(df.all, TE_length = TE_end - TE_start)

# add column span
get_span <- function(values){
  # values is a vector of length 4, containing 25% and 75% 3'UTR length for both samples
  # obtain IQR span from sample1 and sample2 (in base pairs)
  return(max(values) - min(values))
}

df.all <- mutate(df.all, sample1_IQR = sample1_length_q0.25 - sample1_length_q0.75,
                 sample2_IQR = sample2_length_q0.25 - sample2_length_q0.75)
df.all$SPAN <- apply(dplyr::select(df.all, sample1_length_q0.25, sample1_length_q0.75,
                                   sample2_length_q0.25, sample2_length_q0.75),
                     1,
                     get_span)

# df.all <- na.omit(df.all, cols = c(diff))
# genes sign. in both
df.tmp <- filter(df.all, gene_name %in% na.omit(intersect(df.auc$shortening, df.scapa$shortening)))  %>% 
  dplyr::select(TE_id, gene_name, coverage_sample1, cpm_sample1, coverage_sample2, cpm_sample2, auc, diff, SPAN) %>% 
  filter(cpm_sample1 >= 2, cpm_sample2 >= 2, SPAN > 200) %>% arrange(-diff)
write.table(df.tmp,
            file = paste0(outdir, '_genes_sig_both.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

# genes sign. in scapa
df.tmp <- filter(df.all, gene_name %in% setdiff(df.scapa$shortening, df.auc$shortening)) %>% 
  dplyr::select(TE_id, gene_name, TE_length, coverage_sample1, cpm_sample1, coverage_sample2, cpm_sample2, auc, diff, SPAN) %>% 
  arrange(cpm_sample1) %>% arrange(cpm_sample2) %>% arrange(-diff)
write.table(df.tmp,
            file = paste0(outdir, '_genes_sig_scapa.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

# genes sign. in AUC
df.tmp <- filter(df.all, gene_name %in% na.omit(setdiff(df.auc$shortening, df.scapa$shortening))) %>% 
  dplyr::select(TE_id, gene_name, TE_length, coverage_sample1, cpm_sample1, coverage_sample2, cpm_sample2, auc, diff, SPAN) %>% 
  filter(cpm_sample1 >= 2, cpm_sample2 >= 2, SPAN > 200) %>% arrange(-auc, -SPAN)
write.table(df.tmp,
            file = paste0(outdir, '_genes_sig_auc.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)


  

#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { stop("Package 'optparse' required!\nExecution aborted.") }

## ------------------------------------------------------------------------
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

## ------------------------------------------------------------------------
description <- "Perform AUC analysis. \n"
author <- "Author: Dominik Burri"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: dominik.burri@unibas.ch"
version <- "Version: 1.0.0 (09-FEB-2021)"
requirements <- c("optparse", "dplyr", "ggplot2")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version, requirements_txt, sep="\n")
notes <- "Notes: No additional notes."

## ------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--analysis-dir",
    action = "store",
    type = "character",
    default = NULL,
    help = "Directory to analysis dir, must contain coverage_3p and randomised_coverage subdirectories.",
    metavar = "dir"
  ),
  make_option(
    "--out-dir",
    action = "store",
    type = "character",
    default = "analysis_out",
    help = "Define output directory. In the supplied directory all results will be written. Default: analysis_out",
    metavar = "dir"
  ),
  make_option(
    "--comparison",
    action = "store",
    type = "character",
    default = NULL,
    help = "Either cell type to analyse (for execution mode cell_state_comparison), 
      respectively comparison to analyse (for execution mode cell_type_comparison, e.g. c1). 
      Default: NULL (no execution)",
    metavar = "char"
  ), 
  make_option(
    "--cpm-threshold",
    action = "store",
    type = "integer",
    default = 10,
    help = "Define minimum expression threshold (Counts per million, CPM) for any given Terminal exon. Default: 10.",
    metavar = "int"
  ),
  make_option(
    "--read-threshold",
    action = "store",
    type = "integer",
    default = NULL,
    help = "Define minimum expression threshold in reads for any given Terminal exon. If this threshold is given, cpm-threshold is ignored. Default: No read threshold.",
    metavar = "int"
  ),
  make_option(
    "--length-threshold",
    action = "store",
    type = "double",
    default = 0,
    help = "Define minimum length of Terminal Exon based on lower quantile. Default: 0 (No filtering on TE length performed)",
    metavar = "double"
  ),
  make_option(
    "--iqr-threshold",
    action = "store",
    type = "integer",
    default = NULL,
    help = "Define minimum interquantile range (IQR) required to separate 1 PAS and multi PAS 3'UTRs. Requires quantile columns created by calculation_auc.py. Default: NULL",
    metavar = "int"
  ),
  make_option(
    "--ALPHA",
    action = "store",
    type = "double",
    default = 0.05,
    help = "Two-tailed probability for considering Terminal exons significant. Default: 0.05",
    metavar = "double"
  ),
  make_option(
    "--out-prefix",
    action = "store",
    type = "character",
    default = "",
    help = "Define file name prefix for selected output files. Default: No prefix.",
    metavar = "char"
  ),
  make_option(
    c("-h", "--help"),
    action="store_true",
    default=FALSE,
    help="Show this information and die."
  )
)

## Parse options
opt_parser <- OptionParser(
  usage=paste(script, "[--help] --analysis-dir [DIR] --out-dir [DIR] 
            --comparison [CT1] --cpm-threshold [INT] --read-threshold [INT] 
            --length-threshold [DOUBLE] --iqr-threshold [INT] --ALPHA [DOUBLE] 
            --out-prefix [STRING]\n"),
  option_list=option_list,
  add_help_option=FALSE,
  description=msg,
  epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=0)

## ------------------------------------------------------------------------
library(ggplot2)
library(dplyr)


## ------------------------------------------------------------------------
# reverse CPM calculation
get_nr_reads <- function(cpm, tot_reads){
  return(round(cpm / 1e6 * tot_reads, digits = 0))
}


## ------------------------------------------------------------------------
analysis_dir <- cli$options[["analysis-dir"]]
out_dir <- cli$options[["out-dir"]]
# if mode cell_state_comparison: cell type, if mode cell_type_comparison: comparison name
ct <- cli$options[["comparison"]]
if(is.null(ct)){
  print("No comparison found. Exit script.")
  quit(save = "no", status = -1)
}

out.prefix <- cli$options[["out-prefix"]]
CPM_THRESHOLD <- cli$options[["cpm-threshold"]]
READ_THRESHOLD <- cli$options[["read-threshold"]]
LENGTH_THRESHOLD <- cli$options[["length-threshold"]]
IQR_THRESHOLD <- cli$options[["iqr-threshold"]]
PROB_TAILS <- cli$options[["ALPHA"]]
dirs <- c("coverage_3p", "randomised_coverage")
# high CPM threshold is used to separate resulting files between low and high CPM for in-depth analysis.
HIGH_CPM_THRESHOLD <- 1000

## ------------------------------------------------------------------------

dfs <- lapply(dirs, function(x) read.table(file.path(analysis_dir, x, paste0("auc_", ct, ".tsv")),
                    sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE))
df <- dfs[[1]]


## ------------------------------------------------------------------------
df <- df %>% mutate(TE_length = TE_end - TE_start)
qlengths <- quantile(df$TE_length, c(LENGTH_THRESHOLD, 0.5))

# only create histogram of TE lengths if not already exists
## ------------------------------------------------------------------------

p <- ggplot(df, aes(x=TE_length)) + 
  geom_histogram(bins=100) + 
  scale_x_log10() +
  geom_vline(xintercept = qlengths['50%'], color = "red") +
  geom_vline(xintercept = qlengths[1], color = "blue") +
  theme_bw() +
  labs(caption = paste0(nrow(df), " TEs plotted\n", 
                        "red line 50%: ", round(qlengths['50%'], digits = 1), "\n", 
                        "blue line ", names(qlengths[1]), ": ", round(qlengths[1], digits = 1)))
ggsave(file.path(out_dir, paste0(out.prefix, "hist_TE_lengths.png")),
      plot = p,
      device = "png", width = 10, height = 8, units = "cm")


## ------------------------------------------------------------------------
p <- ggplot(df, aes(x=auc)) + geom_histogram(bins = 151) + theme_bw() +
  labs(caption = paste0(length(which(!is.na(df$auc))), " TEs plotted"))
ggsave(file.path(out_dir, paste0("hist_", ct, ".png")),
       plot = p, 
       device = "png", width = 10, height = 8, units = "cm")


## ------------------------------------------------------------------------
qs1.cov <- quantile(df$coverage_sample1, probs = c(0.5,0.75))
qs2.cov <- quantile(df$coverage_sample2, probs = c(0.5,0.75))
qs1 <- quantile(df$cpm_sample1, probs = c(0.5,0.75))
qs2 <- quantile(df$cpm_sample2, probs = c(0.5,0.75))


## ------------------------------------------------------------------------
# nr reads in TE regions per sample
reads <- colSums(df[c(7,9)])

p <- ggplot(df, aes(x=coverage_sample1)) + geom_histogram(bins = 100) + 
  scale_x_log10() + 
  geom_vline(xintercept = qs1.cov['50%'], color = "red") + 
  geom_vline(xintercept = qs1.cov['75%'], color = "blue") +
  theme_bw() + 
  labs(caption = paste0(length(which(df$coverage_sample1 != 0)), " TEs plotted; ", 
                        "total: ", reads['coverage_sample1'], "\n",
                        "red line at 50% (", round(qs1.cov['50%'], digits = 1), " reads), ", 
                        "blue line at 75% (", round(qs1.cov['75%'], digits = 1), " reads)"))
ggsave(file.path(out_dir, paste0("coverage_log_", ct, "_sample1.png")),
       plot = p, 
       device = "png", width = 10, height = 8, units = "cm")
p <- ggplot(df, aes(x=coverage_sample2)) + geom_histogram(bins = 100) + 
  scale_x_log10() + 
  geom_vline(xintercept = qs2.cov['50%'], color = "red") + 
  geom_vline(xintercept = qs2.cov['75%'], color = "blue") +
  theme_bw() +
    labs(caption = paste0(length(which(df$coverage_sample2 != 0)), " TEs plotted; ", 
                        "total: ", reads['coverage_sample2'], "\n",
                        "red line at 50% (", round(qs2.cov['50%'], digits = 1), " reads), ", 
                        "blue line at 75% (", round(qs2.cov['75%'], digits = 1), " reads)"))
ggsave(file.path(out_dir, paste0("coverage_log_", ct, "_sample2.png")),
       plot = p,
       device = "png", width = 10, height = 8, units = "cm")

p <- ggplot(df, aes(x=cpm_sample1)) + geom_histogram(bins = 100) + 
  scale_x_log10() + 
  geom_vline(xintercept = qs1['50%'], color = "red") + 
  geom_vline(xintercept = qs1['75%'], color = "blue") +
  theme_bw() + 
  labs(caption = paste0(length(which(df$cpm_sample1 != 0)), " TEs plotted; 10 CPM = ", 
                        get_nr_reads(10, reads['coverage_sample1']), " reads (total: ", reads['coverage_sample1'], ")\n",
                        "red line at 50% (", round(qs1['50%'], digits = 1), " CPM), ", 
                        "blue line at 75% (", round(qs1['75%'], digits = 1), " CPM)"))
ggsave(file.path(out_dir, paste0("cpm_log_", ct, "_sample1.png")),
       plot = p,
       device = "png", width = 10, height = 8, units = "cm")
p <- ggplot(df, aes(x=cpm_sample2)) + geom_histogram(bins = 100) + 
  scale_x_log10() + 
  geom_vline(xintercept = qs2['50%'], color = "red") + 
  geom_vline(xintercept = qs2['75%'], color = "blue") +
  theme_bw() +
    labs(caption = paste0(length(which(df$cpm_sample2 != 0)), " TEs plotted; 10 CPM = ", 
                        get_nr_reads(10, reads['coverage_sample2']), " reads (total: ", reads['coverage_sample2'], ")\n",
                        "red line at 50% (", round(qs2['50%'], digits = 1), " CPM), ", 
                        "blue line at 75% (", round(qs2['75%'], digits = 1), " CPM)"))
ggsave(file.path(out_dir, paste0("cpm_log_", ct, "_sample2.png")),
       plot = p,
       device = "png", width = 10, height = 8, units = "cm")


## ------------------------------------------------------------------------
# Filter TEs by CPM_THRESHOLD and LENGTH_THRESHOLD
if(!is.null(READ_THRESHOLD)){
  df.filtered <- filter(df, coverage_sample1 >= READ_THRESHOLD, coverage_sample2 >= READ_THRESHOLD)
} else{
  df.filtered <- filter(df, cpm_sample1 >= CPM_THRESHOLD, cpm_sample2 >= CPM_THRESHOLD)
}


short_TEs <- df.filtered[df.filtered$TE_length< qlengths[1], "TE_id"]
df.filtered <- filter(df.filtered, ! TE_id %in% short_TEs)

# TODO: keep gene names filtered out and use for filter background distribution

# write out TE ids used for analysis
write.table(pull(df.filtered, TE_id),
              file = file.path(out_dir, paste0(out.prefix, "TEs_filtered_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## ------------------------------------------------------------------------
nr_TEs <- list()
if(!is.null(READ_THRESHOLD)){
  nr_TEs["sample1"] <- nrow(filter(df, coverage_sample1 >= READ_THRESHOLD))
  nr_TEs["sample2"] <- nrow(filter(df, coverage_sample2 >= READ_THRESHOLD))
} else{
  # nr. TEs in sample1 >= CPM_THRESHOLD
  nr_TEs["sample1"] <- nrow(filter(df, cpm_sample1 >= CPM_THRESHOLD))
  # nr. TEs in sample2 >= CPM_THRESHOLD
  nr_TEs["sample2"] <- nrow(filter(df, cpm_sample2 >= CPM_THRESHOLD))
}
# nr. TEs in sample1 and sample2 >= CPM_THRESHOLD
nr_TEs["both_samples"] <- nrow(df.filtered)


## ------------------------------------------------------------------------
p <- ggplot(df.filtered, aes(x=auc)) + 
  geom_histogram(aes(y=..density..), bins = 101, alpha = 0.8) +
    stat_function(fun = dnorm, 
                  args = list(mean = 0.5, 
                              sd = sd(df.filtered$auc)), 
                  n = 301, color="red") +
  theme_bw() + 
  labs(caption = paste0("N(0.5,", round(sd(df.filtered$auc), digits = 4), ")\n",
                        "mean: ", round(mean(df.filtered$auc), digits = 5), "\n",
                        nr_TEs["both_samples"], " TEs plotted\n",
                        nr_TEs["sample1"], " TEs in sample1\n",
                        nr_TEs["sample2"], " TEs in sample2"))
ggsave(file.path(out_dir, paste0(out.prefix, "hist_filtered_", ct, ".png")),
       plot = p,
       device = "png", width = 10, height = 8, units = "cm")


## ------------------------------------------------------------------------
rnorm <- rnorm(nrow(df.filtered), mean = 0.5, sd = sd(df.filtered$auc))
auc <- df.filtered$auc
ratios <- data.frame(dist = seq(0, max(0.5-min(auc), max(auc)-0.5), length.out = 76))
ratios['ratio'] <- 0
ratios['ratio_norm'] <- 0
for(i in 1:nrow(ratios)){
  geq <- length(which(auc >= 0.5 & auc <= 0.5 + ratios[i, 'dist']))
  seq <- length(which(auc <= 0.5 & auc >= 0.5 - ratios[i, 'dist']))
  r <- geq / seq
  ratios[i,'ratio'] <- ifelse(is.na(r), 0, r)
  
  geq <- length(which(rnorm >= 0.5 & rnorm <= 0.5 + ratios[i, 'dist']))
  seq <- length(which(rnorm <= 0.5 & rnorm >= 0.5 - ratios[i, 'dist']))
  r <- geq / seq
  ratios[i,'ratio_norm'] <- ifelse(is.na(r), 0, r)
}


## ------------------------------------------------------------------------
p <- ggplot(ratios, aes(x=dist, y = ratio)) + geom_line() + 
  theme_bw()
ggsave(file.path(out_dir, paste0(out.prefix, "ratio_", ct, ".png")),
       plot = p,
       device = "png", width = 8, height = 6, units = "cm")


## ------------------------------------------------------------------------
df.filtered <- mutate(df.filtered, 
                      mean_cpm = rowMeans(dplyr::select(df.filtered, starts_with("cpm_")), na.rm = TRUE))


## ------------------------------------------------------------------------
p <- ggplot(df.filtered, aes(x=mean_cpm, y = auc)) + 
  geom_point(size = 0.1, color = "black") +
  scale_fill_continuous(type = "viridis") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
  scale_x_log10() + theme_bw()
ggsave(file.path(out_dir, paste0(out.prefix, "scatter_mean_cpm_auc_", ct, ".png")),
       plot = p,
       device = "png", width = 10, height = 6, units = "cm")


## ------------------------------------------------------------------------
p <- ggplot(df.filtered, aes(x=coverage_sample1, y=coverage_sample2, color = auc)) +
  geom_point(alpha=0.8) +
  scale_color_gradient2(midpoint = 0.5, low = "red", high = "blue") +
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
  theme_bw()
ggsave(file.path(out_dir, paste0(out.prefix, "scatter_coverage_", ct, ".png")),
       plot = p,
       device = "png", width = 12, height = 8, units = "cm")

## ------------------------------------------------------------------------
if(!is.null(READ_THRESHOLD)){
  df.rnd <- filter(dfs[[2]], 
                   coverage_sample1 >= READ_THRESHOLD, coverage_sample2 >= READ_THRESHOLD)
} else{
  df.rnd <- filter(dfs[[2]], 
                   cpm_sample1 >= CPM_THRESHOLD, cpm_sample2 >= CPM_THRESHOLD)
}

df.rnd <- df.rnd %>% filter(! TE_id %in% short_TEs) %>% 
  mutate(mean_cpm = rowMeans(dplyr::select(df.filtered, starts_with("cpm_")), na.rm = TRUE)) %>% 
  mutate(mean_cpm_log = log10(mean_cpm))
# bin mean CPM
breaks <- seq(from = min(df.rnd$mean_cpm_log),
    to = (max(df.rnd$mean_cpm_log)+0.1), 
    length.out = 21)

# construct bins with quantile measurements
probs_to_use <- c(PROB_TAILS/2, (1 - PROB_TAILS/2))
quants <- data.frame()
for(i in 1:(length(breaks)-1)){
  # filter df for values in breaks
  df.bin <- filter(df.rnd, mean_cpm_log > breaks[i], mean_cpm_log < breaks[i+1])
  # obtain sample quantiles
  quants <- rbind(quants, quantile(df.bin$auc, probs = probs_to_use))
}
names(quants) <- sapply(probs_to_use, function(x) paste0("p", as.character(x*100)))
# add columns of start and end
quants$start <- 10^breaks[-length(breaks)]
quants$end <- 10^breaks[-1]

# remove bins without entries
ind_to_keep <- !is.na(quants[,1]) | !is.na(quants[,2])
quants <- quants[ind_to_keep,]
# smooth quantiles
quants[,1] <- runmed(quants[,1], k = 5, endrule = "median")
quants[,2] <- runmed(quants[,2], k = 5, endrule = "median")


## ------------------------------------------------------------------------
p <- ggplot(df.rnd, aes(x=auc)) + 
  geom_histogram(aes(y=..density..), bins = 101, alpha = 0.8) +
    stat_function(fun = dnorm, 
                  args = list(mean = 0.5, 
                              sd = sd(df.rnd$auc)), 
                  n = 301, color="red") +
  theme_bw() + 
  labs(caption = paste0("N(0.5,", round(sd(df.rnd$auc), digits = 4), ")\n",
                        "mean: ", round(mean(df.rnd$auc), digits = 5)))
ggsave(file.path(out_dir, paste0(out.prefix, "hist_background_filtered_", ct, ".png")),
       plot = p,
       device = "png", width = 10, height = 8, units = "cm")

p <- ggplot(df.rnd, aes(x=mean_cpm, y = auc)) + 
  geom_point(size = 0.1, color = "black") +
  scale_fill_continuous(type = "viridis") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
  geom_segment(data = quants, 
               aes_string(x="start", xend = "end", 
                          y=names(quants)[1], yend=names(quants)[1]), 
               color = "red") +
  geom_segment(data = quants, 
                 aes_string(x="start", xend = "end", 
                            y=names(quants)[2], yend=names(quants)[2]), 
                 color = "red") +
    scale_x_log10() + theme_bw()
ggsave(file.path(out_dir, paste0(out.prefix, "scatter_background_mean_cpm_auc_", ct, ".png")),
       plot = p,
       device = "png", width = 10, height = 6, units = "cm")



## ------------------------------------------------------------------------
p <- ggplot(df.filtered, aes(x=mean_cpm, y = auc)) + 
  geom_point(size = 0.1, color = "black") + 
  scale_fill_continuous(type = "viridis") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
  geom_segment(data = quants, 
               aes_string(x="start", xend = "end", 
                          y=names(quants)[1], yend=names(quants)[1]), 
               color = "red") +
  geom_segment(data = quants, 
                 aes_string(x="start", xend = "end", 
                            y=names(quants)[2], yend=names(quants)[2]), 
                 color = "red") +
  scale_x_log10() + theme_bw() +
  labs(x='mean CPM', y='AUC') +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 10),
        legend.position = "none")
ggsave(file.path(out_dir, paste0(out.prefix, "quantiles_", ct, ".png")),
       plot = p,
       device = "png", width = 5, height = 5, units = "cm")


## ------------------------------------------------------------------------
tes.shorter <- list()
tes.longer <- list()
for(i in 1:nrow(quants)){
  tmp <- filter(df.filtered, mean_cpm > quants[i,"start"], mean_cpm < quants[i,"end"])
  tes.shorter[[i]] <-filter(tmp, auc > quants[i, 2]) %>% pull(TE_id)
  tes.longer[[i]] <- filter(tmp, auc < quants[i, 1]) %>% pull(TE_id)
}

# save to file
write.table(unlist(tes.shorter),
              file = file.path(out_dir, paste0(out.prefix, "TEs_shorter_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unlist(tes.longer),
              file = file.path(out_dir, paste0(out.prefix, "TEs_longer_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


## ------------------------------------------------------------------------

high_cpm <- quants$start > HIGH_CPM_THRESHOLD

write.table(unlist(tes.shorter[high_cpm]),
              file = file.path(out_dir, paste0(out.prefix, "TEs_shorter_high_cpm_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unlist(tes.shorter[!high_cpm]),
              file = file.path(out_dir, paste0(out.prefix, "TEs_shorter_low_cpm_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(unlist(tes.longer[high_cpm]),
              file = file.path(out_dir, paste0(out.prefix, "TEs_longer_high_cpm_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unlist(tes.longer[!high_cpm]),
              file = file.path(out_dir, paste0(out.prefix, "TEs_longer_low_cpm_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


## ------------------------------------------------------------------------
# Find genes with significant APA events
get_span <- function(values){
  # values is a vector of length 4, containing 25% and 75% 3'UTR length for both samples
  # obtain IQR span from sample1 and sample2 (in base pairs)
  return(max(values) - min(values))
}
req.cols <- c('sample1_length_q0.25', 'sample1_length_q0.75', 
                             'sample2_length_q0.25', 'sample2_length_q0.75')
if(!is.null(IQR_THRESHOLD) & length(which(colnames(df.filtered) %in% req.cols)) == 4){
  # Filter 1PAS Terminal exons by InterQuantileRange (IQR)
  df.filtered <- mutate(df.filtered, sample1_IQR = sample1_length_q0.25 - sample1_length_q0.75,
                        sample2_IQR = sample2_length_q0.25 - sample2_length_q0.75)
  df.filtered$SPAN <- apply(dplyr::select(df.filtered, sample1_length_q0.25, sample1_length_q0.75,
                                           sample2_length_q0.25, sample2_length_q0.75),
                            1,
                             get_span)
  
  df.iqr <- filter(df.filtered, SPAN >= IQR_THRESHOLD)
  
  # Write out TE ids with IQR >= IQR_THRESHOLD and sign. AUC change
  tes.shorter.iqr <- filter(df.iqr, TE_id %in% unlist(tes.shorter)) %>% pull(TE_id)
  write.table(tes.shorter.iqr,
              file = file.path(out_dir, paste0(out.prefix, "TEs_shorter_IQR_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  tes.longer.iqr <- filter(df.iqr, TE_id %in% unlist(tes.longer)) %>% pull(TE_id)
  write.table(tes.longer.iqr,
              file = file.path(out_dir, paste0(out.prefix, "TEs_longer_IQR_", ct, ".tsv")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  ## ------------------------------------------------------------------------
  # Add info about significance in df.filtered
  sign.auc <- c(unlist(tes.shorter), unlist(tes.longer))
  df.filtered$is_sign_auc <- sapply(df.filtered$TE_id, function(x) ifelse(x %in% sign.auc, TRUE, FALSE))
  sign.iqr <- c(tes.shorter.iqr, tes.longer.iqr)
  df.filtered$is_sign_iqr <- sapply(df.filtered$TE_id, function(x) ifelse(x %in% sign.iqr, TRUE, FALSE))
  # Plot quantiles
  cols <- c("TRUE" = "darkgreen", "FALSE" = "grey")
  size <- c("TRUE" = 0.5, "FALSE" = 0.1)
  p <- ggplot(df.filtered, aes(x=mean_cpm, y = auc, color = is_sign_iqr)) + 
    geom_point(aes(size = is_sign_iqr)) + 
    scale_color_manual(values = cols) +
    scale_size_manual(values = size) + 
    geom_segment(data = quants, 
                 aes_string(x="start", xend = "end", 
                            y=names(quants)[1], yend=names(quants)[1]), 
                 color = "red") +
    geom_segment(data = quants, 
                 aes_string(x="start", xend = "end", 
                            y=names(quants)[2], yend=names(quants)[2]), 
                 color = "red") +
    scale_x_log10() + theme_bw() +
    labs(x='mean CPM', y='AUC') +
    theme(axis.text = element_text(size=8),
          axis.title = element_text(size = 10),
          legend.position = "none")
  ggsave(file.path(out_dir, paste0(out.prefix, "quantiles_iqr_", ct, ".png")),
         plot = p,
         device = "png", width = 5, height = 5, units = "cm")
} else if(!is.null(IQR_THRESHOLD) & req.cols != 4){
  print("The desired columns for IQR test are not provided.")
  print(paste("Columns required: ", paste(req.cols, collapse = " ")))
  print(paste("Columns provided: ", paste(colnames(df.filtered), collapse = " ")))
}

# save complete filtered df
write.table(df.filtered,
            file = file.path(out_dir, paste0(out.prefix, "TEs_filtered_complete_", ct, ".tsv")), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#######

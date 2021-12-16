#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { 
  stop("Package 'optparse' required!\nExecution aborted.") 
}

## ------------------------------------------------------------------------
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

## ------------------------------------------------------------------------
description <- "Plot pathway enrichment for multiple pathways as heatmap. \n"
author <- "Author: Dominik Burri"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: dominik.burri@unibas.ch"
version <- "Version: 1.0.0 (27-APR-2021)"
requirements <- c("optparse", "dplyr", "ggplot2", "pheatmap", "RColorBrewer")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version, requirements_txt, sep="\n")
notes <- "Notes: No additional notes."

## ------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--sample-names",
    action = "store",
    type = "character",
    default = NULL,
    help = ("Comma separated list of sample names. List entries are matched with in-files. Default: sample[ind]")
  ),
  make_option(
    "--min-size",
    action = "store",
    type = "integer",
    default = 10,
    help = "Define minimum pathway size. Default: %default",
    metavar = "int"
  ),
  make_option(
    "--max-size",
    action = "store",
    type = "integer",
    default = 2000,
    help = "Define maximum pathway size. Default: %default",
    metavar = "int"
  ),
  make_option(
    "--top-n",
    action = "store",
    type = "integer",
    default = 10,
    help = "Top n pathways to plot. Default: %default",
    metavar = "int"
  ),
  make_option(
    "--out-file",
    action = "store",
    type = "character",
    default = "top_pathways.png",
    help = "Define output file. Default: top_pathways.png"
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
  usage=paste(script, "[--help] [--sample-names] [--min-size] [--max-size] [--top-n] [--out-file] IN-FILES \n"),
  option_list=option_list,
  add_help_option=FALSE,
  description=msg,
  epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=c(2, Inf))

## ------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

## ------------------------------------------------------------------------
in_files <- strsplit(cli$args, ",")
if(is.null(cli$options[["sample-names"]])){
  sample_names <- paste0("sample", 1:length(in_files))
} else{
  sample_names <- strsplit(cli$options[["sample-names"]], ",")[[1]]
}
out_file <- cli$options[["out-file"]]
min_size <- cli$options[["min-size"]]
max_size <- cli$options[["max-size"]]
n <- cli$options[["top-n"]]


## ------------------------------------------------------------------------
dfs <- lapply(in_files, function(x) read.table(x, sep = "\t", stringsAsFactors = FALSE))
column_names <- function(df){
  colnames(df) <- c("term", "description", "observed", 
                    "background", "strength", "FDR", 
                    "matching_protein_ids", "matching_protein_labels")
  return(df)
}
dfs <- lapply(dfs, column_names)

# Filter pathways by size
filter_size <- function(df, min_size, max_size){
  df <- filter(df, background >= min_size, background <= max_size)
  return(df)
}
dfs <- lapply(dfs, filter_size, min_size = min_size, max_size = max_size)

# Combine all dataframes into one
df.all <- select(dfs[[1]], term, description, FDR)
colnames(df.all)[3] <- sample_names[1]
for(i in 2:length(dfs)){
  tmp <- select(dfs[[i]], term, description, FDR)
  colnames(tmp)[3] <- sample_names[i]
  df.all <- full_join(df.all, tmp, by = c('term', 'description'))
}

# Obtain union of top n pathways
get_top_n <- function(df, n){
  terms <- arrange(df, FDR) %>% slice_head(n = n) %>% select(term, description)
  return(terms)
}
top_ns <- lapply(dfs, get_top_n, n = n)


union_terms <- Reduce(function(...) full_join(..., by = c('term', 'description')), top_ns)

# filter dfs for union terms
filter_union_terms <- function(df, terms){
  df <- filter(df, term %in% terms)
  return(df)
}
dfs.filtered <- lapply(dfs, filter_union_terms, terms = union_terms$term)

# log10 transform FDR
dfs.filtered <- lapply(dfs.filtered, function(x) mutate(x, logFDR = log10(FDR)))

# Combine all dataframes into one
df <- union_terms
for(i in 1:length(dfs)){
  tmp <- select(dfs.filtered[[i]], term, logFDR)
  colnames(tmp)[2] <- sample_names[i]
  df <- left_join(df, tmp, by = 'term')
}

# restrict description to m chars
m <- 50
clip_description <- function(string, m){
  if(nchar(string) > m){
    string <- paste0(substr(string, 1, m), "*")
  } 
  return(string)
}
df$description <- sapply(df$description, clip_description, m=m)

# make description unique
ind <- duplicated(df$description)
df[ind, "description"] <- paste0(df[ind, "description"], "_1")

# remove term from columns and add description as row name
rownames(df) <- df$description
df <- select(df, -term, -description)

# only plot intersection of top n
# df <- na.omit(df)

# Plot heatmap
pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(brewer.pal(8, "Reds")),
         na_col = "Grey",
         fontsize = 8,
         legend =TRUE,
         width = 10 / cm(1), height = 12 / cm(1),
         filename = out_file)


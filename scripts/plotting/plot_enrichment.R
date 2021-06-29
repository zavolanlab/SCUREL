#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { stop("Package 'optparse' required!\nExecution aborted.") }

## ------------------------------------------------------------------------
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

## ------------------------------------------------------------------------
description <- "Plot pathway enrichment. \n"
author <- "Author: Dominik Burri"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: dominik.burri@unibas.ch"
version <- "Version: 1.0.0 (27-APR-2021)"
requirements <- c("optparse", "dplyr", "ggplot2")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version, requirements_txt, sep="\n")
notes <- "IN-FILE: file with significant enrichment pathway results, column specification as for string-db export."

## ------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--min-size",
    action = "store",
    type = "integer",
    default = 0,
    help = "Define minimum pathway size. Default: %default"
  ),
  make_option(
    "--max-size",
    action = "store",
    type = "integer",
    default = 1e5,
    help = "Define maximum pathway size. Default: %default"
  ),
  make_option(
    "--top-n",
    action = "store",
    type = "integer",
    default = 10,
    help = "Top n pathways to plot. Default: %default"
  ),
  make_option(
    "--out-file",
    action = "store",
    type = "character",
    default = "top_pathways.png",
    help = "Define output file. Default: top_pathways.png",
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
  usage=paste(script, "[--help] [--min-size] [--max-size] [--top-n] [--out-file] IN-FILE \n"),
  option_list=option_list,
  add_help_option=FALSE,
  description=msg,
  epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=1)

## ------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

## ------------------------------------------------------------------------
in_file <- cli$args[1]
out_file <- cli$options[["out-file"]]
min_size <- cli$options[["min-size"]]
max_size <- cli$options[["max-size"]]
n <- cli$options[["top-n"]]

## ------------------------------------------------------------------------
df <- read.table(in_file, sep = "\t", stringsAsFactors = FALSE)
colnames(df) <- c("term", "description", "observed", 
                  "background", "strength", "FDR", 
                  "matching_protein_ids", "matching_protein_labels") 

# restrict description to m chars
m <- 40
clip_description <- function(string, m){
  if(nchar(string) > m){
    string <- paste0(substr(string, 1, m), "*")
  } 
  return(string)
}
df$description <- sapply(df$description, clip_description, m=m)

# Filter pathways by size
df <- filter(df, background >= min_size, background <= max_size)

# Filter pathways by n
df <- arrange(df, FDR) %>% slice_head(n = n)

# log10 transform FDR
df <- mutate(df, logFDR = log10(FDR))

# make description unique
ind <- duplicated(df$description)
df[ind, "description"] <- paste0(df[ind, "description"], "_1")

# Plot bargraph
p <- ggplot(df, aes(x=reorder(description, -logFDR), y=logFDR)) + 
  geom_col() + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "log10 FDR") +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
ggsave(out_file, plot = p,
       width = 5, height = 3.5, units = "cm",
       scale = 2,
       dpi = 500)



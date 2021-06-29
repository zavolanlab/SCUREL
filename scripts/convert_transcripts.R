#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { stop("Package 'optparse' required!\nExecution aborted.") }

## ------------------------------------------------------------------------
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

## ------------------------------------------------------------------------
description <- "Convert Refseq identifiers to gene names. \n"
author <- "Author: Dominik Burri"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: dominik.burri@unibas.ch"
version <- "Version: 1.0.0 (10-MAY-2021)"
requirements <- c("optparse", "biomaRt")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version, requirements_txt, sep="\n")
notes <- "IN-FILE: file with transcript identifiers to be converted."

## ------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--organism",
    action = "store",
    type = "character",
    default = "hsapiens_gene_ensembl",
    help = "Organism specification, see listDatasets(useMart('ensembl')) for details. Default: hsapiens_gene_ensembl.",
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
  usage=paste(script, "[--help] [--organism ORG] IN-FILE OUT-FILE  \n"),
  option_list=option_list,
  add_help_option=FALSE,
  description=msg,
  epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=2)

## ------------------------------------------------------------------------
library(biomaRt)

## ------------------------------------------------------------------------
in_file <- cli$args[1]
out_file <- cli$args[2]
organism <- cli$options[["organism"]]

out <- read.table(in_file, sep = "\t", stringsAsFactors = FALSE)

# get gene names
new_config <- httr::config(ssl_verifypeer = FALSE)
httr::set_config(new_config, override = FALSE)

conv_ttog <- function(myvector, mytable){
  res <- sapply(myvector, function(x) mytable[mytable[,1] == x, 2][1])
  return(unname(res))
}

# select proper organism!
mart <- useDataset(organism, useMart("ensembl"))
transcripts <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),
                     values = out, mart= mart)

# remove duplicates
gene_names <- unique(transcripts$external_gene_name)

write.table(gene_names, 
            file = out_file,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

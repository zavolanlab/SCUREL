#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { stop("Package 'optparse' required!\nExecution aborted.") }

## ------------------------------------------------------------------------
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

## ------------------------------------------------------------------------
description <- "Get background genome from expressed transcripts. \n"
author <- "Author: Dominik Burri"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: dominik.burri@unibas.ch"
version <- "Version: 1.0.0 (07-MAY-2021)"
requirements <- c("optparse")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version, requirements_txt, sep="\n")
notes <- "IN-FILES: files with gene names or transcript identifiers to be considered."

## ------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--out-file",
    action = "store",
    type = "character",
    default = "identifiers.tsv",
    help = "Define output file. Default: identifiers.tsv",
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
  usage=paste(script, "[--help] [--organism ORG] -o OUT-FILE IN-FILES  \n"),
  option_list=option_list,
  add_help_option=FALSE,
  description=msg,
  epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=c(2, Inf))

## ------------------------------------------------------------------------
in_files <- cli$args
out_file <- cli$options[["out-file"]]
organism <- cli$options[["organism"]]

in_vec <- vector()
for(inf in in_files){
  f <- read.table(inf, sep = "\t", stringsAsFactors = FALSE)
  in_vec <- append(in_vec, f[,1])
}

out <- unique(in_vec)

write.table(out, 
            file = out_file,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

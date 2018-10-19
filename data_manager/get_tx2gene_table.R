library(getopt)

# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "input_type", "t", 1, "character",
  "outfile", "o", 1, "character",
  "gtfFile", "x", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

suppressPackageStartupMessages({library("GenomicFeatures")})
txdb <- makeTxDbFromGFF(opt$gtfFile, format=opt$input_type)
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1] # tx ID, then gene ID
write.table(tx2gene,file = opt$outfile, quote = FALSE, sep = " ",row.names = FALSE,col.names = FALSE)


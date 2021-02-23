# Find technical replicates with different cluster assignment
# and map csSNP performance

library("RMySQL") # mysql connection
library("dplyr") # set operations

installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- paste(installdir, "/curated_SNPs/", sep = "")

# Current known csSNP's
snp_df <-
  read.csv(paste(snp_dir, "02-07-2020_csSNPs.csv", sep = ""), sep = ' ')

# wgsid data for cluster ids
con_wgsid <- dbConnect(RMySQL::MySQL(), group = "wgsid")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'wgsid')
dbDisconnect(con_wgsid)

# Connect to annot3 database for Dutch isolate SNPs
con_annot3 <- dbConnect(RMySQL::MySQL(), group = "annot3")
annot.data <- dbReadTable(conn = con_annot3, name = 'annot3')
dbDisconnect(con_annot3)
# Limit scope to SNPs for now, no indels
snp.data <- droplevels(subset(annot.data, VAR1 == "SNP"))
snp.data <- snp.data[, c('strain', 'run', 'VAR5', 'VAR6', 'proxclad')]
# Filter out unique strain+run combinations
pairs <- unique(snp.data[,c('strain','run')])
pairs <- data.frame(pairs)
# map wgsids to pairs
ij_pairs <-
  inner_join(
    pairs,
    wgsid.data,
    by = c("strain", "run")
  )
# calculate isolates with >1 runs (technical replicates)
dups <- data.frame(table(ij_pairs$strain))
dups <- dups[dups$Freq>1,]
dups_list <- dups$Var1
# leave only technical replicates
ij_pairs <- ij_pairs[ij_pairs$strain %in% dups_list,]
ij_pairs$wgsidnum <- NULL
# find replicates with mismatching wgsid's
for (dup in dups_list) {
  wgsids <- unique(ij_pairs[ij_pairs$strain==dup,"wgsid"])
  size <- ij_pairs[ij_pairs$strain==dup,"Freq"]
  # if more than one distinct wgsids for a replicate
  if (length(wgsids) >1) {
    print(dup)
    print(wgsids)
    print(size)
  }
}

# Mapping resistance-related SNP's to csSNP collection

library("RMySQL") # mysql connection
library("ggplot2") # plotting

installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- paste(installdir, "/curated_SNPs/", sep = "")

# Current known csSNP's
snp_df <-
  read.csv(paste(snp_dir, "02-07-2020_csSNPs.csv", sep = ""), sep = ' ')

# Connect to annot3 database for Dutch isolate SNPs
con_annot3 <- dbConnect(RMySQL::MySQL(), group = "annot3")
annot.data <- dbReadTable(conn = con_annot3, name = 'annot3')
dbDisconnect(con_annot3)
# Limit scope to SNPs for now, no indels
snp.data <- droplevels(subset(annot.data, VAR1 == "SNP"))
# Remove non-resistance columns
resi_snps <-
  snp.data[, c("VAR5", "VAR6", "pza")]
# Filter out csSNP's not in that region
snp_df_pnca <- 
  # e.g. pncA region: 2288681-2289240
  snp_df[snp_df$VAR5 < 2289241 & snp_df$VAR5 > 2288680, ]
# Result df
resi_df <- NULL
# Per csSNP, check if it overlaps with resistance.
for (i in 1:nrow(snp_df_pnca)) {
  print(paste(i, nrow(snp_df_pnca), sep = " / "))
  wgsid <- snp_df_pnca[i, "WGSID"]
  posi <- snp_df_pnca[i, "VAR5"]
  allele <- snp_df_pnca[i, "VAR6"]
  matches <-
    resi_snps[resi_snps$VAR5 == posi & resi_snps$VAR6 == allele, ]
  if (nrow(matches) > 0) {
    str(matches)
    resi_df <-
      rbind(resi_df, data.frame(wgsid, posi, allele, matches))
  }
}
# remove redundant columns
resi_df$posi <- NULL
resi_df$allele <- NULL
resi_df$clade <- NULL
resi_df[is.na(resi_df)] <- ""
# frequencies
table_pnca <- data.frame(table(resi_df[, c("VAR6", "pza")]))#,useNA="always"))
range_nrs <- seq(2288681,2289240,50)
if (tail(range_nrs,n=1) < 2289240) {
  range_nrs <- append(range_nrs, 2289240)
}
resi_df$VAR5 <- as.integer(resi_df$VAR5)

# plot histogram of allele resistance distribution
p <- ggplot(table_pnca, aes(VAR6, Freq, fill = pza)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_minimal() + 
  labs(title = "Distribution of intra-pncA-gene csSNP's", 
       x = 'csSNP allele', y = 'Frequencies')
p

# Save resistance df to file
write.table(
  resi_df,
  paste(snp_dir, "pncA_csSNPs.txt", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  col.names = c(
    "wgsid",
    "VAR5",
    "VAR6",
    "pza resistance"
  ),
  sep = '\t'
)

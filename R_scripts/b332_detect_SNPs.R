# Set working directory
installdir <- "/data/haddadt/NL_SNP"
setwd(dir=installdir)

library("RMySQL") # myswl connection
library("dplyr") # set operations

###################################
# B332 cluster retrieval          #
###################################

# Find small cluster to detect characteristic SNP
# mysql: select * from wgsid where wgsid = "B332";
con_wgsid <- dbConnect(RMySQL::MySQL(), group = "wgsid")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'wgsid')
names(wgsid.data)
dim(wgsid.data)
# Close connection
dbDisconnect(con_wgsid)

# Find subset of B332 WGS ID, n=17
wgsid.data <- droplevels(subset(wgsid.data, wgsid=="B332"))
# Write to csv file, ; as sep
write.csv(wgsid.data, "B332_cluster.csv", row.names = FALSE)
# List of B332 strains
b332_strains <- c(wgsid.data$strain)

###################################
# annot3 isolates retrieval       #
###################################

# Connect to annot3 database for Dutch isolate SNPs
con_annot3 <- dbConnect(RMySQL::MySQL(), group = "annot3")
annot.data <- dbReadTable(conn = con_annot3, name = 'annot3')
names(annot.data)
dim(annot.data)
dbDisconnect(con_annot3)

# Limit scope to SNPs for now, no indels
snp.data <- droplevels(subset(annot.data, VAR1=="SNP"))
write.csv(snp.data, "annot3_SNP_only.csv", row.names = FALSE)

###################################
# Separate sets on strains        #
###################################

# Filter out B332 isolates from annot3 set
b332_subset <- droplevels(subset(snp.data, snp.data$strain %in% b332_strains))
write.csv(b332_subset, "B332_annot3_subset.csv", row.names = FALSE)
# annot3 subset without B332 cluster
annot3_subset <- droplevels(subset(snp.data, !(snp.data$strain %in% b332_strains)))
write.csv(annot3_subset, "annot3_without_B332.csv", row.names = FALSE)

###################################
# Find B332 exclusive muts        #
###################################

# B332 cluster mutations
b332_muts <- data.frame(AL_POS = c(b332_subset$VAR5), 
                        AL_NEW = c(b332_subset$VAR6))
# All annot3 mutations
annot3_muts <- data.frame(AL_POS = c(annot3_subset$VAR5), 
                          AL_NEW = c(annot3_subset$VAR6))
# B332 muts not in rest of annot3 isolates
b332_specific <- setdiff(b332_muts, annot3_muts)

# Calculate frequency of b332 specific muts
# dplyr to count frequency of variable combinations
b332_sub_spec <- droplevels(subset(b332_subset, 
                                  (b332_subset$VAR5 %in% b332_specific$AL_POS) & 
                                    (b332_subset$VAR6 %in% b332_specific$AL_NEW)))
b332_sub_pre_cur <- count(b332_sub_spec, vars = c("VAR5","VAR6"))
# final table, all cluster-unique SNPs identified in every single isolate
b332_sub_cur <- b332_sub_pre_cur[b332_sub_pre_cur$freq==17,]

# Filter out non-synonymous SNPs
b332_silent_snps <- droplevels(subset(b332_subset, 
                                 (b332_subset$VAR5 %in% b332_sub_cur$VAR5) & 
                                   (b332_subset$VAR6 %in% b332_sub_cur$VAR6) &
                                   (b332_subset$snp_type == "synonymous")))
# display location + SNP freq
table(b332_silent_snps$VAR5, b332_silent_snps$VAR6)


# Here I show that B332 has 13 unique strains, but 17 strain+run combs total
#strain_table <- table(c(b332_subset$strain)) # Strains in b332 cluster
#b332_000001_strains <- subset(b332_subset, mut=="000001") # Strains with mut 000001, n=13
#View(table(b332_000001_strains$strain)) # Table, shows 17 total, some strains >1


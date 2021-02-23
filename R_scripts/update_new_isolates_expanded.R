# Updater script for csSNP finding in newly added isolates

library("RMySQL") # mysql connection

installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- paste(installdir, "/curated_SNPs/", sep = "")

# Current known csSNP's
snp_df <-
  read.csv(paste(snp_dir, "02-07-2020_csSNPs.csv", sep = ""), sep = ' ')

###################################
# WGS clusters retrieval          #
###################################
print("Collecting wgs clusters")
# Collecting wgs id's of new isolates, if they're already assigned
con_wgsid <- dbConnect(RMySQL::MySQL(), group = "wgsid")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'wgsid')
dbDisconnect(con_wgsid)

###################################
# annot3 isolates retrieval       #
###################################
print("Collecting annot3 data")
# Connect to annot3 database for Dutch isolate SNPs
con_annot3 <- dbConnect(RMySQL::MySQL(), group = "annot3")
annot.data <- dbReadTable(conn = con_annot3, name = 'annot3')
dbDisconnect(con_annot3)
# Limit scope to SNPs for now, no indels
snp.data <- droplevels(subset(annot.data, VAR1 == "SNP"))

###################################
# strain distance retrieval
print("Collecting strain snp distances")
con_dist <- dbConnect(RMySQL::MySQL(), group = "distance3")
dist.data <- dbReadTable(conn = con_dist, name = 'distance3')
dbDisconnect(con_dist)

###################################
# EDIT THIS DATE for future runs  #
###################################
# Find latest entry date (run id) to target isolates added after that
last_date <- 200702 # 2 juli 2020
snp.data$run <- as.numeric(snp.data$run)
snp.data <- droplevels(subset(snp.data, run > last_date))
# Give standard format to both dataframes
snp.data$VAR5 <- as.character(snp.data$VAR5)
snp.data$VAR6 <- as.character(snp.data$VAR6)
snp_df$VAR5 <- as.character(snp_df$VAR5)
snp_df$VAR6 <- as.character(snp_df$VAR6)
# Filter out unnecessary data
snp.data <- snp.data[, c('strain', 'run', 'VAR5', 'VAR6', 'proxclad')]

# Total new snps found
new_snps <- nrow(snp.data)
# Nr of unique strains added
new_strains <- unique(snp.data$strain)
new_strains <- unique(snp.data[, c('strain', 'run')])
snp.data$found <- NULL

# Collecting statistics
nwgs <- 0 # no wgsid
ncssnp <- 0 # no csSNP found
ma <- 0 # actual wgsid+csSNP match
mm <- 0 # mismatch wgsid / csSNP
totalm <- 0 # total nr of strains that have a csSNP found
result_df <- NULL # df with end results of isolates+csSNP's

# Do a snp check on a per strain basis
#for (i in 1:nrow(new_strains)) {
for (i in 1:5) {
  iso <- new_strains[i, 'strain']
  run <- new_strains[i, 'run']
  id <- paste(iso, run, sep = '_')
  print(id)
  # Check if a wgsid is already assigned
  assigned <-
    unique(wgsid.data[wgsid.data$strain == iso &
                        wgsid.data$run == run, c('wgsid', 'Freq')])
  if (nrow(assigned) > 0) {
    print(paste(assigned, sep = ' '))
    #assigned <- as.list(assigned$wgsid,assigned$Freq)
    assigned <- as.list(assigned[1, ])
  } else {
    print('No wgsid found')
    nwgs <- nwgs + 1
    assigned <- list("None", "0")
  }
  # All snps for the specific strain
  strain_snps <- snp.data[snp.data$strain == iso & snp.data$run == run, ]
  # Iterate known csSNP's
  csSNP_wgsids <- unique(snp_df$WGSID)
  # Is the assigned WGS cluster present in the csSNP collection?
  has_csSNP <- FALSE
  if (assigned %in% csSNP_wgsids) {
    has_csSNP <- TRUE
  } else {
    has_csSNP <- FALSE
  }
  found_wgs <- list()
  mab <- FALSE  # Boolean for match
  for (wgsid in csSNP_wgsids) {
    # csSNP's for the wgsid calculated previously
    wgsid_snps <- snp_df[snp_df$WGSID == wgsid, ]
    # do the isolate snp's and wgs csSNP's have entries in common?
    snp_overlap <-
      nrow(semi_join(strain_snps, wgsid_snps, by = c('VAR5', 'VAR6')))
    # Filter out likely FP's with 1 cluster size
    if (snp_overlap > 0) {
      # add length of cluster size of matched csSNP
      wgsid_size <-
        unique(wgsid.data[wgsid.data$wgsid == wgsid, c('wgsid', 'Freq')])[1, 2]
      if (wgsid_size > 1) {
        found_wgs <- append(found_wgs, c(wgsid, wgsid_size, snp_overlap))
        print(wgsid)
        # if assigned wgs equals csSNP wgs
        if (wgsid == assigned) {
          ma <- ma + 1
          mab <- TRUE
        } else {
          mm <- mm + 1
        }
        print(mab)
      }
      
      
    }
  }
  # replace empty list with None
  if (length(found_wgs) == 0) {
    found_wgs <- "None"
  }
  # remove NULL elements from list
  found_wgs <- found_wgs[lengths(found_wgs) != 0]
  print(found_wgs)
  # If no match found, no csSNP found, else +1 to total
  if (mab == FALSE) {
    ncssnp <- ncssnp + 1
  } else {
    totalm <- totalm + 1
  }
  # Append result df
  result_df <- rbind(
    result_df,
    # iso_run id, isolate, run id, wgsid, wgsid size, csSNP's found + their size, match(TRUE)/mismatch(FALSE)
    data.frame(
      id,
      iso,
      run,
      as.character(assigned[1]),
      as.integer(assigned[2]),
      has_csSNP,
      toString(found_wgs),
      as.character(mab)
    )
  )
}
colnames(result_df) <-
  c(
    "ID",
    "isolate",
    "run",
    "assigned cluster",
    "assigned cluster size",
    "assigned WGS has csSNP",
    "found csSNP's + its size",
    "csSNP matches wgsid"
  )
date <- format(Sys.time(), "%d-%m-%Y")
write.table(
  result_df,
  paste(snp_dir, date, "_updated_isolates.txt_no_singles", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  col.names = c(
    "ID",
    "isolate",
    "run",
    "assigned cluster",
    "assigned cluster size",
    "assigned WGS has csSNP",
    "found csSNP's + its size",
    "csSNP matches wgsid"
  ),
  sep = '\t'
)
print(paste("distinct new isolates: ", length(unique(new_strains$strain)), sep = ''))
print(paste("total nr of new isolate-run pairs: ", nrow(new_strains), sep = ''))
print(paste("no wgsid assigned: ", nwgs, sep = ''))
print(paste("no csSNP found: ", ncssnp, sep = ''))
print(paste("total nr of csSNP's found: ", totalm, sep = ''))
print(paste("wgsid/csSNP matches: ", ma, sep = ''))

###################################
# Filter and benchmark            #

###################################
# Remove (x,1,1) (likely FP's)
filtered_results <- result_df
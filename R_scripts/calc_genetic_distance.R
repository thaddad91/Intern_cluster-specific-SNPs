###################################
# Calculate genetic distance      #
###################################
# Author: Thierry Haddad

# Calculate the genetic distance, based on WGS sequencing,
# between isolates concatenated during the negative set solver.

# c1   size c2
# A037	5	 B601
# B491	4	 B772
# A022	2	 A444
# A528	2	 A656
# A799	2	 B667
# B771	2	 B770
# A507	1	 A428

# Set working directory
print("Setting up dirs")
installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- "/data/haddadt/NL_SNP/curated_SNPs/"

library("RMySQL") # mysql connection

###################################
# WGS clusters retrieval
print("Collecting wgs clusters")
con_wgsid <- dbConnect(RMySQL::MySQL(), group = "wgsid")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'wgsid')
dbDisconnect(con_wgsid)

###################################
# strain distance retrieval
print("Collecting strain snp distances")
con_dist <- dbConnect(RMySQL::MySQL(), group = "distance3")
dist.data <- dbReadTable(conn = con_dist, name = 'distance3')
dbDisconnect(con_dist)

# convert concatenation set into a df
neg_set <- data.frame(
  c("A037", "5", "B601"),
  c("B491", "4", "B772"),
  c("A022", "2", "A444"),
  c("A528", "2", "A656"),
  c("A799", "2", "B667"),
  c("B771", "2", "B770"),
  c("A507", "1", "A428")
)
neg_set <- t(neg_set)

# result df
distances <- NULL

# Iterate over clusters to find average distance
for (i in 1:nrow(neg_set)) {
  # collect strain+run for cluster 1
  row <- neg_set[i,]
  c1 <- row[1]
  size <- as.numeric(row[2])
  c1_strains <-
    wgsid.data[wgsid.data$wgsid == c1[1], c("strain", "run")]
  c1_strains <- paste(c1_strains$strain, c1_strains$run, sep = "_")
  # collect strain+run for cluster 2
  c2 <- row[3]
  c2_strains <-
    wgsid.data[wgsid.data$wgsid == c2[1], c("strain", "run")]
  c2_strains <- paste(c2_strains$strain, c2_strains$run, sep = "_")
  # calculate distances using distance3 db
  sum_dist <- 0
  print(c1)
  for (str_run in c1_strains) {
    distance <-
      dist.data[dist.data$isolate1 == str_run &
                  dist.data$isolate2 == c2_strains,]
    # if no distance metrics where found
    if (nrow(distance) == 0) {
      print(paste("no dist: ", str_run, c2_strains, sep = " "))
      # remove 1 from c1 size for avg calc
      size <- size - 1
    } else {
      # add found distance to sum
      distance <- as.numeric(distance$distance)
      print(paste(str_run,c2_strains,distance,sep = " "))
      sum_dist <- sum_dist + distance
    }
  }
  # calc avg distance from isolates c1 <-> c2
  avg <- sum_dist / size
  if (avg == c1) {
    avg <- "NA"
  }
  distances <- rbind(distances, c(c1, size, c2, avg))
}
colnames(distances) <- c("c1", "size", "c2", "avg dist")

write.table(
  distances,
  paste(snp_dir, "neg_dist.csv", sep = ""),
  sep = ",",
  row.names = FALSE
)

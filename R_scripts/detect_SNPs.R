###################################
# Summary of pipeline             #
###################################
# Author: Thierry Haddad
#
# The function of this pipeline is to identify and quantify
# cluster-specific SNP's for Dutch M. tuberculosis isolates.
# - First, the annotated mutation database Annot3 and WGS cluster
#   database wgsid are pulled and merged on allele and loci.
# - Second, per WGS clusters the presence of csSNP's are being
#   calculated (calc_snps())
# - Third, for those clusters that did not find useable csSNP's,
#   an iterative approach was applied by adding single-isolate clusters
#   to the respective cluster, and retrying the csSNP identification.
# - Lastly, the results are saved in tables/csv formats and some plots are made.

# Set working directory
print("Setting up dirs")
installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- "/data/haddadt/NL_SNP/curated_SNPs/"

library("RMySQL") # mysql connection
library("dplyr") # set operations
library("ggplot2") # plotting

###################################
# WGS clusters retrieval          #
###################################
print("Collecting wgs clusters")
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
# Iterate over clusters through   #
# generalization of functions     #
###################################
print("Extracting wgs clusters from annot3")
# Fuse together annot3 + wgsid
test_annot <- snp.data
test_wgs <- wgsid.data

# inner-join
ij_annot <-
  inner_join(
    test_annot,
    test_wgs,
    by = c("strain", "run"),
    suffix = c(".an", ".wg")
  )

# Scatter plot of WGS cluster size against size frequency
tmp_wgsid <- unique(as.data.frame(wgsid.data[, c("wgsid", "Freq")]))
colnames(tmp_wgsid)[2] <- "size"
tmp_wgsid2 <- as.data.frame(table(tmp_wgsid$size))
# Had to do type conversion, as table() merged two columns wrong
tmp_wgsid2$Var1 <- as.numeric(as.character(tmp_wgsid2$Var1))
# Save plot of WGS cluster size against its frequency
jpeg(file = "/data/haddadt/NL_SNP/curated_SNPs/wgs_size_x_size_frequency.jpeg",
     width = 600,
     height = 300)
ggplot(tmp_wgsid2, aes(x = Var1, y = Freq)) +
  geom_point(
    size = 2,
    shape = 20,
    color = "blue",
    fill = "blue"
  ) +
  labs(title = "Frequency of WGS cluster sizes",
       x = "WGS cluster size", y = "Frequency of size") + xlim(0, 54)
dev.off()

# Gather all cluster names as a vector
clusters <- unique(ij_annot$wgsid)
nr_of_clusters <- length(clusters)
# DF to capture cluster-specific SNPs
snp_df <- NULL
# Vector for clusters with >1 Coll lineage
multi_lin <- list()
# Vector for clusters with 1 isolate, used for solving non-csSNP clusters
single_iso <- NULL
# Vector for clusters w/o result
no_snps <- NULL

calc_snps <- function() {
  print("Calculating csSNPs")
  # DF to capture cluster-specific SNPs
  snp_df <<- NULL
  # Vector for clusters with >1 Coll lineage
  multi_lin <<- list()
  # Vector for clusters with 1 isolate, used for solving non-csSNP clusters
  single_iso <<- NULL
  # Vector for clusters w/o result
  no_snps <<- NULL
  # Iterate over clusters
  c <- 0
  clusters <- c("B752")
  for (cluster in clusters) {
    # Boolean if cluster was solved
    found <- 0
    c <- c + 1
    # Get size of cluster
    cluster_size <-
      nrow(dplyr::count_(ij_annot[ij_annot$wgsid == cluster, ], vars = c("strain", "run")))
    print(paste(cluster, "-", cluster_size, "-", c, "/", nr_of_clusters, sep = ""))
    # Temporary data sets divided by cluster
    tmp_cluster <- droplevels(subset(ij_annot, wgsid == cluster))
    tmp_ref <- droplevels(subset(ij_annot, wgsid != cluster))
    # Identify Coll's clades in cluster
    proxclads <- unique(tmp_cluster$proxclad)
    if (length(proxclads) > 1) {
      # Catch how many clusters >1 lineage
      multi_lin <<- c(multi_lin, cluster)
    }
    proxclad <- as.character(paste(proxclads, collapse = "+"))
    # Find cluster-specific SNPs based on SNP location + allele
    cluster_specific <-
      anti_join(tmp_cluster,
                tmp_ref,
                by = c("VAR5", "VAR6"),
                keep = TRUE)
    if (nrow(cluster_specific) > 0) {
      print(cluster_specific)    
      # Count frequencies of cluster-specific SNPs (csSNP's)
      cluster_curated <-
        cluster_specific[cluster_specific$VAR1 == "SNP", ]
      cluster_curated <-
        dplyr::count_(cluster_curated, vars = c("VAR5", "VAR6"))
      print(cluster_curated)
      # Gather those that equal cluster size (omni-present)
      cluster_curated <-
        cluster_curated[cluster_curated$n == cluster_size, ]
      # Remove redundant frequency
      cluster_curated$n <- NULL
      # If any csSNP at all
      if (nrow(cluster_curated) > 0) {
        for (row in 1:nrow(cluster_curated)) {
          # Append csSNP's to dataframe
          c_wgsid <- cluster
          c_var5 <- cluster_curated[row, "VAR5"]
          c_var6 <- cluster_curated[row, "VAR6"]
          c_proxclad <- proxclad
          snp_df <<-
            rbind(snp_df,
                  data.frame(c_wgsid, c_var5, c_var6, c_proxclad))
        }
        found <- 1
      } else {
        no_snps <<- rbind(no_snps,
                          data.frame(cluster, cluster_size))
      }
    } else {
      no_snps <<- rbind(no_snps,
                        data.frame(cluster, cluster_size))
    }
    # Catch single iso clusters + csSNP found or not
    if (cluster_size == 1) {
      single_iso <<- rbind(single_iso,
                           data.frame(cluster, found))
    }
  }
  # Temp save to test neg function
  print("Saving DF to CSV...")
  date <- format(Sys.time(), "%d-%m-%Y")
  write.table(
    snp_df,
    paste(snp_dir, date, "_snp_df.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    single_iso,
    paste(snp_dir, date, "_single_iso.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    no_snps,
    paste(snp_dir, date, "_no_snps.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE
  )
}
calc_snps()

#####################################
# Expanding clusters without csSNPs #
#####################################
# Here the set of clusters without csSNPs will receive
# a single isolate and be retested for csSNPs, where negative
# isolates have priority over single isolates with csSNPs

# Read the csv's if DF's are empty
date <- format(Sys.time(), "%d-%m-%Y")
if (is.null(snp_df)) {
  print("No snp DF found, reading from CSV...")
  csv_name <-
    as.character(paste(snp_dir, date, "_snp_df.csv", sep = ""))
  snp_df <- read.table(csv_name, header = TRUE, sep = " ")
}
if (is.null(single_iso)) {
  print("No single iso DF found, reading from CSV...")
  csv_name <-
    as.character(paste(snp_dir, date, "_single_iso.csv", sep = ""))
  single_iso <- read.table(csv_name, header = TRUE, sep = " ")
}
if (is.null(no_snps)) {
  print("No non_snps DF found, reading from CSV...")
  csv_name <-
    as.character(paste(snp_dir, date, "_no_snps.csv", sep = ""))
  no_snps <- read.table(csv_name, header = TRUE, sep = " ")
}
# Collect some meta-data for downstream comparison
pre_length <- nrow(snp_df)
pre_iso <- nrow(single_iso)
pre_nosnp <- nrow(no_snps)

# Collect successful concatenations
used <- list()

# Calculate negative cluster specific SNPs (not in snp_df)
neg_csSNPs <- function() {
  print("Solving negative set")
  # Sort negative clusters on size, high->low
  no_snps <<-
    no_snps[order(as.numeric(no_snps$cluster_size), decreasing = TRUE),]
  # Sort single isolates to those without csSNP are used first
  single_iso <<-
    single_iso[order(as.numeric(single_iso$found), decreasing = FALSE), ]
  # Iterate over all clusters
  for (j in 1:nrow(no_snps)) {
    cluster <- no_snps[j, ]$cluster
    print(paste("CLUSTER", cluster, sep = " ---- "))
    # Temporary data sets divided by cluster
    tmp_cluster <- droplevels(subset(ij_annot, wgsid == cluster))
    # Identify Coll's clades in cluster
    proxclads <- unique(tmp_cluster$proxclad)
    proxclad <- as.character(paste(proxclads, collapse = "+"))
    
    # Iterate single isolates for merger
    for (i in 1:nrow(single_iso)) {
      iso_cluster <- single_iso[i,]$cluster
      print(paste(cluster, "testing isolate:", iso_cluster, sep = " "))
      # Reference set, must not contain current cluster wgsid + current isolate wgsid
      tmp_ref <-
        droplevels(subset(ij_annot, wgsid != cluster &
                            wgsid != iso_cluster))
      # Merge the isolate with the cluster
      iso_data <-
        droplevels(subset(ij_annot, wgsid == iso_cluster))
      tmp_cluster2 <- rbind(tmp_cluster, iso_data)
      # Find cluster-specific SNPs based on SNP 
      tmp_ref$VAR5 <- as.character(tmp_ref$VAR5)
      tmp_cluster2$VAR5 <- as.character(tmp_cluster2$VAR5)
      cluster_specific <-
        anti_join(tmp_cluster2,
                  tmp_ref,
                  by = c("VAR5", "VAR6"),
                  keep = TRUE)
      # If any csSNP at all
      if (nrow(cluster_specific) > 0) {
        # Count frequencies of cluster-specific SNPs (csSNP's)
        cluster_curated <-
          cluster_specific[cluster_specific$VAR1 == "SNP", ]
        cluster_curated <-
          dplyr::count_(cluster_curated, vars = c("VAR5", "VAR6"))
        # Get size of new cluster
        cluster_size <-
          nrow(dplyr::count_(ij_annot[ij_annot$wgsid == cluster, ],
                             vars = c("strain", "run")))
        # Account for isolate
        cluster_size <- cluster_size + 1
        # Gather those that equal cluster size (omni-present)
        cluster_curated <-
          cluster_curated[cluster_curated$n == cluster_size, ]
        cluster_curated$n <- NULL
        # If omni-present csSNPs, means cluster+iso merge was successful
        if (nrow(cluster_curated) > 0) {
          print(paste('Solved', cluster, 'with', iso_cluster, sep = ' '))
          for (row in 1:nrow(cluster_curated)) {
            # Append csSNP's to dataframe
            c_wgsid <- cluster
            c_var5 <- cluster_curated[row, "VAR5"]
            c_var6 <- cluster_curated[row, "VAR6"]
            c_proxclad <- proxclad
            snp_df <<-
              rbind(snp_df,
                    data.frame(c_wgsid, c_var5, c_var6, c_proxclad))
          }
          # Remove cluster and single iso from iteration, to
          # avoid being used again afterwards
          single_iso <<-
            single_iso[!single_iso$cluster == iso_cluster, ]
          no_snps <<- no_snps[!no_snps$cluster == cluster, ]
          used[[as.character(cluster)]] <<-
            as.character(iso_cluster)
          # Cluster now has csSNP's, so exit this loop
          break
        } else {
          print(paste(iso_cluster, "failed...", sep = " "))
        }
      }
      #
    }
  }
}
neg_csSNPs()

# Comparison of pre- and post-negative cluster set solving
post_length <- nrow(snp_df)
post_iso <- nrow(single_iso)
post_nosnp <- nrow(no_snps)
print("snp df:")
print(pre_length)
print(post_length)
print("single isolates total:")
print(pre_iso)
print(post_iso)
print("no_snps:")
print(pre_nosnp)
print(post_nosnp)
print("used cluster-isolate pairs:")
print(used)
print('no csSNPs clusters:')
print(no_snps)

write_results <- function() {
  # Write results to files
  date <- format(Sys.time(), "%d-%m-%Y")
  write.table(
    snp_df,
    paste(snp_dir, date, "_csSNPs.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    col.names = c("WGSID", "VAR5", "VAR6", "proxclad")
  )
  write.csv(
    no_snps,
    paste(snp_dir, date, "_no_result_clusters.csv", sep = ""),
    row.names = TRUE,
    quote = FALSE
  )
  write.csv(
    multi_lin,
    paste(snp_dir, date, "_multi_lineages.csv", sep = ""),
    row.names = TRUE,
    quote = FALSE
  )
  
  # nr of SNPs per wgsid
  cluster_freq <- as.data.frame(table(snp_df$c_wgsid))
  # frequency of cluster-specific SNP sizes
  size_freq <- as.data.frame(table(cluster_freq$Freq))
  
  write.table(
    cluster_freq,
    paste(snp_dir, "cluster_frequencies.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    col.names = c("WGSID", "Frequency")
  )
  write.table(
    size_freq,
    paste(snp_dir, "frequency_of_cluster_sizes.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    col.names = c("Unique_SNP_size", "Frequency")
  )
  
  ###################################
  # Plot frequencies of cluster-    #
  # specific SNPs                   #
  ###################################
  
  # let's plot
  x_nr <- as.numeric(size_freq$Var1)
  y_nr <- as.numeric(size_freq$Freq)
  jpeg(file = "/data/haddadt/NL_SNP/curated_SNPs/size_freq.jpeg")
  plot(x_nr,
       y_nr,
       xlab = "Number of cluster-specific SNPs",
       ylab = "Frequency",
       main = "Frequencies of cluster-specific SNP group sizes")
  # Add polynomial curve for visual clarity
  lo <- loess(y_nr ~ x_nr)
  lines(predict(lo), col = 'red', lwd = 2)
  dev.off()
  
  ###################################
  # Comparison of cluster size      #
  # relative to SNP group size      #
  ###################################
  
  #Cluster sizes
  c_sizes <- as.data.frame(table(wgsid.data$wgsid))
  colnames(c_sizes)[2] <- "cluster_size"
  c_sizes <- inner_join(c_sizes, cluster_freq, by = c("Var1"))
  colnames(c_sizes)[1] <- "WGSID"
  colnames(c_sizes)[3] <- "specific_SNPs"
  
  # let's plot
  cx_nr <- as.numeric(c_sizes$cluster_size)
  cy_nr <- as.numeric(c_sizes$specific_SNPs)
  jpeg(file = "/data/haddadt/NL_SNP/curated_SNPs/csize_x_snpsize.jpeg")
  plot(cx_nr,
       cy_nr,
       xlab = "WGSID cluster size",
       ylab = "Nr of cluster-specific SNPs",
       main = "Comparison of WGSID cluster size to number of cluster-specific SNPs")
  dev.off()
}
write_results()
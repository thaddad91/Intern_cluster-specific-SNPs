# Quick script for fixing the csSNP plots

installdir <- "/data/haddadt/NL_SNP"
setwd(dir = installdir)
snp_dir <- paste(installdir, "/curated_SNPs/", sep = "")

library("ggplot2") # plotting
library("numbers") # fibonacci nrs for bins

# Read frequencies from files
cluster_freq <-
  read.csv(paste(snp_dir, "cluster_frequencies.csv", sep = ""), sep = ' ')
cluster_freq <-
  cluster_freq[!cluster_freq$Frequency == 0, ] # Drop zero-levels
size_freq <-
  read.csv(paste(snp_dir, "frequency_of_cluster_sizes.csv", sep = ""),
           sep = ' ')
size_freq <-
  size_freq[!size_freq$Unique_SNP_size == 0, ] # Drop zero-levels

###################################
# Comparison of cluster size      #
# relative to SNP group size      #
###################################

con_wgsid <- dbConnect(RMySQL::MySQL(), group = "wgsid")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'wgsid')
dbDisconnect(con_wgsid)

#Cluster sizes
c_sizes <- as.data.frame(table(wgsid.data$wgsid))
colnames(c_sizes)[2] <- "cluster_size"
colnames(c_sizes)[1] <- "WGSID"
c_sizes <- inner_join(c_sizes, cluster_freq, by = c("WGSID"))
colnames(c_sizes)[3] <- "specific_SNPs"

# let's plot
cx_nr <- log(as.numeric(c_sizes$cluster_size), 2)
bin_cx_nr <- as.data.frame(table(cut(cx_nr, c(0, 1, 2, 3, 4, 5, 6))))
cy_nr <- log(as.numeric(c_sizes$specific_SNPs), 2)
plot(cx_nr,
     cy_nr,
     xlab = "Log2 cluster size",
     ylab = "Log2 number of cluster-specific SNPs",
     main = "Comparison of WGSID cluster size to number of cluster-specific SNPs")
xc_labels <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7")
freq2c <- ggplot(data = bin_cx_nr, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Comparison of cluster size to number of csSNPs",
       x = 'Log2 cluster size', y = 'Log2 number of cluster-specific SNPs') +
  scale_x_discrete(labels = xc_labels)
freq2c

###################################
# Plot frequencies of cluster-    #
# specific SNPs                   #
###################################

# Frequencies of nr of csSNP found
y_nr <- as.numeric(cluster_freq$Frequency)

# Dividing the ranges into bins for clarity, fibonacci ranges
range_nrs <- unique(fibonacci(16, sequence = TRUE))
bin_counts <- as.data.frame(table(cut(y_nr, range_nrs)))
bin_counts$Freq <- log(bin_counts$Freq, 2)
xlabels <-
  list(
    '1-2',
    '2-3',
    '3-5',
    '5-8',
    '8-13',
    '13-21',
    '21-34',
    '34-55',
    '55-89',
    '89-144',
    '144-233',
    '233-377',
    '377-610',
    '610-987'
  )
freq2 <- ggplot(data = bin_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Frequencies of the number of csSNP's found per cluster",
       x = 'binned number of csSNPs found per cluster', y = 'Log2 frequencies') +
  scale_x_discrete(labels = xlabels)
freq2
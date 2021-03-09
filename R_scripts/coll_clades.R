# Updater script for csSNP finding in newly added isolates

library("RMySQL") # mysql connection
library("ggplot2") # plotting

installdir <- "masked_dir"
setwd(dir = installdir)
snp_dir <- paste(installdir, "/curated_SNPs/", sep = "")

# Current known csSNP's
snp_df <-
  read.csv(paste(snp_dir, "02-07-2020_csSNPs.csv", sep = ""), sep = ' ')

# WGS cluster size frequencies
print("Collecting wgs clusters")
# Collecting wgs id's of new isolates, if they're already assigned
# masked several parameters
con_wgsid <- dbConnect("masked")
wgsid.data <- dbReadTable(conn = con_wgsid, name = 'masked')
dbDisconnect(con_wgsid)
wgsid.freq <- unique(as.data.frame(wgsid.data[, c("wgsid", "Freq")]))
range <- ceiling(max(wgsid.freq$Freq)/10)*10
range_nrs <- c(0,1,5,10,15,20,25,30,35,40,45,50,55,60)
xlabels <- c("1","2-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50-55","55-60")
bin_counts <- as.data.frame(table(cut(wgsid.freq$Freq, range_nrs)))
bin_counts$Freq <- log(bin_counts$Freq,2)
freq2 <- ggplot(data = bin_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Frequencies of WGS cluster sizes",
       x = 'Binned WGS cluster sizes', y = 'Frequencies') +
  scale_x_discrete(labels = xlabels)
freq2

###################################
# annot3 isolates retrieval       #
###################################
print("Collecting annot3 data")
# Connect to annot3 database for Dutch isolate SNPs
con_annot3 <- dbConnect("masked")
annot.data <- dbReadTable(conn = con_annot3, name = 'masked')
dbDisconnect(con_annot3)
# Limit scope to SNPs for now, no indels
snp.data <- droplevels(subset(annot.data, VAR1 == "SNP"))

###################################
# EDIT THIS DATE for future runs  #
###################################
# Find latest entry date (run id) to target isolates added after that
last_date <- 200702 # 2 juli 2020
snp.data$run <- as.numeric(snp.data$run)
snp.data <- droplevels(subset(snp.data, run <= last_date))


snp_df2 <- unique(snp_df[,c("WGSID","proxclad")])
clade_freq <- as.data.frame(table(snp_df2$proxclad))
clade_freq$Freq <- log(clade_freq$Freq,2)
freq <- ggplot(data = clade_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Frequencies of Coll clades in clusters with csSNP's",
       x = 'Clade categories', y = 'Frequencies')
freq

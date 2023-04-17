SNP_start <- readRDS("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/start/SNP_POS.rds")
SNP_end <- readRDS("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/end/SNP_POS.rds")

SNP_POS <- c(SNP_start,SNP_end)

#STEP ONE:
sep_MNV_POS <- strsplit(SNP_POS, ":")

#STEP TWO:
mnv_chr <- lapply(sep_MNV_POS, `[[`, 1)
mnv_strtend <- lapply(sep_MNV_POS, `[[`, 2)

library(dplyr)
my_data <- data.table(chr = mnv_chr)

# Use the group_by and summarize functions from dplyr to count the frequency of each unique element
obj <-my_data %>% 
  group_by(chr) %>%
  summarize(rec_count = n()) 
obj <- as.data.frame(obj)
obj$rec_count <- obj$rec_count*1000


snp_chr.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/SNP_chr.rds")
saveRDS(obj, snp_chr.rds)



````

MNV_start <- readRDS("dbSNP-151_data/start/MNV_POS.rds")
MNV_end <- readRDS("dbSNP-151_data/end/MNV_POS.rds")
SNP_POS <- c(MNV_start,MNV_end)

#STEP ONE:
sep_MNV_POS <- strsplit(SNP_POS, ":")

#STEP TWO:
mnv_chr <- lapply(sep_MNV_POS, `[[`, 1)
mnv_strtend <- lapply(sep_MNV_POS, `[[`, 2)

library(dplyr)
my_data <- data.table(chr = mnv_chr)

# Use the group_by and summarize functions from dplyr to count the frequency of each unique element
obj <-my_data %>% 
  group_by(chr) %>%
  summarize(rec_count = n()) 
obj <- as.data.frame(obj)
obj$rec_count <- obj$rec_count*1000



bi_snps_chr <- readRDS("dbSNP-144_data/bi_snps_chr.rds")

non_bi_snps_chr <- readRDS("dbSNP-151_data/non_bi_snps_chr.rds")

non_bi_snps_chr <- data.table::as.data.table(non_bi_snps_chr)

path <- ("GCF_000001405.25.gz.first100krows.vcf")

VariantAnnotation::readVcf(path)



a <- readLines(path, n=2000)

vcfDF <- VariantAnnotation::readVcf(path)

vcfDF <- vcf2df(vcfDF)

bi_snps_MA <- readRDS("dbSNP-150_data/bi_snps_MA.rds")


obj <- vcfDF %>% 
  group_by(FREQ) %>%
  summarize(rec_count = n()) 
obj <- as.data.frame(obj)


res <- prop.test( x= c(non_bi_snp_no_144,non_bi_snp_no_150,non_bi_snps_151,non_bi_snps_155), n= c(snp_no_144,snp_no_150,snps_151,snps_155))

snps_155 <- data.table(names=c("Bi-SNPs","Non-Bi SNPs"), rec_count=c(bi_snps_no_155,non_bi_snps_no_155)) 

#Establishing non-bi & snp numbers from builds:
non_bi_snp_no_144 <- readRDS("dbSNP-144_data/non_bi_snps_NO.rds")
non_bi_snp_no_150 <- readRDS("dbSNP-150_data/non_bi_snps_NO.rds")
non_bi_snps_151 <- readRDS("dbSNP-151_data/start/non_bi_snps_NO.rds") + readRDS("dbSNP-151_data/end/non_bi_snps_NO.rds")
non_bi_snps_155 <- readRDS("dbSNP-155_data/start/non_bi_snps_NO.rds") + readRDS("dbSNP-155_data/mid/non_bi_snps_NO.rds") + readRDS("dbSNP-155_data/end/non_bi_snps_NO.rds")
snp_no_144 <- readRDS("dbSNP-144_data/SNP_NO.rds")
snp_no_150 <- readRDS("dbSNP-150_data/SNP_NO.rds")
snps_151 <- readRDS("dbSNP-151_data/start/SNP_NO.rds") + readRDS("dbSNP-151_data/end/SNP_NO.rds")
snps_155 <- readRDS("dbSNP-155_data/start/SNP_NO.rds") + readRDS("dbSNP-155_data/mid/SNP_NO.rds") + readRDS("dbSNP-155_data/end/SNP_NO.rds")
save.image("~/Documents/GitHub/MungeSumStats/MungeSumStats/dbSNP-155.RData")
snps_155 <- data.table(names=c("Bi-SNPs","Non-Bi SNPs"), rec_count=c(bi_snp_no_144,non_bi_snp_no_144)) 

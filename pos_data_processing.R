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







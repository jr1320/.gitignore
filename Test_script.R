#Creating vcf2df fucntion if MSS integrated function causes error:
#drop_duplicate_rows is commented out as this is causing forderv error.

#Messager function:
messager <- function(..., v = TRUE) {
  if (v) {
    msg <- paste(...)
    message(msg)
  }
}
#vcf2df function:
vcf2df <- function(vcf, 
                   add_sample_names=TRUE,
                   add_rowranges=TRUE,
                   drop_empty_cols=TRUE,
                   unique_cols=TRUE,
                   unique_rows=TRUE,
                   unlist_cols=TRUE,
                   sampled_rows=NULL,
                   verbose=TRUE) {
  requireNamespace("VariantAnnotation")
  requireNamespace("MatrixGenerics")
  
  messager("Converting VCF to data.table.",v=verbose) 
  #### Expand VCF ####
  if (methods::is(vcf,"CollapsedVCF")) {
    messager('Expanding VCF first, so number of rows may increase.',
             v=verbose)
    ## Not all VCFs work with this function
    vcf <- tryCatch({
      VariantAnnotation::expand(x = vcf)
    }, error = function(e){message(e); vcf}) 
    print("vcf expanded. Line 30")
  }  
  #### Automatically determine whether to add sample names ####
  if(is.null(add_sample_names)){
    header <- VariantAnnotation::header(x = vcf)
    samples <- VariantAnnotation::samples(header)
    add_sample_names <- length(samples)!=1
  } 
  #### .anncols function ####
  .anncols = function(anncol,headerstring) {
    anncols <- strsplit(sub("Functional annotations: '",'',
                            headerstring),' \\| ')[[1]]
    dfannempty <- data.frame(matrix(vector(), 0, length(anncols),
                                    dimnames=list(c(), anncols)),
                             stringsAsFactors=FALSE)
    yy <- data.frame(
      suppressWarnings(
        do.call(
          rbind,
          c(dfannempty,lapply(lapply(anncol,`[`,1),
                              function(x){strsplit(x,'\\|')[[1]]})
          )
        )
      ),
      stringsAsFactors=FALSE)
    yy <- data.frame(lapply(yy,type.convert))
    colnames(yy) <- paste("ANN",anncols,sep="_")
    return(yy)
  }
  print("Line 62")
  #### convert to data.table #### 
  t1 <- Sys.time()
  #### Get rowranges only if available ####
  if("rowRanges" %in% methods::slotNames(vcf)){
    gr <- MatrixGenerics::rowRanges(vcf)
  } else {
    gr <- NULL
  }
  df <- data.table::data.table(
    ID = names(gr),
    if(add_rowranges) granges_to_dt(gr = gr) else NULL,
    DF_to_dt(DF = VariantAnnotation::info(vcf)) 
  )
  remove(gr)
  #### Parse ANN column ####
  if('ANN' %in% colnames(df)) {
    dfann <- .anncols(
      anncol = df$ANN,
      headerstring = VariantAnnotation::info(
        VariantAnnotation::header(vcf)
      )['ANN',]$Description
    )
    df <- df[,colnames(df)!="ANN"]
    df <- cbind(df,dfann)
  }
  ##### Convert geno data ##### 
  ## Works better for VCFHeader
  if(methods::is(vcf,"VCFHeader")){
    tmp <- data.table::as.data.table(VariantAnnotation::geno(vcf), 
                                     keep.rownames = "name")
    # tmp <- DF_to_dt(DF = VariantAnnotation::geno(x))    
  } else {
    ## Works better for VCF
    n <- names(VariantAnnotation::geno(vcf))
    #### Avoid parsing redundant columns ####
    if("ID" %in% colnames(df)) n <- n[n!="ID"]
    tmp <- lapply(n,function(col) {
      coldat <- VariantAnnotation::geno(vcf)[[col]]  
      #### Convert 3D matrix --> 2D ####
      ## This sometimes happens by accident with 
      ## VariantAnnotation::writeVcf
      if(length(dim(coldat))==3){
        rn <- rownames(coldat)
        cn <- rep(colnames(coldat), dim(coldat)[3])
        if(length(unique(colnames(coldat)))==1){
          coldat <- matrix(data = coldat[,,1, drop=FALSE])
          rownames(coldat) <- rn
          colnames(coldat) <- cn[1]
        } else {
          coldat <- matrix(data = coldat,
                           nrow = dim(coldat)[1], 
                           ncol = dim(coldat)[2] * dim(coldat)[3])
          rownames(coldat) <- rn
          colnames(coldat) <- cn
        }  
      }
      print("Line 120")
      #### Drop empty columns within each matrix ####
      if(isTRUE(drop_empty_cols)){
        for(i in ncol(coldat)){
          if(all(is.na(coldat[,i])) || 
             all(coldat[,i]==".") ||
             all(coldat[,i]=="")){
            coldat <- coldat[,-i]
          }
        }
      } 
      #### Check if empty ####
      if(ncol(coldat)==0 | nrow(coldat)==0){
        return(NULL)
      } else {
        ## keeps colnames unchanged
        data.table::as.data.table(coldat)
      } 
    }) ## end lapply
    
    #### Remove NULL entries ####
    nulls <- unname(unlist(lapply(tmp,is.null)))
    tmp <- tmp[!nulls]
    n <- n[!nulls]
    ## Each element can potentially have >1 column 
    ncols <- unlist(lapply(tmp,ncol))
    tmp <- do.call(cbind, tmp)
    #### Add column names ####
    if(!is.null(tmp)){
      if(isTRUE(add_sample_names)){
        colnames(tmp) = paste(rep(n, times = ncols), 
                              colnames(tmp),sep = "_") 
      } else {
        colnames(tmp) <- rep(n, times = ncols)
      } 
    }
  } 
  df <- cbind(df, tmp) 
  remove(tmp)
  
  ## ----- Post-processing ----- ####
  #### Only keep unique rows ####
  #### Remove duplicated columns ####
  if(isTRUE(unique_cols)){
    drop_duplicate_cols(dt = df)
  } 
  #### Remove any remaining empty columns #####
  if(isTRUE(drop_empty_cols)){
    remove_empty_cols(sumstats_dt = df,
                      sampled_rows = sampled_rows,
                      verbose = verbose)
  }
  #### Unlist columns inplace ####
  if(isTRUE(unique_rows) && isFALSE(unlist_cols)){
    messager("Must set unlist_cols=TRUE to use unique_rows=TRUE.",
             "Setting unlist_cols=TRUE.",v=verbose)
    unlist_cols <- TRUE
  }
  if(isTRUE(unlist_cols)){ 
    unlist_dt(dt = df,
              verbose = verbose)
  } 
  #### Remove duplicated rows #### 
  if(isTRUE(unique_rows)){
    # df <- drop_duplicate_rows(dt = df, 
    # verbose = verbose)
  }  
  #### Report ####
  if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
  messager("VCF data.table contains:",
           formatC(nrow(df),big.mark = ","),"rows x",
           formatC(ncol(df),big.mark = ","),"columns.",
           v=verbose)
  return(df)
}

#granges_to_dt
granges_to_dt  <- function(gr) {
  if(is.null(gr)) return(gr)
  #### Convert metadata ####
  DF <- GenomicRanges::elementMetadata(gr)
  #### Combine metadata with ranges data ####
  if( ncol(DF) > 0) { 
    meta <- DF_to_dt(DF = DF)
    DT <- data.table::data.table(
      chr=as.vector(GenomicRanges::seqnames(gr)), 
      start=GenomicRanges::start(gr),
      end=GenomicRanges::end(gr), 
      meta)
  } else {
    DT <- data.table::data.table(
      chr=as.vector(GenomicRanges::seqnames(gr)),
      start=GenomicRanges::start(gr),
      end=GenomicRanges::end(gr))
  }
  return(DT)
}

#DF_to_dt:
DF_to_dt <- function(DF){ 
  #### Check if DF is empty ####
  # DF <- cbind(DF,dummy=NA)
  if(nrow(DF)==0 | ncol(DF)==0) return(data.table::data.table())
  m <- mapply(DF, 
              FUN=function(s){ 
                # s <- DF[["REF"]]
                # s <- DF[["ALT"]]
                # s <- DF[[1]]
                if(methods::is(s,"DNAStringSet") ){
                  s <- as.character(s)
                } else if(methods::is(s,"DNAStringSetList")){
                  s <- IRanges::CharacterList(s)
                  s <- Biostrings::unstrsplit(s, sep=",")
                } else if(methods::is(s,"NumericList")){
                  s <- as.numeric(s)
                } else if(methods::is(s,"list")){
                  s <- unlist(s)
                } else {
                  s <- as.vector(s)
                }
                #### Check if empty ####
                if(all(is.na(s)) || 
                   all(s==".") || 
                   all(s=="")){
                  return(NULL)
                } else {
                  return(s)
                } 
              }) 
  data.table::as.data.table(m)
}

#drop_duplicate_cols:
drop_duplicate_cols <- function(dt){
  dups <- which(duplicated(names(dt)))
  if(length(dups)>0){
    messager("Dropping",length(dups),"duplicate column(s).")
    dt[,  which(duplicated(names(dt))):= NULL] 
  } 
}

#remove_empty_cols:
remove_empty_cols <- function(sumstats_dt, 
                              sampled_rows=NULL,
                              verbose=TRUE){
  #### Check for empty columns ####
  messager("Checking for empty columns.",v=verbose)
  empty_cols <- check_empty_cols(
    sumstats_dt = sumstats_dt,
    sampled_rows = sampled_rows,
    verbose = FALSE
  )
  #### Remove empty columns #####
  if(length(empty_cols)>0) {
    messager("Removing",length(empty_cols),"empty columns.",v=verbose)
    sumstats_dt[,(names(empty_cols)):=NULL] 
  }
} 

#check_empty_cols:
check_empty_cols <- function(sumstats_dt,
                             sampled_rows = NULL, 
                             verbose = TRUE) {
  if (is.null(sampled_rows)) {
    sampled_rows <- nrow(sumstats_dt)
  } else {
    sampled_rows <- min(sampled_rows, nrow(sumstats_dt))
  }
  empty_cols <- vapply(
    colnames(sumstats_dt), function(x) {
      dt <- utils::head(sumstats_dt, sampled_rows)
      (sum(unlist(dt[[x]]) != ".") == 0) |
        (sum(!is.na(unlist(dt[[x]]))) == 0) |
        (sum(unlist(dt[[x]]) != 0) == 0)
    },
    FUN.VALUE = logical(1)
  )
  #sometimes NA values appear, they aren't empty so change these to false
  empty_cols <- ifelse(is.na(empty_cols),FALSE,empty_cols)
  empty_cols <- empty_cols[empty_cols]
  messager(length(empty_cols), "empty column(s) detected.",
           v=verbose)
  return(empty_cols)
}

#unlist_dt:
unlist_dt <- function(dt,
                      verbose=TRUE) {
  .SD <- NULL
  cols <- names(dt)[ unlist(lapply(dt, methods::is,"list")) ]
  if(length(cols)>0){
    messager("Unlisting",length(cols),"columns.",v=verbose)
    dt[ , (cols) := lapply(.SD,unlist), .SDcols = cols]
  } 
}

#drop_duplicate_rows:
drop_duplicate_rows <- function(dt,
                                verbose=TRUE){
  nrows_old <- nrow(dt)
  dt <- unique(dt)
  nrows_new <- nrow(dt)
  dropped <-  nrows_old-nrows_new
  if(dropped>0){
    messager("Dropped",formatC(dropped,big.mark = ","),
             "duplicate rows.",v=verbose)
  }
  return(dt)
}


#Read-in VCF:
path <-(path)
vcf<- VariantAnnotation:: readVcf(path,"hg19")


# default start and end:
increment <- 1200000
strt1 <- 1
end1 <-  0 + increment

# Calculate the number of iterations
iterations <- ceiling(nrow(vcf) / increment)


#Setting objects to 0
SNP_NO <- 0
INDEL_NO <- 0
MNV_NO <- 0
non_bi_snps_NO <- 0
bi_snps_NO <- 0

bi_snps_POS <- NULL
INDEL_POS<- NULL
non_bi_snps_POS <- NULL
SNP_POS <- NULL
MNV_POS <- NULL

INDEL_MA <- 0
SNP_MA <- 0
MNV_MA <- 0
non_bi_snps_MA <- 0
bi_snps_MA <- 0



# Loop through the data.frame in chunks

while (strt1<nrow(vcf)){
  
  chunk <- vcf[strt1:end1,]
  
  chunk <- vcf2df(chunk,unlist_cols = FALSE)

  ### ANALYSIS
  #Create new column for position with chr, start & end:
  chunk[,pos_chr:= paste0(chr,":",start,"-",end)]
  
  #Identifying INDELs:
  INDELS <- (chunk[VC=="DIV"])
  INDEL_NO <- INDEL_NO + length(unique(INDELS$ID))
  INDEL_POS <- c(INDEL_POS, unique(INDELS$pos_chr))
  #Count number that are G5A and/or G5 =TRUE for minor allele frequency 
  #INDEL_G5 <- INDEL_G5 + nrow(INDELS[G5==TRUE])
  #INDEL_G5A <-  (INDEL_G5A + nrow(INDELS[G5A==TRUE]))
  #INDEL_MA <- c(INDEL_G5,INDEL_G5A)
  remove(INDELS)
  
  #Identifying SNPs:
  SNPs <- (chunk[VC=="SNV"])
  SNP_NO <- SNP_NO + length(unique(SNPs$ID))
  SNP_POS <- c(SNP_POS, unique(SNPs$pos_chr))
  #SNP_G5 <- SNP_G5 + nrow(SNPs[G5==TRUE])
  #SNP_G5A <- SNP_G5A + nrow(SNPs[G5A==TRUE])
  #SNP_MA <- c(SNP_G5,SNP_G5A)
  remove(SNPs)
  
  #Identifying MNV:
  MNVs <- (chunk[VC=="MNV"])
  MNV_NO <- MNV_NO + length(unique(MNVs$ID))
  MNV_POS <- c(MNV_POS, unique(MNVs$pos_chr)) 
  #MNV_G5 <- MNV_G5 + nrow(MNVs[G5==TRUE])
  #MNV_G5A <- MNV_G5A + nrow(MNVs[G5A==TRUE])
  #MNV_MA <- c(MNV_G5,MNV_G5A)
  remove(MNVs)
  
  #Identifying non-biallelic SNPs:
  non_bi_snps <- chunk[nchar(ALT)>1 & VC=="SNV",]
  non_bi_snps_NO <- non_bi_snps_NO + length(unique(non_bi_snps$ID))
  non_bi_snps_POS <-  c(non_bi_snps_POS, unique(non_bi_snps$pos_chr)) 
  #non_bi_snps_G5 <- non_bi_snps_G5 + nrow(non_bi_snps[G5==TRUE])
  #non_bi_snps_G5A <- non_bi_snps_G5A + nrow(non_bi_snps[G5A==TRUE])
  #non_bi_snps_MA <- c(non_bi_snps_G5,non_bi_snps_G5A)
  remove(non_bi_snps)
  
  #Identifying biallelic SNPs:
  bi_snps <- chunk[nchar(ALT)==1 & VC=="SNV",]
  bi_snps_NO <- bi_snps_NO + length(unique(bi_snps$ID))
  bi_snps_POS <-  c(bi_snps_POS, unique(bi_snps$pos_chr)) 
  #bi_snps_G5 <- bi_snps_G5 + nrow(bi_snps[G5==TRUE])
  #bi_snps_G5A <- bi_snps_G5A + nrow(bi_snps[G5A==TRUE])
  #bi_snps_MA <- c(bi_snps_G5,bi_snps_G5A)
  remove(bi_snps)
  
  
  # Update the start for the next iteration
  strt1 <- strt1 + increment
  end1 <- end1 + increment
  
  if(end1>nrow(vcf)){
    end1 <- nrow(vcf)
  }
  
}

#Remove objects:
remove(chunk)


#Print results:
print(SNP_NO)
print(INDEL_NO)
print(MNV_NO)
print(non_bi_snps_NO)
print(bi_snps_NO)


### EXPORTING DATA:


SNP_NO.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/SNP_NO.rds")
saveRDS(SNP_NO, "SNP_NO.rds")
SNP_POS.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/SNP_POS.rds")
saveRDS(SNP_POS, file = "SNP_POS.rds")
SNP_MA.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/SNP_MA.rds")
saveRDS(SNP_MA, file = "SNP_MA.rds")


INDEL_NO.rds <-("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/INDEL_NO.rds")
saveRDS(INDEL_NO, "INDEL_NO.rds")
INDEL_POS.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/INDEL_POS.rds")
saveRDS(INDEL_POS, "INDEL_POS.rds")
INDEL_MA.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/INDEL_MA.rds")
saveRDS(INDEL_MA, "INDEL_MA.rds")


MNV_NO.rds <-("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/MNV_NO.rds")
saveRDS(MNV_NO, "MNV_NO.rds")
MNV_POS.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/MNV_POS.rds")
saveRDS(MNV_POS, "MNV_POS.rds")
MNV_MA.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/MNV_MA.rds")
saveRDS(MNV_MA, "MNV_MA.rds")


non_bi_snps_NO.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/non_bi_snps_NO.rds")
saveRDS(non_bi_snps_NO, "non_bi_snps_NO.rds")
non_bi_snps_POS.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/non_bi_snps_POS.rds")
saveRDS(non_bi_snps_POS, "non_bi_snps_POS.rds")
non_bi_snps_MA.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/non_bi_snps_MA.rds")
saveRDS(non_bi_snps_MA, "non_bi_snps_MA.rds")


bi_snps_NO.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/bi_snps_NO.rds")
saveRDS(bi_snps_NO, "bi_snps_NO.rds")
bi_snps_POS.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/bi_snps_POS.rds")
saveRDS(bi_snps_POS, "bi_snps_POS.rds")
bi_snps_MA.rds <- ("/rds/general/user/jr1320/home/../projects/neurogenomics-lab/live/Projects/dbsnp_non_bi_allelic_snps/results/dbSNP-151/bi_snps_MA.rds")
saveRDS(bi_snps_MA, "bi_snps_MA.rds")

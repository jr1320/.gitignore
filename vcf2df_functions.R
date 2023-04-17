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
  ## if (methods::is(vcf,"CollapsedVCF")) {
  ##   messager('Expanding VCF first, so number of rows may increase.',
  ##         v=verbose)
  ## Not all VCFs work with this function
  ##   vcf <- tryCatch({
  ##     VariantAnnotation::expand(x = vcf)
  ##   }, error = function(e){message(e); vcf})
  ##    print("vcf expanded. Line 30")
  ##  }
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
    if("ID" %in% colnames(df))
      n <- n[n!="ID"]
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
    #if you find anyhtin in the geno slot...
    if(length(tmp)>0){
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
  }
  #if you find anyhtin in the geno slot...
  if(length(tmp)>0)
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
  not_uniq_cols <- c()
  if(isTRUE(unique_rows) && isFALSE(unlist_cols)){
    # messager("Not using columns that are lists to check for unique rows",
    #  v=verbose)
    not_uniq_cols <- names(df)[ unlist(lapply(df, methods::is,"list")) ]
    #keep RS
    not_uniq_cols <- not_uniq_cols[!not_uniq_cols %in% 'RS']
  }
  if(isTRUE(unlist_cols)){
    unlist_dt(dt = df,
              verbose = verbose)
  }
  #### Remove duplicated rows ####
  if(isTRUE(unique_rows)){
    #first make RS col a single value instead of a list
    func <- function(x) as.character(paste0(unlist(x),collapse='-')[[1]])
    df[, RS:= lapply(RS,func),]
    df[, RS:= as.character(RS),]
    df <- drop_duplicate_rows(dt = df,
                              verbose = verbose,
                              ignore_cols = not_uniq_cols)
  }
  #### Report ####
  if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
  #messager("VCF data.table contains:",
  # formatC(nrow(df),big.mark = ","),"rows x",
  #  formatC(ncol(df),big.mark = ","),"columns.",
  #  v=verbose)
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
    # messager("Dropping",length(dups),"duplicate column(s).")
    dt[,  which(duplicated(names(dt))):= NULL]
  }
}
#remove_empty_cols:
remove_empty_cols <- function(sumstats_dt,
                              sampled_rows=NULL,
                              verbose=TRUE){
  #### Check for empty columns ####
  #  messager("Checking for empty columns.",v=verbose)
  empty_cols <- check_empty_cols(
    sumstats_dt = sumstats_dt,
    sampled_rows = sampled_rows,
    verbose = FALSE
  )
  #### Remove empty columns #####
  if(length(empty_cols)>0) {
    # messager("Removing",length(empty_cols),"empty columns.",v=verbose)
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
  # messager(length(empty_cols), "empty column(s) detected.",
  #      v=verbose)
  return(empty_cols)
}
print("Line 311")
#unlist_dt:
unlist_dt <- function(dt,
                      verbose=TRUE) {
  .SD <- NULL
  cols <- names(dt)[ unlist(lapply(dt, methods::is,"list")) ] 
  if(length(cols)>0){
    # messager("Unlisting",length(cols),"columns.",v=verbose)
    dt[ , (cols) := lapply(.SD,unlist), .SDcols = cols]
  }
}
print("Line 322")
#drop_duplicate_rows:
drop_duplicate_rows <- function(dt,
                                verbose=TRUE,
                                ignore_cols = c('CAF','TOPMED')){
  nrows_old <- nrow(dt)
  unique_cols <- colnames(dt)
  dt <- unique(dt,by=unique_cols[!(unique_cols %in% ignore_cols)])
  nrows_new <- nrow(dt)
  dropped <-  nrows_old-nrows_new
  if(dropped>0){
    #messager("Dropped",formatC(dropped,big.mark = ","),
    #"duplicate rows.",v=verbose)
  }
  return(dt)
}


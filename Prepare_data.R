#
# Symmetric sodium encryption and decryption of a datafile using cyphr.
# https://github.com/ropensci/cyphr
#
# idea: 
# 1) We generate a secure pwd and deliver safely to end users
# 3) APP automatically checks if a file ends with ".cyphr.bin" and decrypts it directly when reading (see decrypt_datafile)
library(cyphr)
library(jsonlite)
library(GenomicRanges)
library(kinship2)
library(dplyr)
library(tidyr)
pwd="Unlock_HICF2"

# load (gziped) data file
loadData <- function(dataF, append="", header=T, sep="\t", source=NA) {
  dataF = as.character(paste0(dataF, append))
  #message("Loading ", dataF)
  if ( endsWith(dataF, ".gz") ) {
    dat_file = gzfile(dataF)
    dat = read.table(dat_file, header=header, sep=sep, na.strings=c("na","NA","."), stringsAsFactors = F)
  } else {
    dat = read.table(dataF, header=header, sep=sep, na.strings=c("na","NA","."), stringsAsFactors = F)
  }
  if ( ! is.na(source)) {
    dat$source = source
  }
  return(dat)
}

# encrypt an R object
encrypt_datafile = function(myobj, outf, pwd) {
  k = cyphr::key_sodium(sodium::hash(charToRaw(pwd)))
  # encrypt data file
  result <- tryCatch({
    cyphr::encrypt(saveRDS(myobj, outf, compress = TRUE), k)
    return(1)
  }, error=function(cond) {
    #message("Encryption error. Unable to generate ", outf)
    if (file.exists(outf)) { file.remove(outf) }
    return(0)
  })
  
  #return (result)
}

# decrypt to an r object
decrypt_datafile = function(inf, outf, pwd) {
  #pwd = .rs.askForPassword("Enter password") # asks for pwd in RStudio
  k = cyphr::key_sodium(sodium::hash(charToRaw(pwd)))
  d = tryCatch({
    cyphr::decrypt(readRDS(inf), k)
  }, error=function(cond) {
    message("Could not decrypt! pwd wrong?")
    if (file.exists(outf)) { file.remove(outf) }
    return(NULL)
  })
  return (d)
}



args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2) {
  stop("Usage as follow:
       Prepare_data.R output_dir config_file.tsv/project_dir
       You can provide either a project folder OR a config file.
       If project folder is provided the script expects data to be distributed in standard sub-folders
       config_file is a tab-separated file containing the following columns:
       - var2reg: var2reg idx file (var2reg files paths and sampleIDs are read from this)
       - ROH: directory of ROH files, file names must match sample ID in var2reg idx (sampleID.txt)
       - ExpansionHunter: directory of ExpHunter json files, file names must match sample ID in var2reg idx (sampleID.json)", 
       call.=FALSE)
}

output_dir <- args[1]
if (file_test("-d", args[2])) {
  from_config = FALSE
  project_dir = args[2]
  
  message("No config file provided. Reading from standard folder")
  message("Output folder: ", output_dir)
  
} else if (file_test("-f", args[2])) {
  config <- read.table(args[2], sep="\t", header=T, stringsAsFactors = F)
  from_config = TRUE
  
  message("Reading from config file ", args[2])
  message("Output folder: ", output_dir)
} else {
  stop("Provided config file or project directory not found")
}

if (from_config==TRUE) {
  message("Reading from config not implemented yet")

} else {
  
  var2reg_dir <- paste0(project_dir,"/var2reg")
  roh_dir <- paste0(project_dir,"/ROH")
  roh_suffix <- "_ROH_full"
  ExpHunter_dir <- paste0(project_dir,"/Expansion_hunter")

  idx_file <- paste0(var2reg_dir, "/V2.var2reg.idx.tsv.gz")
  idx_df <- loadData(idx_file)
  message("Loaded var2reg idx file containing ", nrow(idx_df), " dataset")
  
  #samplesID <- gsub("V2\\.|\\.var2reg.vars.tsv.gz","",files)
  
  message("Loading and preparing data")
  saved_files <- 0
  failed_files <- 0
  good_peds <- 0
  pb <- txtProgressBar(min=0, max=nrow(idx_df),style = 3)
  for (n in 1:nrow(idx_df)) {
    setTxtProgressBar(pb,value = n)
    
    #convert idx df to list
    newlist <- as.list(idx_df[n,])
    
    #split samples ids into vectors
    newlist$all_samples <- strsplit(newlist$all_samples, ",")[[1]]
    newlist$affected_samples <- strsplit(newlist$affected_samples, ",")[[1]]
    newlist$unaffected_samples <- setdiff(newlist$all_samples, newlist$affected_samples)
    newlist$n_affected <- length(newlist$affected_samples)
    newlist$n_unaffected <- length(newlist$unaffected_samples)
    
    #Set some values for the app
    newlist$values_affected <- seq(1,newlist$n_affected)
    names(newlist$values_affected) <- seq(1,newlist$n_affected)
    newlist$values_affected <- c("NOT_ALLOWED" = 0, newlist$values_affected)
    newlist$values_unaffected <- seq(0, newlist$n_unaffected)
    
    #load gado data
    gado_score <- loadData(newlist$gado_file)
    newlist$gado90 <- quantile(gado_score$Zscore, .9)[[1]]
    
    #load comphet data
    comphet_df <- loadData(paste0(newlist$comphet_file,".gz"))
    comphet_df$ID <- paste0("CompHet_", comphet_df$rec_id)
    comphet_df$Class <- "PASS"
    newlist$comphet_df <- comphet_df
    
    #load variants data and set pop AD to zero when missing
    variants_df <- loadData(paste0(newlist$variant_file,".gz"))
    variants_df$max_pop_af[is.na(variants_df$max_pop_af)] <- 0
    variants_df[is.na(variants_df)] <- 0
    variants_df$Class <- "PASS"
    variants_ranges <- GRanges(seqnames = variants_df$chr, ranges = IRanges(variants_df$start, end=variants_df$end),ID=variants_df$rec_id)
    newlist$variants_df <- variants_df
    newlist$variants_ranges <- variants_ranges

    #Create segregation df
    seg_df1 <- as.data.frame(
      comphet_df %>% select(rec_id,num_aff) %>% 
        mutate(hom_aff=0, hom_unaff=0, het_aff=0, het_unaff=0) %>%
        dplyr::rename(comphet_aff = num_aff) )
    seg_df2 <- as.data.frame(
      variants_df %>% select(rec_id,hom_aff,hom_unaff,het_aff,het_unaff) %>%
        mutate(comphet_aff = 0) )
    newlist$segregation_df <- rbind(seg_df1, seg_df2)
    
    #Load genes table. Gene scores are moved to a separate table
    #Missing scores are set to zero
    genes_df <- loadData(paste0(newlist$gene_file,".gz"))
    genes_df$Class <- "PASS"
    genes_scores <- as.data.frame(genes_df %>% select(gene,gado_zscore,exomiser_gene_pheno_score,pLI_exac,pLI_gnomad,Class) %>% distinct() %>% arrange(desc(gado_zscore)))
    genes_df <- as.data.frame(genes_df %>% select(gene,inh_model,variants_n,variants,Class))
    
    genes_with_comphet <- as.data.frame(
      comphet_df %>% select(gene,rec_id) %>% group_by(gene) %>% 
        mutate(variants=paste(rec_id,collapse = ","), variants_n=n(), inh_model="comphet", Class="PASS") %>% 
        select(gene,variants,variants_n,inh_model,Class) %>% 
        distinct()
    )
    genes_df <- rbind(genes_df, genes_with_comphet)
    genes_df <- as.data.frame(genes_df %>% separate_rows(variants, sep=","))
    genes_scores$pLI_gnomad[is.na(genes_scores$pLI_gnomad)] <- 0
    genes_scores$pLI_exac[is.na(genes_scores$pLI_exac)] <- 0
    newlist$genes_df <- genes_df
    newlist$genes_scores <- genes_scores
    
    #Load Expansion Hunter jsons and ROH data
    newlist$ROH_data <- NULL
    newlist$ROH_ranges <- list()
    newlist$ExpHunter <- list()
    for (s in newlist$all_samples) {
      ROH_file <- paste0(roh_dir, "/",s,roh_suffix,".txt")
      ROH_df <- loadData(ROH_file)
      newlist$ROH_ranges[[s]] <- GRanges(
        seqnames = ROH_df$Chromosome, 
        ranges = IRanges(ROH_df$Start, end=ROH_df$End), 
        ID = paste0("ROH_", 1:nrow(ROH_df)),
        value = ROH_df$End - ROH_df$Start)
      
      ROH_df$Sample <- s
      newlist$ROH_data <- rbind(newlist$ROH_data,ROH_df)

      ExpHunter_file <- paste0(ExpHunter_dir, "/", s, ".json")
      newlist$ExpHunter[[s]] <- read_json(ExpHunter_file)
    }
    
    #Prepare ExpHunter table
    exphunter <- NULL
    for (sID in names(newlist$ExpHunter)) {
      #message(sID)
      locus_vars <- NULL
      for (gene in names(newlist$ExpHunter[[sID]]$LocusResults)) {
        #message(gene)
        single_vars <- NULL
        for (varID in names(newlist$ExpHunter[[sID]]$LocusResults[[gene]]$Variants)) {
          mydf <- as.data.frame(newlist$ExpHunter[[sID]]$LocusResults[[gene]]$Variants[[varID]])
          #colnames(mydf) <- col_names
          if (!is.null(mydf$Genotype)) {single_vars <- rbind(single_vars,mydf)}
        }
        single_vars$LocusId <- newlist$ExpHunter[[sID]]$LocusResults[[gene]]$LocusId
        single_vars$Coverage <- newlist$ExpHunter[[sID]]$LocusResults[[gene]]$Coverage
        locus_vars <- rbind(locus_vars,single_vars)
      }
      locus_vars$sampleID <- sID
      locus_vars$sex <- newlist$ExpHunter[[sID]]$SampleParameters$Sex
      exphunter <- rbind(exphunter, locus_vars)
    }
    exphunter$group <- "UNAFFECTED"
    exphunter$group[exphunter$sampleID %in% newlist$affected_samples] <- "AFFECTED"
    newlist$ExpHunter$PEDIGREE <- exphunter
    
    #Prepare ROH intersection for affected samples
    intersect_ROH <- newlist$ROH_ranges[[newlist$affected_samples[1]]]
    for (sID in newlist$affected_samples[-1]) {
      intersect_ROH <- GenomicRanges::intersect(intersect_ROH, newlist$ROH_ranges[[sID]])
    }
    newlist$ROH_ranges$AFFECTED_SHARED <- intersect_ROH
    newlist$ROH_ranges$AFFECTED_SHARED$value <- width(newlist$ROH_ranges$AFFECTED_SHARED)
    newlist$ROH_ranges$AFFECTED_SHARED$ID <- paste0("ROH_", 1:length(newlist$ROH_ranges$AFFECTED_SHARED))
    
    #Load ped file into kinship2 format
    ped_df <- loadData(newlist$pedigree_file, header=F)
    ped_df$V5[ped_df$V5==1] <- "male"
    ped_df$V5[ped_df$V5==2] <- "female"
    ped_df$V5[ped_df$V5==0] <- "unknown"
    ped_df$V6[ped_df$V6 == 0] <- 1
    
    tryCatch({
      fixed_ped <- with(ped_df, fixParents(id=V2, dadid = V3, momid=V4, sex=V5, missid=0))
      newlist$ped <- with(ped_df, pedigree(id=V2, dadid = V3, momid = V4, affected = V6, sex = V5, missid = 0))
      good_peds <- good_peds + 1
    }, error=function(cond) {
      newlist$ped <- NA
    })
      
    #Save encrypted object
    out_file <- paste0(output_dir, "/", newlist$pedigree, ".RData")
    save_results <- encrypt_datafile(newlist,outf=out_file,pwd = pwd)
    if (save_results == 1) {
      saved_files = saved_files + 1
    } else {
      failed_files = failed_files + 1
    }
  }
  
  message("\nProcess completed\n",saved_files, " saved\n", failed_files, " failed\n", good_peds, " correct pedigree loaded")
}


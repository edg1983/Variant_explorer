## DESCRIPTION ------------------
# Cohort varan post-process script
# Author: Edoardo Giacopuzzi

# Pre-process var2reg results, ROH data and Expansion Hunter data
# Save a list as single S3 object ready to be loaded into Variant Explorer app
# optparse lib must be installed in one of the locations accessible to R
# you can then specify an optional lib location to provide other libs

## DECLARE VARS ------------
libraries <- c("data.table","dplyr","GenomicRanges","jsonlite","kinship2", "optparse","Rlabkey","tidyr")
VERSION <- "v1.0"
labkey_url <- "https://labkey-embassy.gel.zone/labkey/"
gel_data_v <- "main-programme_v10_2020-09-03"
run_date <- format(Sys.time(),"%Y%m%d_%H%M%S")

## FUNCTIONS ---------------

# load (gziped) data file
loadData <- function(dataF, skipchar="#", header=T, sep="\t", source=NA) {
  #message("Loading ", dataF)
  if ( endsWith(dataF, ".gz") ) {
    dat_file = paste0('zgrep -v "', skipchar, '" ', dataF)
  } else {
    dat_file = paste0('grep -v "', skipchar, '" ', dataF)
  }
  dat = fread(cmd=dat_file, header=header, sep=sep, na.strings=c("na","NA","."), stringsAsFactors = F)
  if ( ! is.na(source)) {
    dat$source = source
  }
  return(dat)
}

# save R object
saveData = function(myobj, outf, mode="overwrite") {
  result <- tryCatch({
    saveRDS(myobj, outf, compress = TRUE)
    return(1)
  }, error=function(cond) {
    if (file.exists(outf)) { file.remove(outf) }
    return(0)
  })
}

# check if a lib is installed
InstalledPackage <- function(package) {
  available <- suppressMessages(suppressWarnings(
    sapply(
      package,
      require,
      quietly = TRUE,
      character.only = TRUE,
      warn.conflicts = FALSE
    )
  ))
  missing <- package[!available]
  if (length(missing) > 0)
    return(FALSE)
  return(TRUE)
}

# load libraries 
loadLibraries <- function(libs) {
  for (l in libs) {
    if (InstalledPackage(l)) {
      suppressMessages(suppressWarnings(
        sapply(
          l,
          library,
          quietly = TRUE,
          character.only = TRUE,
          warn.conflicts = FALSE
        )
      ))
    } else {
      stop(l, " library is not installed", call. = FALSE)
    }
  }
}

# check file exists and print message
checkFileExists <- function(files, mode="o", exit="warn") {
  for (file in files) {
    if (mode == "o") {
      if (!file.exists(file)) {
        if (exit == "warn") {
          warning(file, " not found. Skipped!",call. = FALSE, immediate. = TRUE)
          return(FALSE)
        } else if (exit == "stop") {
          stop(file, " not found",call. = FALSE)
        }
      } else {
        return(TRUE)
      }
    } else if (mode == "w") {
      if (file.exists(file)) {
        if (exit == "rename") {
          warning(file, " already present. Old file renamed to .", run_date, call. = FALSE, immediate. = TRUE)
          file.rename(file, paste0(file,".",run_date))
        } 
      }
    }
  }
}

# load idx files from config and make them to list
readFromConfig <- function(config, tag) {
  outlist <- NULL
  if (!is.null(config[[tag]]$index_file) & !config[[tag]]$index_file=="") {
    if (checkFileExists(c(config[[tag]]$index_file),"o","warn")) {
      df <- loadData(config[[tag]]$index_file, header = config[[tag]]$has_header)
      outlist <- makeNamedList(df, config[[tag]]$sampleID_col, config[[tag]]$file_col)
    } 
  }
  return(outlist)
}

# Make a named list from 2 cols in a df, one as names, one as values
makeNamedList <- function(df,names_col,values_col) {
  outlist <- as.list(df[[values_col]])
  names(outlist) <- df[[names_col]]
  return(outlist)
}

## INITIALIZE --------------
loadLibraries("optparse")

option_list <- list(
  make_option(c("-c","--config"), type = "character", default = NA,
              help = "Config file (JSON)"),
  make_option(c("-i","--index"), type = "character", default = NA,
              help = "var2reg index file (.tsv.gz)"),
  make_option(c("-d","--dataset_version"), type = "character", default = run_date,
              help = "A name for this version of the dataset. Default to the current date"),
  make_option(c("-o","--output"), type = "character", default = "processed_data",
              help = "output folder where processed data is saved"),
  make_option(c("-w","--overwrite"), action = "store_true", default=FALSE,
              help = "set this option to overwrite existing output files"),
  make_option(c("-k","--use_labkey"), action = "store_true", default=FALSE,
              help = "When running in GEL, set this option to get bam / roh files locations from LabKey"),
  make_option(c("-l","--lib_path"), type = "character", default = NA, 
              help = "additional path for R libraries")
)
args <- parse_args(OptionParser(option_list = option_list))

message(
  "=======================================
Cohort_VARAN POST PROCESSING ", VERSION, "\n",
  "----------------------------------------
var2reg index: ", args$index, "\n",
  "Config file: ", args$config, "\n",
  "Get from LabKey: ", args$use_labkey, "\n",
  "Output folder: ", args$output, "\n",
  "Dataset version: ", args$dataset_version, "\n",
  "======================================="
)

# If an additional lib path is specified this is added to libPaths
if (!is.na(args$lib_path)) {
  if (!file_test("-d", args$lib_path)) { 
    warning("You have specified a non existing folder as additional libpath", immediate. = T)
  } else {
    .libPaths(c(args$lib_path, .libPaths()))
  }
}
message("Loading libraries...")
loadLibraries(libraries)

# LabKey may be used when in GEL environment to get path of data files
# If you set --use_labkey, the script will create bam_file and roh_file tabs from labkey
# At the moment we read ExpHunter data from json while GEL only provides vcf, so this is skipped
if (args$use_labkey) {
  exphunter_files <- NULL
  message("LABKEY option active - Get BAM/ROH locations from LabKey DB")
  message("Eventual config file will be ignored")
  #loadLibraries("Rlabkey")
  labkey.setDefaults(baseUrl=labkey_url)
  message("Performing LabKey query...")
  query_bam <- 'SELECT Platekey, File_path FROM genome_file_paths_and_types 
                WHERE File_sub_type == "BAM"'
  
  suppressMessages(suppressWarnings(
  bam_df <- labkey.executeSql(folderPath=paste0("/main_programme/", gel_data_v),
                                 schemaName="lists",
                                 colNameOpt="rname",
                                 sql = query_bam,
                                 maxRows = 1e+06 )
  ))
  if (is.null(bam_files) | nrow(bam_files) == 0) {
    warning("Unable to get bam files from LabKey table. BAM and ROH files will not be loaded!", immediate. = T)
    bam_files <- NULL
    roh_files <- NULL
  } else {
    message(nrow(bam_files), " samples loaded from LabKey table")
    bam_files <- makeNamedList(bam_df, "platekey", "file_path")
    
    roh_df <- bam_files
    roh_df$file_path <- gsub("Assembly.*","",roh_df$file_path)
    roh_df$file_path <- paste0(roh_df$file_path, "Variations/", roh_df$platekey, ".ROH.bed")
    roh_files <- makeNamedList(roh_df, "platekey", "file_path")
    
    roh_cols <- list(chrom = "V1", start = "V2", stop = "V3")
  }
}

if (args$overwrite) { write_mode <- "overwrite" } else { write_mode <- "rename" }

if (is.na(args$index)) { stop("You must specify an index file") }

dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
if (!file_test("-d", args$output)) { 
  stop("Output folder: ", args$output, "\nThe folder does not exists and cannot be created", call. = FALSE)
}

checkFileExists(args$index,"o", "stop")

releaseID <- args$dataset_version
output_dir <- args$output
idx_file <- args$index 
  
## READ DATA FROM CONFIG --------------
#Config is loaded only if use labkey is false otherwise it is ignored
if (!args$use_labkey) {
  if (is.na(args$config)) {
    warning("You have specified neither LabKey or a config file. No BAM / ROH / ExpHunter files will be loaded",
            immediate. = T) 
    bam_files <- NULL
    roh_files <- NULL
    exphunter_files <- NULL
  } else {
    message("Loading BAM / ROH locations from files specified in config file")
    file_exists <- checkFileExists(args$config,"o","stop")
    config <- read_json(args$config)
    bam_files <- readFromConfig(config,"BAM")
    roh_files <- readFromConfig(config,"ROH")
    exphunter_files <- readFromConfig(config,"EXPHUNTER")
    roh_cols <- list(
      chrom = config$ROH$ROH_file_structure$chr_col,
      start = config$ROH$ROH_file_structure$start_col,
      stop = config$ROH$ROH_file_structure$stop_col
    )
  }
}
# TODO make segregation col names configurable

## PROCESS DATA ----------------
start_time <- Sys.time()
idx_df <- loadData(idx_file)
saved_files <- 0
failed_files <- 0
good_peds <- 0
total <- nrow(idx_df)
message("Loaded var2reg idx file containing ", total, " dataset\n")
message("#### START PROCESSING ####")

for (n in 1:nrow(idx_df)) {
  #convert idx file line to list
  newlist <- as.list(idx_df[n,])
  newlist$releaseID <- releaseID
  message(Sys.time(), " #### file ",n, " ", newlist$pedigree, " --- ", round((n/total) * 100, 2), " %")
  
  #Set output filename
  out_file <- paste0(output_dir, "/", newlist$pedigree, ".RData")
  checkFileExists(c(out_file),"w", write_mode)
  
  #Split samples IDs into vectors 
  newlist$all_samples <- strsplit(newlist$all_samples, ",")[[1]]
  newlist$affected_samples <- strsplit(newlist$affected_samples, ",")[[1]]
  
  #Check mandatory input files were present
  checkFileExists(c(newlist$gado_file,
                    newlist$gene_file,
                    newlist$comphet_file,
                    newlist$variant_file,
                    newlist$pedigree_file,
                    newlist$known_variant_file
                    ),"o", "stop")
  
  #GADO data ----------
  gado_score <- loadData(newlist$gado_file)
  newlist$gado90 <- quantile(gado_score$Zscore, .9)[[1]]
  
  #COMPHET data -------
  comphet_df <- loadData(newlist$comphet_file)
  #if (nrow(comphet_df) == 0) {next} #skip loading data if comphet empty
  comphet_df$ID <- paste0("CompHet_", comphet_df$rec_id)
  comphet_df$Class <- "PASS"
  newlist$comphet_df <- comphet_df
  
  #VARIANTS data ------
  variants_df <- loadData(newlist$variant_file)
  variants_df$Class <- "PASS"
  variants_df$internal_id <- paste(variants_df$chr,variants_df$start,variants_df$end,variants_df$ref,variants_df$alt, sep="_")
  variants_ranges <- GRanges(seqnames = variants_df$chr, ranges = IRanges(variants_df$start, end=variants_df$end),ID=variants_df$rec_id)
  newlist$variants_df <- variants_df
  newlist$variants_ranges <- variants_ranges
  
  #retain only samples with genotypes
  #it is possible that additional samples are present in the PED but not in the VCF
  #we want samples IDs only for samples that are actually in the VCF
  all_samples_GT <- paste("GT", newlist$all_samples, sep ="_")
  
  tot_vars <- nrow(newlist$variants_df)
  samples_with_GT <- NULL 
  for (sID in all_samples_GT) {
    if (sum(is.na(newlist$variants_df[[sID]])) < tot_vars) {
      samples_with_GT <- c(samples_with_GT, gsub("GT_","",sID))
    }
  }
  newlist$all_samples <- samples_with_GT
  newlist$affected_samples <- intersect(newlist$affected_samples, newlist$all_samples)
  
  newlist$unaffected_samples <- setdiff(newlist$all_samples, newlist$affected_samples)
  newlist$n_affected <- length(newlist$affected_samples)
  newlist$n_unaffected <- length(newlist$unaffected_samples)
  
  #These are values used by Variant Explorer for segregation filters
  newlist$values_affected <- seq(0,newlist$n_affected)
  newlist$values_unaffected <- seq(0, newlist$n_unaffected)

  #fill hom_ / het_ columns with genotypes counts
  #var2reg counts het / hom only if they are above quality threshold
  #we want full counts for the app segregation filter instead, since GQ filter can be set separately
  affected_cols <- which(colnames(newlist$variants_df) %in% paste("GT",newlist$affected_samples, sep="_"))
  newlist$variants_df$het_aff <- rowSums(newlist$variants_df[,..affected_cols] == 1,na.rm = T)
  newlist$variants_df$hom_aff <- rowSums(newlist$variants_df[,..affected_cols] == 2,na.rm = T)
  if (length(newlist$unaffected_samples) > 0) {
    unaffected_cols <- which(colnames(newlist$variants_df) %in% paste("GT",newlist$unaffected_samples, sep="_"))
    newlist$variants_df$het_unaff <- rowSums(newlist$variants_df[,..unaffected_cols] == 1,na.rm = T)
    newlist$variants_df$hom_unaff <- rowSums(newlist$variants_df[,..unaffected_cols] == 2,na.rm = T)
  } else {
    newlist$variants_df$het_unaff <- 0
    newlist$variants_df$hom_unaff <- 0
  }
  
  #SEGREGATION data ------------
  if (nrow(newlist$comphet_df) > 0) {
    seg_df1 <- as.data.frame(
      newlist$comphet_df %>% select(rec_id,num_aff) %>% 
      mutate(hom_aff=0, hom_unaff=0, het_aff=0, het_unaff=0, sup_dnm=0) %>%
      dplyr::rename(comphet_aff = num_aff) )
  } else {
    seg_df1 <- NULL
  }
  seg_df2 <- as.data.frame(
      newlist$variants_df %>% select(rec_id,hom_aff,hom_unaff,het_aff,het_unaff,sup_dnm) %>%
      mutate(comphet_aff = 0) )
  newlist$segregation_df <- rbind(seg_df1, seg_df2)
  
  #GENES data -------------------- 
  #Genes scores are moved to a separate table (genes_scores)
  genes_df <- loadData(paste0(newlist$gene_file))
  genes_df$Class <- "PASS"
  genes_scores <- as.data.frame(genes_df %>% 
    select(gene,gado_zscore,exomiser_gene_pheno_score,pLI_exac,pLI_gnomad,GDI_phred,EDS,RVIS,Class) %>% 
      distinct() %>% 
      arrange(desc(gado_zscore)))
  genes_df <- as.data.frame(genes_df %>% select(gene,inh_model,variants_n,variants,Class))
  
  #if there are comphet, we copy the genes and vars to genes table to ensure all genes are represented
  if (nrow(comphet_df) > 0) {
    genes_with_comphet <- as.data.frame(
      comphet_df %>% select(gene,rec_id) %>% group_by(gene) %>% 
        mutate(variants=paste(rec_id,collapse = ","), variants_n=n(), inh_model="comphet", Class="PASS") %>% 
        select(gene,variants,variants_n,inh_model,Class) %>% 
        distinct() )
  } else {
    genes_with_comphet <- NULL
  }
  
  genes_df <- rbind(genes_df, genes_with_comphet)
  genes_df <- as.data.frame(genes_df %>% separate_rows(variants, sep=","))
  newlist$genes_df <- genes_df
  newlist$genes_scores <- genes_scores
  
  #PED files ---------------
  #Read PED into kinship2 format, singleton cannot be displayed
  ped_df <- loadData(newlist$pedigree_file, header=F)
  ped_df$V5[ped_df$V5==1] <- "male"
  ped_df$V5[ped_df$V5==2] <- "female"
  ped_df$V5[ped_df$V5==0] <- "unknown"
  ped_df$V6[ped_df$V6 == 0] <- 1
  
  #Load ped file and auto fix, if failed the ped cannot be displayed
  tryCatch({
    fixed_ped <- with(ped_df, fixParents(id=V2, dadid = V3, momid=V4, sex=V5, missid=0))
    fixed_ped <- merge(fixed_ped, ped_df[,c("V2","V6")], by.x="id", by.y="V2", all.x=T)
    fixed_ped$V6[is.na(fixed_ped$V6)] <- 1
    newlist$ped <- with(fixed_ped, pedigree(id=id, dadid = dadid, momid = momid, sex = sex, affected = V6, missid = 0))
    good_peds <- good_peds + 1
  }, error=function(cond) {
    newlist$ped <- NA
  })
  
  #KNOWN VARS data ----------------
  newlist$known_vars <- loadData(newlist$known_variant_file)
  newlist$known_vars <- newlist$known_vars %>% separate_rows(known_ids, sep=",")
  newlist$known_vars$internal_id <- paste(newlist$known_vars$chr,newlist$known_vars$start,newlist$known_vars$end,newlist$known_vars$ref,newlist$known_vars$alt, sep="_")
  newlist$known_clinvar <- newlist$known_vars[grep("CV[0-9]+",newlist$known_vars$known_ids,perl = T),]
  newlist$known_cosmic <- newlist$known_vars[grep("COSV[0-9]+",newlist$known_vars$known_ids,perl = T),]
  
  #ROH data --------------
  newlist$ROH_data <- NULL
  newlist$ROH_ranges <- list()
  
  if (inherits(roh_files, "list")) {
    for (s in newlist$all_samples) {
      if (!is.null(roh_files[[s]])) {
        ROH_file <- roh_files[[s]]
        if (checkFileExists(ROH_file,"o","warn")) {
          ROH_df <- loadData(ROH_file)
          newlist$ROH_ranges[[s]] <- GRanges(
            seqnames = ROH_df[[roh_cols$chrom]], 
            ranges = IRanges(ROH_df[[roh_cols$start]], end=ROH_df[[roh_cols$stop]]), 
            ID = paste0("ROH_", 1:nrow(ROH_df)),
            value = ROH_df[[roh_cols$stop]] - ROH_df[[roh_cols$start]])
          
          ROH_df$Sample <- s
          newlist$ROH_data <- rbind(newlist$ROH_data,ROH_df)
        }
      }
    }
    
    #Prepare ROH intersection for affected samples
    if (newlist$n_affected > 1) {
      affected_withROH <- intersect(newlist$affected_samples, names(newlist$ROH_ranges))
      intersect_ROH <- newlist$ROH_ranges[[affected_withROH[1]]]
      if (!is.null(affected_withROH)) {
        for (sID in affected_withROH[-1]) {
          intersect_ROH <- GenomicRanges::intersect(intersect_ROH, newlist$ROH_ranges[[sID]])
        }
        newlist$ROH_ranges$AFFECTED_SHARED <- intersect_ROH
        newlist$ROH_ranges$AFFECTED_SHARED$value <- width(newlist$ROH_ranges$AFFECTED_SHARED)
        newlist$ROH_ranges$AFFECTED_SHARED$ID <- paste0("ROH_", 1:length(newlist$ROH_ranges$AFFECTED_SHARED))
      }
    }
  }
  
  #EXPANSION HUNTER (jsons) --------------
  newlist$ExpHunter <- NULL
  
  if (inherits(exphunter_files, "list")) {
    for (s in newlist$all_samples) {
      if (!is.null(exphunter_files[[s]])) {
        ExpHunter_file <- exphunter_files[[s]]
        if (checkFileExists(ExpHunter_file,"o","warn")) {
          newlist$ExpHunter[[s]] <- read_json(ExpHunter_file)
        }
      }
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
        if (!is.null(single_vars)) {
          single_vars$LocusId <- newlist$ExpHunter[[sID]]$LocusResults[[gene]]$LocusId
          single_vars$Coverage <- newlist$ExpHunter[[sID]]$LocusResults[[gene]]$Coverage
          locus_vars <- rbind(locus_vars,single_vars)
        }
      }
      locus_vars$sampleID <- sID
      locus_vars$sex <- newlist$ExpHunter[[sID]]$SampleParameters$Sex
      exphunter <- rbind(exphunter, locus_vars)
    }
    exphunter$group <- "UNAFFECTED"
    exphunter$group[exphunter$sampleID %in% newlist$affected_samples] <- "AFFECTED"
    newlist$ExpHunter$PEDIGREE <- exphunter
  }

  #BAM FILES locations --------------
  newlist$BAM_files <- NULL
  if(inherits(bam_files, "list")) {
    newlist$bam_files <- bam_files[newlist$all_samples]
  }
  
  #SAVE DATA --------------------
  #Save the list containig processed data into an RDS object
  save_results <- saveData(newlist,outf=out_file)
  if (save_results == 1) {
    saved_files = saved_files + 1
  } else {
    failed_files = failed_files + 1
  }
}

## CLOSE MESSAGE ----------------------
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time <- paste0(round(elapsed_time,2), " ", units(elapsed_time))
message(
"\n=======================================
FINISHED PROCESSING in ", elapsed_time, "\n",
"----------------------------------------
N input cases from index: ", total, "\n",
"N saved: ", saved_files, "\n",
"N failed: ", failed_files, "\n",
"N with ped: ", good_peds, "\n",
"======================================="
)
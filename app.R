# VARIANT EXPLORER
# Author: Edoardo Giacopuzzi
# Contact: edoardo.giacopuzzi@well.ox.ac.uk
# Explore and filter annotated variants from var2reg

# Input are (encrypted) RData objects created with Prepare_Rdata_object.R     
# Each object contains data from VARAN V2 and var2reg, ROH data and Exp Hunter data

## CONSTANTS --------------------
APP_VERSION <- "1.2.4"
vis_cols <- c("gene","chr","start","end","ref","alt","var_type","consequence","known_ids","max_pop_af","cohort_af")

## FUNCTIONS --------------------------
`%nin%` = Negate(`%in%`)

combine_df <- function (x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

# load libraries 
loadLibraries <- function(libs) {
  message("Loading libraries...")
  for (l in libs) {
    suppressMessages(suppressWarnings(
      sapply(
        l,
        library,
        quietly = TRUE,
        character.only = TRUE,
        warn.conflicts = FALSE
      )
    ))
  }
}

readGMT <- function(input_file) {
  out_list <- list()
  con  <- file(input_file, open = "r")
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myvector <- (strsplit(oneLine, "\t"))
    out_list[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
  }
  close(con) 
  return(out_list)
}

getResource <- function(input_file, parent_dir=resource_dir) {
  resource_file <- paste0(parent_dir,"/",input_file)
  return(resource_file)
}

decrypt_datafile = function(inf, pwd) {
  #pwd = .rs.askForPassword("Enter password") # asks for pwd in RStudio
  k = cyphr::key_sodium(sodium::hash(charToRaw(pwd)))
  d = tryCatch({
    cyphr::decrypt(readRDS(inf), k)
  }, error=function(cond) {
    #RV$decrypt_status = "Could not decrypt! pwd wrong?"
    #if (file.exists(outf)) { file.remove(outf) }
    return(0)
  })
  return (d)
}

computeNormZ <- function(gene,score) {
  if (gene %in% names(genes_dist)) {
    return(genes_dist[[gene]](score))
  } else {
    return(NA)
  }
}

#Retrieve clinical impact and disease from ClinVar
getClinvarInfo <- function(ids, prefix="CV") {
  #pb <- txtProgressBar(min = 0, max = length(ids),style = 3)
  clinvar_info <- list()
  ids <- gsub(prefix,"",ids)
  #i = 0
  for (n in ids) {
    #i = i + 1
    #setTxtProgressBar(pb, i)
    myvar <- fromJSON(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=",n,"&retmode=json"))
    clinvar_info$significance <- c(clinvar_info$significance,
                                   myvar$result[[as.character(n)]]$trait_set$trait_name)
    clinvar_info$phenotype <- c(clinvar_info$phenotype,
                                myvar$result[[as.character(n)]]$clinical_significance$description)
  }
  clinvar_info$significance <- paste(clinvar_info$significance, collapse="; ")
  clinvar_info$phenotype <- paste(clinvar_info$phenotype, collapse="; ")
  return(clinvar_info)
}

makeIGVxml <- function(region, affected_samples, unaffected_samples, VCF_files, BAM_files, SV_files) {
  header<- paste('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
                 paste0('<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" locus="', region, '">'),
                 '<Resources>', sep="\n")
  
  BAM_panels <- NULL
  VCF_panels <- NULL
  SV_panels <- NULL
  for (sample in affected_samples) {
    BAM_panel <- paste(paste0('<Panel name="', sample,'_panel">'),
                   paste0('<Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" id="', BAM_files[[sample]], '_coverage" name="', sample, ' Coverage" snpThreshold="0.2" visible="true"/>'),
                   paste0('<Track clazz="org.broad.igv.sam.SpliceJunctionTrack" id="', BAM_files[[sample]], '_junctions" visible="false"/>'),
                   paste0('<Track clazz="org.broad.igv.sam.AlignmentTrack" id="', BAM_files[[sample]], '" name="', sample, ' - affected" visible="true"/>'),
                   '</Panel>',
                   sep="\n")
    VCF_panel <- paste(paste0('<Panel name="', sample,' VCF">'),
                       paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', 
                             VCF_files$single_vcf[[sample]], 
                             '" name="', sample, 'VCF - affected" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>'),
                       '</Panel>',
                       sep="\n")
    SV_panel <- paste(paste0('<Panel name="', sample,' SV VCF">'),
                       paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', 
                              SV_files$single_vcf[[sample]], 
                              '" name="', sample, 'SV - affected" visible="true"/>'),
                       '</Panel>',
                       sep="\n")
    BAM_panels <- c(BAM_panels, BAM_panel)
    VCF_panels <- c(VCF_panels, VCF_panel)
    SV_panels <- c(SV_panels, SV_panel)
  }
  for (sample in unaffected_samples) {
    BAM_panel <- paste(paste0('<Panel name="', sample,'_panel">'),
                   paste0('<Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" id="', BAM_files[[sample]], '_coverage" name="', sample, ' Coverage" snpThreshold="0.2" visible="true"/>'),
                   paste0('<Track clazz="org.broad.igv.sam.SpliceJunctionTrack" id="', BAM_files[[sample]], '_junctions" visible="false"/>'),
                   paste0('<Track clazz="org.broad.igv.sam.AlignmentTrack" id="', BAM_files[[sample]], '" name="', sample, ' - unaffected" visible="true"/>'),
                   '</Panel>',
                   sep="\n")
    VCF_panel <- paste(paste0('<Panel name="', sample,' VCF">'),
                       paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', 
                              VCF_files$single_vcf[[sample]], 
                              '" name="', sample, 'VCF - unaffected" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>'),
                       '</Panel>',
                       sep="\n")
    SV_panel <- paste(paste0('<Panel name="', sample,' SV VCF">'),
                      paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', 
                             SV_files$single_vcf[[sample]], 
                             '" name="', sample, 'SV - unaffected" visible="true"/>'),
                      '</Panel>',
                      sep="\n")
    BAM_panels <- c(BAM_panels, BAM_panel)
    VCF_panels <- c(VCF_panels, VCF_panel)
    SV_panels <- c(SV_panels, SV_panel)
  }
  BAM_panels <- paste(BAM_panels, collapse="\n")
  VCF_panels <- paste(VCF_panels, collapse="\n")
  SV_panels <- paste(SV_panels, collapse="\n")
  
  resource_bam <- paste('<Resource path="', unlist(BAM_files), '"/>', collapse="\n", sep="")
  if (is.null(VCF_files$family_vcf)) {
    resource_vcf <- paste('<Resource path="', unlist(VCF_files$single_vcf), '"/>', collapse="\n", sep="")
  } else {
    resource_vcf <- paste0('<Resource path="', VCF_files$family_vcf, '"/>')
    VCF_panels <- paste('<Panel name="VCF_panel">',
                        paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', VCF_files$family_vcf, '" name="Family VCF" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>'),
                        '</Panel>',
                        sep="\n")
  }
  if (is.null(SV_files$family_vcf)) {
    resource_sv <- paste('<Resource path="', unlist(SV_files$single_vcf), '"/>', collapse="\n", sep="")
  } else {
    resource_sv <- paste0('<Resource path="', SV_files$family_vcf, '"/>')
    SV_panels <- paste('<Panel name="SV_panel">',
                       paste0('<Track clazz="org.broad.igv.variant.VariantTrack" id="', SV_files$family_vcf, '" name="SV family VCF" visible="true"/>'),
                       '</Panel>',
                       sep="\n")
  }
  
  close <- paste('<Panel height="61" name="FeaturePanel" width="1235">',
                 '<Track clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" visible="true"/>',
                 '<Track clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;426.0;255,255,255;0,0,178" fontSize="10" height="35" id="hg38_genes" name="Gene" visible="true"/>',
                 '</Panel>',
                 '</Session>', 
                 sep="\n")
  
  xml_text <- paste(header, 
                    resource_bam,
                    resource_vcf, 
                    resource_sv, 
                    '</Resources>', 
                    VCF_panels, 
                    SV_panels, 
                    BAM_panels,
                    close,
                    sep="\n")
  return(xml_text)
}

## CHECK AND LOAD PACKAGES ----------------
# Install missing packages
# pkgs "doParallel", "foreach" needed for cohort mode, now disabled
packages.cran <- c("tibble","R.utils","RSQLite","data.table","cyphr","shiny","shinyBS","stringr", "DT", "dplyr", "plotly", "kinship2", "tidyr", "shinydashboard", "gridExtra", "ggplot2", "jsonlite", "ontologyIndex")
new.packages <- packages.cran[!(packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages)) {install.packages(new.packages)}

#Packages from Bioconductor
packages.bioc <- c("GenomicRanges")
if(!("GenomicRanges" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install("GenomicRanges")
}

#Packages from github
packages.github <- c("scattermore","shinycssloaders")
if(!("scattermore" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
    devtools::install_github('exaexa/scattermore', upgrade = FALSE)
}
if(!("shinycssloaders" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
    devtools::install_github('daattali/shinycssloaders',upgrade = FALSE)
}

# Load libraries
#options(repos = BiocManager::repositories())
#getOption("repos")
loadLibraries(c(packages.cran, packages.bioc, packages.github))

## LOAD MODULES -----------------------
source("modules/plotModule.R")
source("modules/downloadModule.R")
source("modules/segregationModule.R")
source("modules/intersectBedModule.R")
source("modules/GQModule.R")
source("modules/filtersModule.R")
source("modules/SQliteModule.R")

## ENVIRONMENT CONFIG ----------------------------
## Read config files
app_settings <- read_json("App_configuration.json")
filters_settings <- read_json("Filters_settings.json")

#Set data dirs from config
data_dir <- app_settings$data_dir
resource_dir <- app_settings$resource_dir
user_dir <- app_settings$user_dir
Coverage_dir <- gsub("@resource_dir", resource_dir, app_settings$coverage_dir)
PanelApp_dir <- gsub("@resource_dir", resource_dir, app_settings$PanelApp_dir)
GeneLists_dir <- gsub("@resource_dir", resource_dir, app_settings$GeneLists_dir)
HPO_dir <- gsub("@resource_dir", resource_dir, app_settings$HPO_dir)
GREENDB_file <- gsub("@resource_dir", resource_dir, app_settings$GREENDB)
compute_cohort_GADO <- app_settings$compute_cohort_GADO

##Set axes options for plots
plot_axes <- list()
for (n in names(app_settings$plot_axes)) {
  plot_axes[[n]] <- unlist(app_settings$plot_axes[[n]], use.names = T)
}

##Plots formatting styles
format1 <- theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
                 axis.text.y = element_text(size=12))

##consequence groups for various variant categories
reg_vars <- app_settings$var_groups$consequence_groups$reg_vars 
exonic_vars <- app_settings$var_groups$consequence_groups$exonic_vars

##var_type groups for various variant categories
sv_vars <- app_settings$var_groups$var_type_groups$sv_vars
small_vars <- app_settings$var_groups$var_type_groups$small_vars

##Set segregation columns names
segregation_cols <- unlist(app_settings$segregation_cols, use.names = T)

##Set 3-letter to 1-letter aa codes
aa_codes <- list(
  Ala= "A",
  Arg= "R",
  Asn= "N",
  Asp= "D",
  Cys= "C",
  Glu= "E",
  Gln= "Q",
  Gly= "G",
  His= "H",
  Ile= "I",
  Leu= "L",
  Lys= "K",
  Met= "M",
  Phe= "F",
  Pro= "P",
  Ser= "S",
  Thr= "T",
  Trp= "W",
  Tyr= "Y",
  Val= "V"
)

##Set reactive objects
RV <- reactiveValues(
  cohort_files = NULL,
  notifications = list(),
  tasks = list(),
  data = 0,
  custom_genes = data.frame(gene=character(), source=character(), stringsAsFactors = F),
  customBed_ranges = FALSE,
  filters_summ_genes = data.frame(),
  filters_summ_vars = data.frame(),
  selected_vars_region = "NONE",
  accepted_reg_db = NULL,
  cov_file = NULL,
  saved_vars = NULL,
  filtered_applied = FALSE,
  messages = list(userguide = messageItem(from="Variant Explorer",
                                          message = "Variant explorer user guide", 
                                          href="https://variant-explorer.readthedocs.io/en/latest") )
)

RV_cohort <- reactiveValues()

##Put X scroll bar on top for data tables
css <- HTML(
  "#vars_results_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
  transform:rotateX(180deg);
  }
  #vars_results_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
  transform:rotateX(180deg);
  }
  #comphet_results_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
  transform:rotateX(180deg);
  }
  #comphet_results_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
  transform:rotateX(180deg);
  }
  #cohort_vars_df > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
  transform:rotateX(180deg);
  }
  #cohort_comphet_df > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
  transform:rotateX(180deg);
  }"
)

#Get list of data objects in the data directory
files <- list.files(data_dir, pattern = ".RData")
samplesID <- gsub("\\.RData[.enc]*","",files, perl = T)

## LOAD SUPPORTING FILES ------------------------
##Load hg38 chrom sizes
message("Load chromosome size")
chr_sizes_file <- getResource("hg38_chromSizes.txt")
chr_sizes <- read.table(chr_sizes_file,header=F, sep="\t", stringsAsFactors = F)
genome_size <- sum(as.numeric(chr_sizes$V2))

##Load gene names and IDs
message("Load gene IDs")
genes_info_file <- getResource("/hgnc_complete_set.txt")
genes_info <- read.table(genes_info_file, header=T, sep="\t", stringsAsFactors = F)

##Load genes coordinates
message("Load genes coordinates")
genes_bed_file <- gzfile(getResource("gencode.v33.basic.genes.bed.gz"))
genes_bed <- read.table(genes_bed_file, header=F, sep="\t",stringsAsFactors = F)

##Load pathways and GO groups (suppose .gmt file from MSigDB)
message("Load GO and pathways annotations")
gene_anno <- list()
gene_anno[["pathways"]] <- readGMT(getResource("c2.cp.v7.0.symbols.gmt"))
gene_anno[["GO_BP"]] <- readGMT(getResource("c5.bp.v7.0.symbols.gmt"))
gene_anno[["GO_MF"]] <- readGMT(getResource("c5.mf.v7.0.symbols.gmt"))
gene_anno[["GO_CC"]] <- readGMT(getResource("c5.cc.v7.0.symbols.gmt"))

## Load GTeX median expression
message("Load GTeX data")
GTeX_file <- gzfile(getResource("GTEx_v8_median_TPM.gct.gz"))
GTeX_data <- read.table(GTeX_file,sep="\t",header=T, stringsAsFactors = F)

## Load PanelApp index table
#PanelApp index table contains an id column with panel id number
#Each panel file is then names GRCh38_Panel_[id].bed
message("Load PanelApp panels")
PanelApp_file <- getResource("PanelApp_index.tsv")
PanelApp_data <- read.table(PanelApp_file, sep="\t", header = T, stringsAsFactors = F)

PanelApp_genes <- data.frame(entity_name=character(), confidence_level=character(), panel_idx=character(), stringsAsFactors = F)

n <- 0
for (id in PanelApp_data$id) {
    filename <- paste0(PanelApp_dir,"/GRCh38_Panel_", id, ".bed")
    if (file_test("-f", filename)) {
      mydf <- read.table(filename, sep="\t", header=T, stringsAsFactors = F, comment.char = "@")
      if (nrow(mydf) > 0) {
        n <- n + 1
        mydf$panel_idx <- id
        mydf <- mydf[,c("entity_name","confidence_level","panel_idx")]
        PanelApp_genes <- rbind(PanelApp_genes, mydf)
      }
    }
}
message(n, " panels loaded")

## Load additional genes list
message("Load gene lists")
#Expect a GeneLists_index.tsv file in the genelist folder
#This file must contain an id column and this ids mus match lists file name
#Example: list1.tsv, id=list1
geneLists_file <- getResource("GeneLists_index.tsv")
geneLists_data <- read.table(geneLists_file, sep="\t", header = T, stringsAsFactors = F)
geneLists_data <- geneLists_data[order(geneLists_data$id),]

geneLists_genes <- data.frame(entity_name=character(), genelist_idx=character(), stringsAsFactors = F)

n <- 0
for (id in geneLists_data$id) {
    n <- n + 1
    filename <- paste0(GeneLists_dir,"/", id, ".tsv")
    mydf <- data.frame(entity_name=scan(filename, what="", sep="\n",quiet = T), genelist_idx = id, stringsAsFactors = F)
    geneLists_genes <- rbind(geneLists_genes, mydf)
}
message(n, " gene lists loaded")

## Load ClinVar pathogenic genes
# Tab-separated file with gene symbol and list of diseases
message("Load ClinVar data")
ClinVar_file <- getResource("clinvar_path_likelypath_20190704.tsv")
ClinVar_genes <- read.table(ClinVar_file, sep="\t", header = T, stringsAsFactors = F)

##Load GADO distribution
message("Load HICF2 GADO distribution")
GADO_file <- getResource("GADO_distribution.tsv.gz")
if (file.exists(GADO_file)) {
  gado_distribution <- fread(GADO_file, sep="\t", header=T)
} else {
  gado_distribution <- NULL
}

#Possible future implementations will include cohort normalized GADO
#However this is quite slow and maybe not that useful
#message("Load cohort normalized GADO score")
#GADO_file <- getResource("GADO_cohortNormalized_score.RData")
#if (file.exists(GADO_file)) {
#  genes_dist <- readR`DS(GADO_file)
#} else {
#  message("Cohort normalized GADO file not found. Generating a new file...")
#  genes_dist <- gado_distribution %>% group_by(Hgnc) %>% group_map(~ ecdf(.x$Zscore))
#  names(genes_dist) <- levels(as.factor(gado_distribution$Hgnc))
#  saveRDS(genes_dist, file=GADO_file)  
#}

##Load HPO data
message("Load HPO profiles")
HPO_obo <- get_OBO(paste0(HPO_dir, "/hp.obo"))
HPO_obo <- as.data.frame(HPO_obo$name)
HPO_obo$HPO_id <- rownames(HPO_obo)
colnames(HPO_obo)[1] <- "HPO_name"
HPO_genes <- read.table(paste0(HPO_dir, "/genes_to_phenotype.txt"), sep="\t", header=F, stringsAsFactors = F) %>% select(V2,V3,V4)
colnames(HPO_genes) <- c("gene","HPO_id","HPO_name")
cohort_HPO_file <- getResource("cohort_HPO_terms.tsv")
if (file.exists(cohort_HPO_file)) {
  cohort_HPO <- read.table(cohort_HPO_file, sep="\t", header=T, stringsAsFactors = F)
  cohort_HPO <- cohort_HPO %>% separate_rows(HPO, sep=",")
} else {
  cohort_HPO <- NULL
}

#Load previously saved vars if present
saved_vars <- paste0(user_dir, "/Saved_variants.csv")
if ( file.exists(saved_vars) ) {
  message("Found previously saved vars")
  RV$saved_vars <- read.table(saved_vars, sep=",",header=T, stringsAsFactors = F, as.is=T)
}

## USER INTERFACE ------------------------
ui <- dashboardPage(
    dashboardHeader(
        title = paste0("VarExplorer v", APP_VERSION),
        dropdownMenuOutput("MessageMenu"),
        dropdownMenuOutput("NotificationMenu")
    ),
    dashboardSidebar(
        #drop-down list of cases
        selectInput("CaseCode", h3("Case code:"), choices = sort(base::unique(samplesID))),
        
        #Password and decrypt
        passwordInput("pwd","Password:"),
        actionButton("decrypt_button",label = "Load data"),
        
        #Custom genomic regions
        fileInput(inputId = "custom_bed",label = "Region BED:", multiple=FALSE, accept=".bed", placeholder = "region bed file"),
        
        sidebarMenu(id = "tabs",
            menuItem("Variants overview", tabName = "overview", icon = icon("th")),
            menuItem("Create custom genes list", tabName = "custom_genes_list_builder", icon = icon("th")),
            menuItem("Filters settings", icon = icon("th"),
              menuSubItem("Overview", tabName = "filters_overview_tab"),
              menuSubItem("Variants", tabName = "variants_filters_tab"),
              menuSubItem("Segregation", tabName = "segregation_filters_tab"),
              menuSubItem("Genes", tabName = "genes_filters_tab"),
              menuSubItem("Regions and ROH", tabName = "regions_filters_tab")
            ),
            menuItem("Filter explorer", tabName = "filter_explorer", icon = icon("th")),
            menuItem("Results", icon = icon("th"),
                menuSubItem("Filtered genes", tabName = "filter_results_genes"),
                menuSubItem("Filtered variants", tabName = "filter_results_variants"),
                menuSubItem("PanelApp and gene lists", tabName = "gene_lists")
            ),
            menuItem("Gene details", tabName = "gene_details", icon = icon("th")),
            menuItemOutput("exphunter_menu"),
            menuItem("Known variants", tabName= "known_variants", icon = icon("th")),
            menuItem("Saved variants", tabName= "preferred_vars", icon = icon("th")),
            menuItemOutput("coverage_menu")
            #menuItem("Cohort analysis", tabName= "cohort_analysis", icon = icon("th"))
        )
    ),
    dashboardBody(
      tags$head(tags$style(css)),
        tabItems(
            tabItem(tabName = "overview",
                    fluidRow(
                      column(2, h4("Case:")), column(2, verbatimTextOutput("activeCaseID")),
                      column(3, h4("Data release:")), column(5, verbatimTextOutput("releaseID"))),
                    fluidRow(
                        box(title = "Pedigree", width=6, status = "primary", solidHeader = TRUE,
                            withSpinner(plotOutput("ped"))
                        ),
                        box(title = "Variant types", width=6, status = "primary", solidHeader = TRUE,
                            withSpinner(plotOutput("Var_type_plot"))
                        )
                    ),
                    fluidRow(
                      box(title = "Phenotype information", width = 12, status = "primary", solidHeader = T, collapsible = T,
                          fluidRow(
                            column(4, h4("Disease:"), verbatimTextOutput("Disease")),
                            column(8, tableOutput("case_HPO_terms"))
                          )
                      )
                    ),
                    fluidRow(
                        box(title = "Variant consequences", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                            withSpinner(plotOutput("Var_consequence_plot"))
                        )
                    ),
                    fluidRow(
                        box(title = "Variant AF", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                            withSpinner(plotOutput("Var_PopAF_plot"))
                        )
                    ),
                    fluidRow(
                        box(title = "ROH distribution", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                            plotOutput("ROH_total_plot"),
                            uiOutput("ROH_plot_select"),
                            plotOutput("ROH_chrom_plot")
                        )
                    )
                    
            ),
            tabItem(tabName = "custom_genes_list_builder",
              DT::dataTableOutput("custom_genes_table"),
              fluidRow(column(12, align="center",
                actionButton("custom_genes_reset", "Reset list"),
                bsTooltip("custom_genes_reset", "Empty and reset the gene list"),
                actionButton("custom_genes_remove", "Remove selected")
              )),
              hr(),
              box(title = "Manual input", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                fluidRow(
                  column(6, 
                         textAreaInput("custom_genes_txt", "Gene list:", placeholder = "Official gene symbols one per line", rows = 10, resize = "vertical"),
                         verbatimTextOutput("custom_genes_list_length"),
                         actionButton("custom_genes_load_txt", "Add genes") 
                         ),
                  column(6, 
                         fileInput(inputId = "custom_genes_file",label = "Gene list file:", multiple=FALSE, accept="text/plain", placeholder = "gene list txt file"),
                         verbatimTextOutput("custom_genes_file_txt")
                        )
                )
              ),
              box(title = "Load from PanelApp", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                fluidRow(
                  column(3, selectInput("panelapp_confidence", "Min confidence level:", choices = c("green" = 3, "amber" = 2), selected = "green",multiple = FALSE)),
                  column(3, actionButton("panelapp_genes_load", "Add panelapp genes")),
                  column(3, actionButton("panelapp_reset", "Reset selection")),
                  column(3, verbatimTextOutput("panelapp_n_genes"))
                ),
                hr(),
                DT::dataTableOutput("panelapp_selection_table")
              ),
              box(title = "Load from ClinVar", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                  fluidRow(
                    column(3, actionButton("clinvar_genes_load", "Add clinvar genes")),
                    column(3, actionButton("clinvar_select_all", "Select all")),
                    column(3, actionButton("clinvar_reset", "Reset selection")),
                    column(3, verbatimTextOutput("clinvar_n_genes"))
                  ),
                  hr(),
                  DT::dataTableOutput("clinvar_selection_table")
              ),
              box(title = "Load from gene lists", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                  fluidRow(
                    column(4, actionButton("genelists_genes_load", "Add gene list genes")),
                    column(4, actionButton("genelists_reset", "Reset selection")),
                    column(4, verbatimTextOutput("genelists_n_genes"))
                  ),
                  hr(),
                  DT::dataTableOutput("genelists_selection_table")
              ),
              box(title = "Load from HPO terms", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                  fluidRow(
                    column(3, actionButton("hpo_genes_load", "Add HPO genes")),
                    column(3, actionButton("hpo_select_all", "Select all")),
                    column(3, actionButton("hpo_reset", "Reset selection")),
                    column(3, verbatimTextOutput("hpo_n_genes"))
                  ),
                  hr(),
                  DT::dataTableOutput("hpo_selection_table")
              ),
              box(title = "Load from HPO profile", width = 12, status = "primary", solidHeader = T, collapsed = T, collapsible = T,
                  fluidRow(
                    column(3, actionButton("hpo_profile_load", "Add genes HPO profile")),
                    column(3, actionButton("hpo_profile_reset", "Reset")),
                    column(3, verbatimTextOutput("hpo_profile_n_genes"))
                  ),
                  hr(),
                  fluidRow(
                    column(6,
                         textAreaInput("hpo_profile_txt", "HPO IDs:", placeholder = "HPO IDs one per line", rows = 10, resize = "vertical"),
                         uiOutput("hpo_n_min"),
                         column(6, actionButton("hpo_profile_set", "Load terms from patient profile")),
                         column(6, actionButton("hpo_profile_get", "Get genes for this profile"))
                         ),
                    column(6, DT::dataTableOutput("hpo_profile_table"))
                  )
              )
            ),
            tabItem(tabName = "filters_overview_tab", 
                    fluidRow(column(3,actionButton("set_filters","Configure filters"), align="center"),
                             column(3,actionButton("apply_filters","Apply filters"), align="center"),
                             column(3,downloadObjUI("get_json_filters", label = "Save filters settings"),
                                    br(),
                                    actionButton("latest_filter_button",label = "Load latest used filters"),
                                    align="center"),
                             column(3,
                                    fileInput(inputId = "load_filters",label = "Load filters:", multiple=FALSE, accept=".json", placeholder = "json file"),
                                    align="center")),
                    hr(),
                    h3("Variants filters"),
                    DT::dataTableOutput("vars_filters_table"),
                    h3("Segregation filters"),
                    DT::dataTableOutput("segregation_filters_table"),
                    h3("Gene scores filters"),
                    DT::dataTableOutput("genes_filters_table"),
                    h3("Genomic regions filters"),
                    DT::dataTableOutput("regions_filters_table")
            ),
            tabItem(tabName = "variants_filters_tab",
                    fluidRow(column(12, align="center", actionButton("next_variants_filters", "Save and next"))),
                    hr(),
                    uiOutput("vars_filters_UI")
            ),
            tabItem(tabName = "segregation_filters_tab", 
                    fluidRow(column(12, align="center",actionButton("next_segregation_filters", "Save and next"))),
                    hr(),
                    box(title = "Segregation filters", id = "segregation_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, width=12,
                        fluidRow(
                          column(6, textOutput("Total_affected")),
                          column(6, textOutput("Total_unaffected")) ),
                        uiOutput("segregation_controls"),
                        "Number of affected / unaffected carrying the variant with the given genotype",
                        "Homozygous and heterozygous conditions evaluated as AND, comphet as OR",
                        uiOutput("GQfilter_controls")
                    ) 
            ),
            tabItem(tabName = "genes_filters_tab",
                    fluidRow(column(12, align="center",actionButton("next_genes_filters", "Save and next"))),
                    hr(),
                    uiOutput("genes_filters_UI")
            ),
            tabItem(tabName = "regions_filters_tab",
                    fluidRow(column(12,align="center",actionButton("next_regions_filters", "Save and finish"))),
                    hr(),
                    box(title = "Genomic regions filter", id = "genomic_filters_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE, width = 12, 
                        h3("ROH regions"),
                        uiOutput("ROH_sample_select"),
                        uiOutput("ROH_filters_UI"),
                        h3("Custom BED regions"),
                        uiOutput("custom_bed_check"))
            ),
            tabItem(tabName = "filter_explorer",
                    uiOutput("filters_explorer_tab")
            ),
            tabItem(tabName = "filter_results_genes",
                    uiOutput("genes_results")
            ),
            tabItem(tabName = "filter_results_variants",
                    uiOutput("variants_results")
            ),
            tabItem(tabName = "gene_lists",
                    uiOutput("gene_lists_results")
            ),
            tabItem(tabName = "gene_details",
                    uiOutput("gene_detail_ui")
            ),
            tabItem(tabName = "expansion_hunter",
                    uiOutput("exphunter_selection"),
                    DT::dataTableOutput("exphunter_variants"),
                    plotOutput("exphunter_plot", width = "100%", height = "800px")
            ),
            tabItem(tabName = "known_variants", 
                    box(title = "ClinVar known variants", id = "known_clinvar_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, width=12,
                        fluidRow(
                          column(6, h4("ClinVar significance for selected variant:", verbatimTextOutput("clinvar_var_sig"))),
                          column(6, h4("ClinVar phenotype for selected variant:", verbatimTextOutput("clinvar_var_pheno")))
                        ),
                        br(),
                        DT::dataTableOutput("known_clinvar_tab")
                    ),
                    box(title = "COSMIC known variants", id = "known_cosmic_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = TRUE, width=12,
                        DT::dataTableOutput("known_cosmic_tab")
                    ) 
            ),
            tabItem(tabName = "preferred_vars",
                    fluidRow(column(3, actionButton("preferred_vars_remove", label = "Remove selected variants")),
                             column(3, actionButton("preferred_vars_reset", label = "Remove all variants")),
                             column(6, fileInput(inputId = "preferred_vars_file",label = "Preferred variants:", multiple=FALSE, accept=".csv", placeholder = "preferred vars csv file"), align="center")
                    ),
                    br(),
                    DT::dataTableOutput("preferred_vars_tab")
            ),
            tabItem(tabName = "coverage_explorer",
                    fluidRow(
                        column(2,selectInput("chr_coverage", "Chromosome: ", choices = paste("chr", c(1:22,"X","Y","M"), sep=""))),
                        column(3,uiOutput("gene_cov_select")),
                        column(2,textInput("cov_region_pad", label = "Flanking region (kb): ", value = "500")),
                        column(3, sliderInput("smooth_factor",
                                              label = "smooth", 
                                              min = 1, max = 20, step = 1, value = 1)),
                        column(2, textOutput("smooth_dimension")) ),
                    fluidRow(column(12, align="center",actionButton("plot_coverage","Plot coverage"))),
                    withSpinner(plotlyOutput("coverage_plot"))
            ),
            tabItem(tabName = "cohort_analysis",
                    fluidRow(
                      column(6, 
                             selectInput(inputId = "cohort_samples", label = "samples to analyze", choices = samplesID, size = 10, selectize=F, multiple=T)
                             ),
                      column(4,
                             fluidRow(actionButton("apply_cohort", "Cohort analysis"), align="center"),
                             br(),
                             textInput("cohort_threads","N threads: ", value=4, placeholder="N parallel threads"),
                             br(),
                             fluidRow(downloadObjUI("save_cohort", label = "Save cohort results"), align="center"),
                             br(),
                             fluidRow(verbatimTextOutput("cohort_exit_status"))
                             )
                    ),
                    br(),
                    box(title = "Cohort genes", id = "cohort_genes_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = TRUE, width=12,
                        DT::dataTableOutput("cohort_genes_df")),
                    box(title = "Cohort variants", id = "cohort_vars_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = TRUE, width=12,
                        DT::dataTableOutput("cohort_vars_df")),
                    box(title = "Cohort comphet", id = "cohort_comphet_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = TRUE, width=12,
                        DT::dataTableOutput("cohort_comphet_df"))
                    
                    )
        )
    )
)

## SERVER -------------------------------
server <- function(input, output, session) {
  ## Notifications, tasks and messages -------------------- 
  output$NotificationMenu <- renderMenu({
    dropdownMenu(type = "notifications", .list = RV$notifications)
  })
  
  output$MessageMenu <- renderMenu({
    dropdownMenu(type = "messages", .list = RV$messages)
  })
  
  ## Side bar reactive ----------------------------
  #Use in UI: menuItemOutput("recOpt") 
  output$exphunter_menu <- renderMenu({
    req(inherits(RV$data, "list"))
    if(!is.null(RV$data$ExpHunter))
      menuItem("Expansion Hunter", tabName= "expansion_hunter", icon = icon("th"))
  })
  
  output$coverage_menu <- renderMenu({
    if( !is.null(Coverage_dir) & length(list.files(Coverage_dir, pattern = ".bed.gz")) > 0)
      menuItem("Explore coverage", tabName= "coverage_explorer", icon = icon("th"))
  })
  
  ## Load Files --------------------
  observeEvent(input$decrypt_button, {
      #Remember that df in pre-processed object are generated with fread so slicing works differently [,..indexes]
      updateTabItems(session, "tabs", "overview")
    
      #Reset relevant values when a new sample is loaded
      #RV$messages <- list()
      #RV$custom_genes <- data.frame(gene=character(), source=character(), stringsAsFactors = F)
      #RV$customBed_ranges <- FALSE
      RV$notifications <- list()
      RV$cov_file <- NULL
      RV$custom_genes_txt_length <- 0
      RV$vars_pass_GQ <- NULL
      RV$vars_pass_BED <- NULL
      RV$vars_pass_ROH <- NULL 
      RV$vars_pass_segregation <- NULL
      RV$vars_pass_filters <- NULL
      RV$vars_filters_df <- NULL
      RV$segregation_filters_df <- NULL
      RV$genes_filters_df <- NULL
      RV$regions_filters_df <- NULL
      RV$filters_json <- NULL
      RV$filtered_applied <- FALSE
      
      showModal(
        modalDialog(paste0("loading ", input$CaseCode, "...\nPlease wait"), footer = NULL)
      )
      #message("loading ", paste0(data_dir,"/",input$CaseCode,".RData[.enc]"))
      #Load data from plain or encrypted object
      if (file_test("-f", paste0(data_dir,"/",input$CaseCode,".RData.enc"))) {
          RV$data <- decrypt_datafile(paste0(data_dir,"/",input$CaseCode,".RData.enc"), pwd = input$pwd)
      } else if (file_test("-f", paste0(data_dir,"/",input$CaseCode,".RData"))) {
          RV$data <- readRDS(paste0(data_dir,"/",input$CaseCode,".RData"))
      }
      if (inherits(RV$data, "list")) {
          RV$data$variants_df <- as.data.frame(RV$data$variants_df %>% replace_na(app_settings$fill_na$fill_na_vars))
          RV$data$variants_df$Class <- "PASS"
          
          RV$data$comphet_df <- as.data.frame(RV$data$comphet_df)
          RV$data$comphet_df$Class <- "PASS"
          
          RV$data$segregation_df <- as.data.frame(RV$data$segregation_df)
          RV$data$segregation_df$sup_dnm[RV$data$segregation_df$sup_dnm < 0] <- 0
          
          genes_scores_valid <- RV$data$genes_scores %>% 
            filter(!is.na(gado_zscore)) %>%
            arrange(dplyr::desc(gado_zscore)) %>%
            mutate(gado_perc=1-(row_number()/nrow(.)))
          genes_scores_na <- RV$data$genes_scores %>% 
            filter(is.na(gado_zscore)) %>%
            mutate(gado_perc=NA)
          genes_scores_perc <- rbind(genes_scores_valid,genes_scores_na)
          
          RV$data$genes_scores <- genes_scores_perc %>% relocate(gado_perc, .after="gado_zscore")
          RV$data$genes_scores <- as.data.frame(RV$data$genes_scores %>% replace_na(app_settings$fill_na$fill_na_genes))
          RV$data$genes_scores$Class <- "PASS"
          
          #Compute cohort normalized GADO is slow, so yo can turn it on or off from app config
          if (compute_cohort_GADO) {
            RV$data$genes_scores$cohort_norm_Z <- apply(RV$data$genes_scores,1,function(x) computeNormZ(x["gene"],x["gado_zscore"]))
          }
          
          tryCatch({
            RV$data$ROH_data$ROHClass <- cut(RV$data$ROH_data$Length_bp, 
                                             breaks = c(0,500000,2000000,max(RV$data$ROH_data$Length_bp)), 
                                             labels = c("small (< 500kb)","medium (500kb-2Mb)","large (>= 2Mb)"))
          },
          error=function(cond) {
            RV$data$ROH_data <- NULL
            RV$notifications[["ROH"]] <- notificationItem(
              text = "Error loading ROH dimension, ROH plot will not be available in overview",
              icon = icon("exclamation-circle"),
              status = "warning")
          })
         
          RV$notifications[["decrypt"]] <- notificationItem(
              text = paste0("Loaded data for ", input$CaseCode),
              icon = icon("check-circle"),
              status = "success")
          RV$notifications[["pedigree"]] <- notificationItem(
              text = paste0("Individuals: ", RV$data$n_all_samples),
              icon = icon("check-circle"),
              status = "success")
          RV$notifications[["variants"]] <- notificationItem(
              text = paste0("variants loaded: ", length(base::unique(RV$data$variants_df$var_id))),
              icon = icon("check-circle"),
              status = "success")
          RV$notifications[["genes"]] <- notificationItem(
              text = paste0("distinct genes: ", length(base::unique(RV$data$genes_scores$gene))),
              icon = icon("check-circle"),
              status = "success")
          
          RV$GQ_cols_all <- which(colnames(RV$data$variants_df) %in% paste("GQ", RV$data$all_samples, sep="_"))
          RV$GQ_cols_affected <- which(colnames(RV$data$variants_df) %in% paste("GQ", RV$data$affected_samples, sep="_"))
          RV$maxGQ <- max(RV$data$variants_df[,RV$GQ_cols_all], na.rm = T)
          
      } else {
          RV$notifications[["decrypt"]] <- notificationItem(
              text = "Error loading data! Missing file or wrong password!",
              icon = icon("exclamation-circle"),
              status = "danger")
      }
      removeModal()
  })
  
  observeEvent(input$custom_bed, {
      req(input$custom_bed)
      tryCatch({
      bed_df <- read.table(input$custom_bed$datapath, sep="\t", header=F, stringsAsFactors = F)
  
      #ID is set to ROH_nline for compatibility with the bed filtering module
      RV$customBed_ranges <- GRanges(seqnames = bed_df$V1, ranges = IRanges(bed_df$V2, end=bed_df$V3), ID = paste("ROH",1:nrow(bed_df),sep="_"))
      
      RV$notifications[["custom_bed"]] <- notificationItem(
          text = paste0(nrow(bed_df), " regions loaded from BED"),
          icon = icon("check-circle"),
          status = "success")
      } , error=function(cond) {
          RV$notifications[["custom_bed"]] <- notificationItem(
              text = paste0("Failed loading custom BED file"),
              icon = icon("exclamation-circle"),
              status = "danger")
      })
  })
    
  ## Data reactive tables configuration -------------------------
  
  filters_summ_genes <- reactive ({
      tot_genes <- nrow(RV$data$genes_scores)
      PASS_count <- length(RV$genes_pass_filters)
      filters_summ_genes <- data.frame(
          Filter=c("Genes filters"), 
          PASS=PASS_count/tot_genes, 
          FILTERED=(tot_genes-PASS_count)/tot_genes)
      filters_summ_genes <- gather(filters_summ_genes, key="Class", value="Count", PASS:FILTERED)
      filters_summ_genes$Count <- as.numeric(filters_summ_genes$Count)
      filters_summ_genes
  })
  
  filters_summ_vars <- reactive ({
      tot_vars <- RV$data$variants_df %>% select(var_id) %>% distinct() %>% nrow()
      comphet_seg_vars <- RV$data$comphet_df %>% filter(
        rec_id %in% RV$vars_pass_segregation) %>% 
        gather(key="Variant",value="VarID",v1:v2) %>% select(VarID)
      PASS_counts <- c(
          RV$data$variants_df %>% filter(rec_id %in% RV$vars_pass_filters$vars) %>% select(var_id) %>% distinct() %>% nrow(),
          RV$data$variants_df %>% filter(rec_id %in% RV$vars_pass_GQ) %>% select(var_id) %>% distinct() %>% nrow(),
          RV$data$variants_df %>% filter(rec_id %in% c(RV$vars_pass_segregation, comphet_seg_vars$VarID)) %>% select(var_id) %>% distinct() %>% nrow(),
          RV$data$variants_df %>% filter(rec_id %in% RV$vars_pass_ROH) %>% select(var_id) %>% distinct() %>% nrow(),
          RV$data$variants_df %>% filter(rec_id %in% RV$vars_pass_BED) %>% select(var_id) %>% distinct() %>% nrow()
      )
      filters_summ_vars <- data.frame(Filter=c("variants","GQ","segregation","ROH","custom BED"),
          PASS=PASS_counts/tot_vars, 
          FILTERED=(tot_vars-PASS_counts)/tot_vars)
      
      filters_summ_vars <- gather(filters_summ_vars, key="Class", value="Count", PASS:FILTERED)
      filters_summ_vars$Count <- as.numeric(filters_summ_vars$Count)
      filters_summ_vars
  })
  
  PanelApp_panels_df <- reactive ({
      panels_idx <- PanelApp_genes$panel_idx[PanelApp_genes$entity_name %in% RV$filtered_genes_list]
      as.data.frame(PanelApp_data %>% filter(id %in% panels_idx))
  })
  
  geneLists_df <- reactive ({
      geneLists_idx <- geneLists_genes$genelist_idx[geneLists_genes$entity_name %in% RV$filtered_genes_list]
      as.data.frame(geneLists_data %>% filter(id %in% geneLists_idx))
  })
  
  ClinVar_df <- reactive ({
      as.data.frame(ClinVar_genes %>% filter(gene %in% RV$filtered_genes_list))
  }) 
  
  ## Overview tab ---------------------------
  
  output$releaseID <- renderText({
    req(inherits(RV$data, "list"))
    if (is.null(RV$data$releaseID)) {
      return("No release ID provided")
    } else {
      return(RV$data$releaseID)
    }
  })
  
  output$activeCaseID <- renderText({
    req(inherits(RV$data, "list"))
    return(RV$data$pedigree)
  })
  
  output$Total_affected <- renderText({
      paste0("\tTotal number of affected individuals in this pedigree: ", RV$data$n_affected)
  })
  
  output$Total_unaffected <- renderText({
      paste0("\tTotal number of unaffected individuals in this pedigree: ", RV$data$n_unaffected)
  })
  
  output$ped <- renderPlot({
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded"))
      shiny::validate(need(!is.null(RV$data$ped), "No PED data found for this sample"),
                      need(RV$data$n_all_samples > 1, "SINGLETON"))
      plot.pedigree(RV$data$ped, mar = c(5, 3, 5, 3))
  })
  
  output$Var_consequence_plot <- renderPlot({
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded"))
      var_data <- RV$data$variants_df %>% select(var_id, consequence) %>% distinct()
      ggplot(var_data, aes(x=consequence)) + geom_bar()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_sqrt()
  })
  
  output$Disease <- renderText({
    shiny::validate(need(inherits(RV$data, "list"), "Please load a case"),
                    need(cohort_HPO, "No phenotype data present"))
    disease <- base::unique(cohort_HPO$Disease[cohort_HPO$CaseID == RV$data$pedigree])
  })
  
  output$case_HPO_terms <- renderTable(align = "c", rownames = F, striped = T, {
    shiny::validate(need(inherits(RV$data, "list"), "Please load a case"),
                    need(cohort_HPO, "No phenotype data present"))
    HPO_ids <- base::unique(cohort_HPO$HPO[cohort_HPO$CaseID == RV$data$pedigree])
    HPO_table <- HPO_obo[HPO_obo$HPO_id %in% HPO_ids,]  
  })
  
  output$Var_type_plot <- renderPlot({
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded"))
      ggplot(RV$data$variants_df, aes(x=var_type)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_sqrt()
  })
  
  output$Var_PopAF_plot <- renderPlot({
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded"))
      AF_data <- RV$data$variants_df %>% select(var_id, max_pop_af, cohort_af, var_type) %>% distinct()
      AF_data$PopAF_level <- cut(AF_data$max_pop_af, breaks=c(0,0.005,0.01,0.05,1), labels = c("VERY_RARE\n(<= 0.005)","RARE\n(<= 0.01)","LOW_FREQ\n(<= 0.05)", "COMMON"), include.lowest = T)
      AF_data$CohortAF_level <- cut(AF_data$cohort_af, breaks=c(0,0.005,0.01,0.05,1), labels = c("VERY_RARE\n(<= 0.005)","RARE\n(<= 0.01)","LOW_FREQ\n(<= 0.05)", "COMMON"), include.lowest = T)
      
      cohortAF_p <- ggplot(AF_data, aes(x=CohortAF_level, fill=CohortAF_level)) + geom_bar() + scale_y_sqrt() + scale_fill_brewer(palette="Set1") + labs(x="", title = "Cohort AF distribution") + theme(legend.position = "none")
      PopAF_p <- ggplot(AF_data, aes(x=PopAF_level, fill=PopAF_level)) + geom_bar() + scale_y_sqrt() + scale_fill_brewer(palette="Set1") + labs(x="", title = "Max pop AF distribution") + theme(legend.position = "none")
      AF_scatter <- ggplot(AF_data, aes(x=max_pop_af, y=cohort_af, color=var_type)) + geom_point() + scale_y_sqrt() + scale_fill_brewer(palette="Set1") + labs(title = "AF comparison") + theme(legend.position = "none")
      plot_layout <- rbind(c(1,1,2,2),c(NA,3,3,NA))
      grid.arrange(cohortAF_p, PopAF_p, AF_scatter, layout_matrix = plot_layout)
      
  })
  
  output$ROH_total_plot <- renderPlot({
      shiny::validate(need(RV$data$ROH_data, "No ROH data found for this sample"))
      Total_ROH <- aggregate(RV$data$ROH_data$Length_bp, list(RV$data$ROH_data$Sample, RV$data$ROH_data$ROHClass), FUN=sum)
      Total_ROH$pctROH <- Total_ROH$x/genome_size
      ggplot(Total_ROH, aes(x=Group.1,y=pctROH,fill=Group.2)) + 
          geom_bar(stat="identity", position = position_dodge(0.9)) + 
          geom_hline(yintercept = 0.1, linetype="dashed") + 
          labs(x="sample", y="fraction of genome within ROH", fill="ROH size") + 
          format1
  })
  
  output$ROH_plot_select <- renderUI ({
      selectInput("ROH_plot_sample", label = "Select sample:", choices = RV$data$all_samples, multiple = FALSE, selected = RV$data$affected_samples[1])
  })
  
  output$ROH_chrom_plot <- renderPlot({
    shiny::validate(need(RV$data$ROH_data, "No ROH data found for this sample"))
      Total_ROH_bychr <- aggregate(RV$data$ROH_data$Length_bp, list(RV$data$ROH_data$Sample, RV$data$ROH_data$ROHClass, RV$data$ROH_data$Chromosome), FUN=sum)
      Total_ROH_bychr <- merge(Total_ROH_bychr, chr_sizes, by.x="Group.3", by.y="V1")
      Total_ROH_bychr$pctROH <- Total_ROH_bychr$x / Total_ROH_bychr$V2
      ggplot(Total_ROH_bychr[Total_ROH_bychr$Group.1 == input$ROH_plot_sample,], aes(x=Group.3,y=pctROH,fill=Group.2)) + 
          geom_bar(stat="identity", position = position_dodge(0.9)) + 
          geom_hline(yintercept = 0.1, linetype="dashed") + 
          labs(x="chromosome", y="fraction of chromosome within ROH", fill="ROH size", title=paste0("ROH distribution by chromosome for ", input$ROH_plot_sample)) + 
          format1
  })
  
  ## Custom genes tab ------------------------------
  ### MANUAL INPUT ###
  observeEvent(input$custom_genes_reset, {
    RV$custom_genes <- data.frame(gene=character(), source=character(), stringsAsFactors = F)
  })
  
  observeEvent(input$custom_genes_remove, {
    shiny::validate(need(input$custom_genes_table_rows_selected, "Select a gene"))
    RV$custom_genes <- RV$custom_genes[-input$custom_genes_table_rows_selected,]
  })
  
  output$custom_genes_table <- DT::renderDataTable(selection="multiple", {
    RV$custom_genes
  })
  
  observeEvent(input$custom_genes_txt, {
    genes <- unlist(strsplit(input$custom_genes_txt, "\n"))
    RV$custom_genes_txt_length <- length(genes) 
  })
  
  output$custom_genes_list_length <- renderText({
    paste0(RV$custom_genes_txt_length, " genes in the list")
  })
  
  observeEvent(input$custom_genes_load_txt, {
    genes <- unlist(strsplit(input$custom_genes_txt, "\n"))
    mydf <- data.frame(gene=genes, source=rep("manual_input", length(genes)))
    RV$custom_genes <- rbind(RV$custom_genes, mydf) %>% distinct()
  })
  
  observeEvent(input$custom_genes_file, {
    req(input$custom_genes_file)
    tryCatch({
      genes <- scan(input$custom_genes_file$datapath,what="",sep="\n")
      RV$custom_genes_n_loaded <- length(genes)
      values <- paste(c(input$custom_genes_txt,
                        paste(genes, collapse="\n") ), collapse="\n")
      updateTextAreaInput(session, "custom_genes_txt", value = values) 
    }, error=function(cond) {
      RV$notifications[["custom_file"]] <- notificationItem(
        text = paste0("Failed loading custom genes list"),
        icon = icon("exclamation-circle"),
        status = "danger")
    })
    
  })
  
  output$custom_genes_file_txt <- renderText({
    paste0(RV$custom_genes_n_loaded, " genes loaded from file")
  })
  
  ### PANELAPP SELECTION ###
  output$panelapp_selection_table <- DT::renderDataTable(selection="multiple", {
    PanelApp_data
  })
  
  panelapp_selected_genes <- reactive ({
    shiny::req(input$panelapp_selection_table_rows_selected)
    panelID <- PanelApp_data[input$panelapp_selection_table_rows_selected, "id"]
    genes_df <- PanelApp_genes %>% 
      filter(panel_idx %in% panelID, confidence_level >= input$panelapp_confidence) %>%  
      left_join(., PanelApp_data[,c("id","name")], by = c("panel_idx" = "id")) %>%
      select(entity_name, name)
    colnames(genes_df) <- c("gene","source")
    genes_df <- genes_df[order(genes_df$source, genes_df$gene),]
  })
  
  output$panelapp_n_genes <- renderText({
    paste0(nrow(panelapp_selected_genes()), " genes to be imported")
  })
  
  panelapp_selection_proxy <- DT::dataTableProxy("panelapp_selection_table")
  
  observeEvent(input$panelapp_reset, {
    panelapp_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$panelapp_genes_load, {
    shiny::req(nrow(panelapp_selected_genes()) > 0)
    RV$custom_genes <- rbind(RV$custom_genes, panelapp_selected_genes()) %>% distinct()
    panelapp_selection_proxy %>% selectRows(NULL)
  })
  
  ### CLINVAR SELECTION ###
  output$clinvar_selection_table <-  DT::renderDataTable(selection="multiple", {
    ClinVar_genes 
  })
    
  clinvar_selection_proxy <- DT::dataTableProxy("clinvar_selection_table")
  
  output$clinvar_n_genes <- renderText({
    paste0(nrow(clinvar_selected_genes()), " genes to be imported")
  })
  
  observeEvent(input$clinvar_reset, {
    clinvar_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$clinvar_genes_load, {
    shiny::req(nrow(clinvar_selected_genes()) > 0)
    RV$custom_genes <- rbind(RV$custom_genes, clinvar_selected_genes()) %>% distinct()
    clinvar_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$clinvar_select_all, {
    clinvar_selection_proxy %>% selectRows(input$clinvar_selection_table_rows_all)
  })
  
  clinvar_selected_genes <- reactive ({
    shiny::req(input$clinvar_selection_table_rows_selected)
    genes <- ClinVar_genes[input$clinvar_selection_table_rows_selected, "gene"]
    genes_df <- data.frame(gene=genes, source=rep("ClinVar", length(genes)))
    genes_df <- genes_df[order(genes_df$gene),]
  })
  
  ### GENE LISTS SELECTION ###
  output$genelists_selection_table <- DT::renderDataTable(selection="multiple", {
    geneLists_data
  })
  
  genelists_selection_proxy <- DT::dataTableProxy("genelists_selection_table")
  
  genelists_selected_genes <- reactive ({
    shiny::req(input$genelists_selection_table_rows_selected)
    panelID <- geneLists_data[input$genelists_selection_table_rows_selected, "id"]
    genes_df <- geneLists_genes %>% 
      filter(genelist_idx %in% panelID) %>%  
      left_join(., geneLists_data[,c("id","list_name")], by = c("genelist_idx" = "id")) %>%
      select(entity_name, list_name)
    colnames(genes_df) <- c("gene","source")
    genes_df <- genes_df[order(genes_df$source, genes_df$gene),]
  })
  
  output$genelists_n_genes <- renderText({
    paste0(nrow(genelists_selected_genes()), " genes to be imported")
  })
  
  observeEvent(input$genelists_reset, {
    genelists_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$genelists_genes_load, {
    shiny::req(nrow(genelists_selected_genes()) > 0)
    RV$custom_genes <- rbind(RV$custom_genes, genelists_selected_genes()) %>% distinct()
    genelists_selection_proxy %>% selectRows(NULL)
  })
  
  ### HPO GENES SELECTION ###
  output$hpo_selection_table <-  DT::renderDataTable(selection="multiple", {
    HPO_genes 
  })
  
  hpo_selection_proxy <- DT::dataTableProxy("hpo_selection_table")
  
  output$hpo_n_genes <- renderText({
    paste0(nrow(hpo_selected_genes()), " genes to be imported")
  })
  
  observeEvent(input$hpo_reset, {
    hpo_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$hpo_genes_load, {
    shiny::req(nrow(hpo_selected_genes()) > 0)
    RV$custom_genes <- rbind(RV$custom_genes, hpo_selected_genes()) %>% distinct()
    hpo_selection_proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$hpo_select_all, {
    hpo_selection_proxy %>% selectRows(input$hpo_selection_table_rows_all)
  })
  
  hpo_selected_genes <- reactive ({
    shiny::req(input$hpo_selection_table_rows_selected)
    genes <- HPO_genes[input$hpo_selection_table_rows_selected, "gene"]
    hpo_ids <- HPO_genes[input$hpo_selection_table_rows_selected, "HPO_id"]
    genes_df <- data.frame(gene=genes, source=hpo_ids)
    genes_df <- genes_df[order(genes_df$source, genes_df$gene),]
  })
  
  ### HPO PROFILE SELECTION ###
  observeEvent(input$hpo_profile_load,{
    shiny::req(nrow(hpo_profile_genes()) > 0)
    hpo_profile_df <- data.frame(gene=hpo_profile_genes()$gene, source="HPO_profile", stringsAsFactors = F)
    RV$custom_genes <- rbind(RV$custom_genes, hpo_profile_df) %>% distinct()
  })
  
  observeEvent(input$hpo_profile_reset, {
    updateTextAreaInput(session, "hpo_profile_txt", value = "")
  })
  
  output$hpo_profile_n_genes <- renderText({
    paste0(nrow(hpo_profile_genes()), " genes to be imported")  
  })
  
  observeEvent(input$hpo_profile_set, {
    req(inherits(RV$data, "list"), cohort_HPO)
    HPO_ids <- base::unique(cohort_HPO$HPO[cohort_HPO$CaseID == RV$data$pedigree])
    updateTextAreaInput(session, "hpo_profile_txt", value = paste(HPO_ids, collapse="\n"))  
  })
  
  observeEvent(input$hpo_profile_get, {
    genes <- HPO_genes[input$hpo_selection_table_rows_selected, "gene"]
    hpo_ids <- HPO_genes[input$hpo_selection_table_rows_selected, "HPO_id"] 
  })
  
  output$hpo_n_min <- renderUI({
    n_hpos <- length(hpo_terms_in_profile())
    selectInput("hpo_n_min_select", "Select genes associated to at least N HPOs:", choices = seq(0,n_hpos), selected = 1, multiple = FALSE)  
  }) 
  
  hpo_terms_in_profile <- reactive({
    unlist(strsplit(input$hpo_profile_txt, "\n"))  
  })
  
  hpo_profile_genes <- reactive({
    HPO_genes %>% 
      filter(HPO_id %in% hpo_terms_in_profile()) %>% 
      group_by(gene) %>% 
      mutate(N_HPOs = n()) %>% 
      select(gene, N_HPOs) %>%
      filter(N_HPOs >= input$hpo_n_min_select) %>%
      distinct()
  })
  
  output$hpo_profile_table <- DT::renderDataTable(selection="none", {
    req(nrow(hpo_profile_genes()) > 0)
    hpo_profile_genes()
  })
  
  
  ## Filters ---------------------------------
    ## Filters Overview ----------------------------
    #Apply the filters and go to gene result page
    observeEvent(input$apply_filters, {
      shiny::validate(need(inherits(RV$data, "list"), FALSE))
      req(RV$vars_filters_df, RV$segregation_filters_df, RV$genes_filters_df, RV$regions_filters_df)
      
      showModal(
        modalDialog("Filtering variants... Please wait", footer = NULL)
      )
      
      #VARIANTS FILTER
      RV$vars_pass_filters <- callModule(getPASSVars_filters, "variants_filters", filters_settings$VARIANTS, RV$data$variants_df, RV$data$comphet_df)
      message("PASS VARS FILTERS ", length(RV$vars_pass_filters$vars))
      
      #GENE SCORES FILTER
      RV$genes_pass_filters <- callModule(getPASSGenes_filters, "genes_filters", filters_settings$GENES, RV$data$genes_scores)
      message("PASS GENES FILTERS ", length(RV$genes_pass_filters))
      
      #SEGREGATION FILTER
      RV$vars_pass_segregation <- callModule(segregationModule, "segregation", segregation_df = RV$data$segregation_df, cols_names = segregation_cols)
      message("PASS SEGREGATION ", length(RV$vars_pass_segregation))
      
      #ROH FILTER
      #call the bed regions module and get rec_id for variants in the ROH regions
      if (!is.null(RV$data$ROH_ranges)) {
        #message("call bed module")
        RV$vars_pass_ROH <- callModule(bedfilterModule,"ROH_filter",variants_ranges=RV$data$variants_ranges, bed_ranges=RV$data$ROH_ranges, multiplier=1000)
      } else {
        RV$vars_pass_ROH <- RV$data$variants_df$rec_id
      }
      message("PASS ROH ", length(RV$vars_pass_ROH))
      
      #BED REGIONS FILTER
      #call the bed regions module and get rec_id for variants in the custom BED regions
      if (inherits(RV$customBed_ranges, "GRanges")) {
        RV$vars_pass_BED <- callModule(bedfilterModule,"bed_filter",variants_ranges=RV$data$variants_ranges, bed_ranges=RV$customBed_ranges)
      } else {
        RV$vars_pass_BED <- RV$data$variants_df$rec_id
      }
      message("PASS BED ", length(RV$vars_pass_BED))
      
      #Get rec_id for variants passing the GQ filter
      RV$vars_pass_GQ <- callModule(GQfilterModule, "GQ_filter", 
                                 variants_df = RV$data$variants_df, 
                                 GQ_cols = RV$GQ_cols_all, 
                                 affected_cols = RV$GQ_cols_affected,
                                 exclude_var_type = sv_vars)
      message("PASS GQ ", length(RV$vars_pass_GQ))
      
      req(RV$vars_pass_GQ, RV$vars_pass_BED, RV$vars_pass_ROH, RV$vars_pass_segregation, RV$vars_pass_filters)
      RV$effective_filters <- rbind(RV$vars_filters_df, RV$segregation_filters_df, RV$genes_filters_df, RV$regions_filters_df) #TABLES OF FILTERS
      
      # Update tables based on filters
      ## VARIANTS
      #Get the final list of vars in accepted comphet intersecting filters and segregation
      comphet_single_vars <- RV$data$comphet_df %>% filter(
        rec_id %in% intersect(RV$vars_pass_segregation,RV$vars_pass_filters$comphet)) %>% 
        gather(key="Variant",value="VarID",v1:v2) %>% select(VarID)
      
      #Get final list of accepted variants
      #the final list of segregating accepted vars is equal to:
      #single vars passing filters and segregation
      #vars part of a comphet passing all filters + segregation
      accepted_vars_list <- base::unique(intersect(intersect(intersect(intersect(
        RV$vars_pass_ROH, RV$vars_pass_GQ), 
        RV$vars_pass_BED),
        RV$vars_pass_filters$vars),
        RV$vars_pass_segregation))
      
      accepted_vars_list <- base::unique(c(accepted_vars_list, comphet_single_vars$VarID))
      
      #Update Class column (PASS/FILTER)
      RV$data$variants_df <- RV$data$variants_df %>% 
        mutate(Class = ifelse(rec_id %in% accepted_vars_list, "PASS","FILTER"))
      
      ## COMPHET
      RV$filtered_vars_list <- RV$data$variants_df$rec_id[RV$data$variants_df$Class == "PASS"]
      RV$data$comphet_df <- RV$data$comphet_df %>% 
        mutate(Class = ifelse(
          v1 %in% RV$filtered_vars_list & 
          v2 %in% RV$filtered_vars_list &
          rec_id %in% intersect(RV$vars_pass_segregation,RV$vars_pass_filters$comphet), "PASS", "FILTER"))
      
      ## GENES
      RV$filtered_vars_list <- base::unique(c(
        RV$filtered_vars_list, 
        RV$data$comphet_df$rec_id[RV$data$comphet_df$Class == "PASS"])) 
    
      RV$filtered_genes_list <- RV$data$genes_df %>% 
        filter(variants %in% RV$filtered_vars_list, gene %in% RV$genes_pass_filters) %>%
        pull(gene) %>% base::unique()
      
      RV$data$genes_scores <- RV$data$genes_scores %>% 
        mutate(Class = ifelse(gene %in% RV$filtered_genes_list, "PASS", "FILTER"))
      
      RV$filtered_applied <- TRUE
      write_json(RV$filters_json, path=paste0(user_dir, "/Last_filters.json"), pretty=T, auto_unbox=T)
      
      removeModal()
      updateTabItems(session, "tabs", "filter_results_genes")
    }) 
    
    #Load latest filters if they are present
    observeEvent(input$latest_filter_button,{
      latest_file <- paste0(user_dir,"/Last_filters.json")
      
      if (file.exists(latest_file)) {
        tryCatch({
          filters_json <- read_json(latest_file)
          #load settings from the json file for each filter module
          callModule(loadSettings_filters, "variants_filters", filters_settings$VARIANTS, filters_json$VARIANTS)
          callModule(loadSettings_segregation, "segregation", filters_json$SEGREGATION)
          callModule(loadSettings_GQ, "GQ_filter", filters_json$GQ)
          callModule(loadSettings_filters, "genes_filters", filters_settings$GENES, filters_json$GENES)
          callModule(loadSettings_regions, "ROH_filter", filters_json$ROH)
          callModule(loadSettings_regions, "bed_filter", filters_json$BED)
  
          #update filter settings summary data frames
          RV$vars_filters_df <- callModule(getDF_filters, "variants_filters", filters_settings$VARIANTS)
          df1 <- callModule(getDF_segregation, "segregation")
          df2 <- callModule(getDF_GQ, "GQ_filter")
          RV$segregation_filters_df <- rbind(df1, df2)
          RV$genes_filters_df <- callModule(getDF_filters, "genes_filters", filters_settings$GENES)
          df1 <- callModule(getDF_regions, "bed_filter", "custom BED")
          df2 <- callModule(getDF_regions, "ROH_filter", "ROH")
          RV$regions_filters_df <- rbind(df1, df2) 
          
          #Create a success notification
          RV$notifications[["filters_file"]] <- notificationItem(
            text = paste0("Filters configuration loaded correctly"),
            icon = icon("check-circle"),
            status = "success")
          
          #Move to variants tab for filter inspection
          updateTabItems(session, "tabs", "variants_filters_tab")
        }, error=function(cond) {
          RV$notifications[["filters_file"]] <- notificationItem(
            text = paste0("Failed loading filters configuration"),
            icon = icon("exclamation-circle"),
            status = "danger")
        })
      
        message_text <- paste0("Load latest filters (", strftime(file.info(latest_file)$ctime,format = "%Y-%m-%d %H:%M:%S"),")")
        showModal(
          modalDialog(message_text)
        )
      } else {
        showModal(
          modalDialog("No latest filters file found in user data folder")
        )
      }
    })
    
    #Move to variant filters tab
    observeEvent(input$set_filters, {
      updateTabItems(session, "tabs", "variants_filters_tab")
    })
    
    #Constantly update a json of filters settings and prepare object to save
    toListenFilters <- reactive({
      list(RV$vars_filters_df, RV$segregation_filters_df, RV$genes_filters_df, RV$regions_filters_df)
    })
    
    observeEvent( toListenFilters(), {
      req(inherits(RV$data, "list"))
      variants_json <- callModule(getJSON_filters, "variants_filters", filters_settings$VARIANTS)
      segregation_json <- callModule(getJSON_segregation, "segregation")
      GQ_json <- callModule(getJSON_GQ, "GQ_filter")
      genes_json <- callModule(getJSON_filters, "genes_filters", filters_settings$GENES)
      ROH_json <- callModule(getJSON_regions, "ROH_filter")
      BED_json <- callModule(getJSON_regions, "bed_filter")
      
      RV$filters_json <- list(
        APP_VERSION = APP_VERSION,
        DATA_VERSION = RV$data$releaseID,
        VARIANTS = variants_json,
        SEGREGATION = segregation_json,
        GQ = GQ_json,
        GENES = genes_json,
        ROH = ROH_json,
        BED = BED_json)
    
      callModule(downloadObj, id="get_json_filters", output_prefix="Effective_filters.json", output_data=RV$filters_json)
    })
    
    #Load filter setting from json file
    observeEvent(input$load_filters, {
      req(input$load_filters)
      tryCatch({
        filters_json <- read_json(input$load_filters$datapath)
        #load settings from the json file for each filter module
        callModule(loadSettings_filters, "variants_filters", filters_settings$VARIANTS, filters_json$VARIANTS)
        callModule(loadSettings_segregation, "segregation", filters_json$SEGREGATION)
        callModule(loadSettings_GQ, "GQ_filter", filters_json$GQ)
        callModule(loadSettings_filters, "genes_filters", filters_settings$GENES, filters_json$GENES)
        callModule(loadSettings_regions, "ROH_filter", filters_json$ROH)
        callModule(loadSettings_regions, "bed_filter", filters_json$BED)
        
        #force module controls update
        #callModule(observeFilters, "variants_filters", 
        #           filters_settings = filters_settings$VARIANTS, 
        #           variants_df = RV$data$variants_df, 
        #           na_values = app_settings$fill_na$fill_na_vars)
        
        #update filter settings summary data frames
        RV$vars_filters_df <- callModule(getDF_filters, "variants_filters", filters_settings$VARIANTS)
        df1 <- callModule(getDF_segregation, "segregation")
        df2 <- callModule(getDF_GQ, "GQ_filter")
        RV$segregation_filters_df <- rbind(df1, df2)
        RV$genes_filters_df <- callModule(getDF_filters, "genes_filters", filters_settings$GENES)
        df1 <- callModule(getDF_regions, "bed_filter", "custom BED")
        df2 <- callModule(getDF_regions, "ROH_filter", "ROH")
        RV$regions_filters_df <- rbind(df1, df2) 
        
        #Create a success notification
        RV$notifications[["filters_file"]] <- notificationItem(
          text = paste0("Filters configuration loaded correctly"),
          icon = icon("check-circle"),
          status = "success")
        
        #Move to variants tab for filter inspection
        updateTabItems(session, "tabs", "variants_filters_tab")
      }, error=function(cond) {
        RV$notifications[["filters_file"]] <- notificationItem(
         text = paste0("Failed loading filters configuration"),
          icon = icon("exclamation-circle"),
          status = "danger")
      })
      
    })
    
    #Data tables summarizing filters settings
    output$vars_filters_table <- DT::renderDataTable(selection="none", {
      shiny::validate(need(RV$vars_filters_df,"Variants filters not configured"))
      RV$vars_filters_df
    })
    
    output$segregation_filters_table <- DT::renderDataTable(selection="none", {
      shiny::validate(need(RV$segregation_filters_df,"Segregation filters not configured"))
      RV$segregation_filters_df
    })
    
    output$genes_filters_table <- DT::renderDataTable(selection="none", {
      shiny::validate(need(RV$genes_filters_df,"Genes filters not configured"))
      RV$genes_filters_df  
    })
    
    output$regions_filters_table <- DT::renderDataTable(selection="none", {
      shiny::validate(need(RV$regions_filters_df,"Regions filters not configured"))
      RV$regions_filters_df
    })
    
    ## Variants filters tab --------------------------------
    
    observeEvent(input$next_variants_filters, {
      RV$vars_filters_df <- callModule(getDF_filters, "variants_filters", filters_settings$VARIANTS)
      updateTabItems(session, "tabs", "segregation_filters_tab")
    })
    
    output$vars_filters_UI <- renderUI({ 
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded" ))
      filtersVariantsUI("variants_filters", filters_settings$VARIANTS, RV$data$variants_df, app_settings$fill_na$fill_na_vars, filters_settings$TOOLTIPS)
    })
    
      #observeEvent(input$decrypt_button, {
      #  callModule(observeFilters, "variants_filters", 
      #             filters_settings = filters_settings$VARIANTS, 
      #             variants_df = RV$data$variants_df, 
      #             na_values = app_settings$fill_na$fill_na_vars)
      #})

    ## Segregation filters tab -----------------------------------
    observeEvent(input$next_segregation_filters, {
      df1 <- callModule(getDF_segregation, "segregation")
      df2 <- callModule(getDF_GQ, "GQ_filter")
      RV$segregation_filters_df <- rbind(df1, df2)
      updateTabItems(session, "tabs", "genes_filters_tab")
    })
    
    output$segregation_controls <- renderUI({
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded" ))
      segregationUI("segregation", choices_affected = RV$data$values_affected, choices_unaffected=RV$data$values_unaffected)
    })
    
    output$GQfilter_controls <- renderUI({
      GQfilterUI("GQ_filter",maxGQ = RV$maxGQ, defaultGQ=10)    
    })
    
    ## Genes filters tab -----------------------
    
    observeEvent(input$next_genes_filters, {
      RV$genes_filters_df <- callModule(getDF_filters, "genes_filters", filters_settings$GENES)
      updateTabItems(session, "tabs", "regions_filters_tab")
    })
    
    output$genes_filters_UI <- renderUI({ 
      shiny::validate(need(inherits(RV$data, "list"), "No data loaded" ))
      filtersVariantsUI("genes_filters", filters_settings$GENES, RV$data$genes_scores, app_settings$fill_na$fill_na_genes, filters_settings$TOOLTIPS)
    })
    
    #callModule(observeFilters, "genes_filters", 
    #           filters_settings = filters_settings$GENES,
    #           variants_df = RV$data$genes_scores, 
    #           na_values = app_settings$fill_na$fill_na_genes)
    
    ## Regions Filters tab ------------------------
    
    observeEvent(input$next_regions_filters, {
      df1 <- callModule(getDF_regions, "bed_filter", "custom BED")
      df2 <- callModule(getDF_regions, "ROH_filter", "ROH")
      RV$regions_filters_df <- rbind(df1, df2)  
      updateTabItems(session, "tabs", "filters_overview_tab")
    })
    
    output$custom_bed_check <- renderUI({
      shiny::validate(need(inherits(RV$customBed_ranges, "GRanges"), "No custom BED"))
      bedcontrolUI("bed_filter", label = "custom regions", )
    })
    
    #output$ROH_sample_select <- renderUI ({
    #    selectInput("ROH_sample", label = "Select sample:", choices = names(RV$data$ROH_ranges), multiple = FALSE, selected = "AFFECTED_SHARED")
    #})
    
    output$ROH_filters_UI <- renderUI({ 
      shiny::validate(need(inherits(RV$data, "list"), "No ROH data loaded" ))  
      shiny::validate(need(inherits(RV$data$ROH_ranges, "list"), "No ROH regions loaded"))
      bedcontrolUI("ROH_filter", label = "ROH regions", ranges=RV$data$ROH_ranges, slider=TRUE, samples=TRUE, selected_sample = RV$data$all_samples[1])
    })
    observeEvent(RV$data, {
      callModule(configureSlider, "ROH_filter", ranges=RV$data$ROH_ranges, scale_value=1000)  
    })
    
  ## Filter Explorer tab ----------------------------------
  output$filters_explorer_tab <- renderUI({
    shiny::validate(need(RV$filtered_applied, "No filters applied, please set filters first"))
    tagList(
    h3("Summary of filters effect"),
    plotOutput("filters_funnel", width="100%"),
    
    box(title = "Filters impact detail", id = "filters_impact_detail", status = "primary", solidHeader = TRUE,
        collapsible = TRUE, collapsed = TRUE, width = 12,
      fluidRow(
        column(6,withSpinner(plotOutput("summary_variants_filters", width = "100%"))), 
        column(4),
        column(4,withSpinner(plotOutput("summary_genes_filters", width = "100%")))
      ),
      withSpinner(plotSelectedUI("variants_barplot", variables=list("x"=plot_axes[["variants_bar_options"]]), plotly=FALSE)),    
    ),
    
    box(title = "Genes filters scatter", id = "gene_filters_scatter_box", status = "primary", solidHeader = TRUE,
        collapsible = TRUE, collapsed = TRUE, width = 12,
        withSpinner(plotSelectedUI("genes_scatter", variables=list("x"=plot_axes[["genes_axes_options"]],"y"=plot_axes[["genes_axes_options"]]), plotly=TRUE)),
    ),
    
    box(title = "Variants filters scatter", id = "gene_filters_scatter_box", status = "primary", solidHeader = TRUE,
        collapsible = TRUE, collapsed = TRUE, width = 12,
        withSpinner(plotSelectedUI("variants_scatter", variables=list("x"=plot_axes[["variants_axes_options"]],"y"=plot_axes[["variants_axes_options"]]), set_limits=c("x","y"), plotly=FALSE)),
    )
    )
  })
  
  output$filters_funnel <- renderPlot({
    comphet_segregation <- RV$data$comphet_df %>% filter(
      rec_id %in% RV$vars_pass_segregation) %>% 
      gather(key="Variant",value="VarID",v1:v2) %>% select(VarID)
    segregation_vars <- unique(c(RV$vars_pass_segregation, comphet_segregation$varID))
    
    var_groups <- c("exonic","regulatory","structural")
    var_steps <- list(
      "BED" = RV$vars_pass_BED,
      "ROH" = RV$vars_pass_ROH,
      "GQ" = RV$vars_pass_GQ,
      "Variants" = RV$vars_pass_filters$vars,
      "Segregation" = segregation_vars
    )
    
    tot_vars <- list(
      step = NULL,
      group = NULL,
      value = NULL
    )
    steps_vars <- list(
      step = NULL,
      group = NULL,
      value = NULL
    )
    
    diagrams <- list()
    for (group in names(var_groups)) {
      switch(group,
             "exonic" = {df <- RV$data$variants_df %>% filter(consequence %in% exonic_vars)},
             "regulatory" = {df <- RV$data$variants_df %>% filter(consequence %in% reg_vars)},
             "structural" = {df <- RV$data$variants_df %>% filter(var_type %in% sv_vars)},
      )
      tot_vars$step = c(tot_vars$step, "Unfiltered")
      tot_vars$group = c(tot_vars$group, group)
      tot_vars$value = c(tot_vars$value, nrow(df))
      
      for (step in var_steps) {
        df <- df %>% filter(rec_id %in% var_steps[[step]])
        steps_vars$step = c(steps_vars$step, step)
        steps_vars$group = c(steps_vars$group, group)
        steps_vars$value = c(steps_vars$value, nrow(df))
      }
      counts_df <- bind_rows(tot_vars,steps_vars)
    }
    
    ggplot(counts_df, aes(x=step,y=value,label=value)) + 
      geom_bar(stat="identity") + geom_label() + 
      labs(x="Filter steps", y="N records") +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      facet_wrap(~group, ncol = 2)
    
  })  
  
  output$summary_variants_filters <- renderPlot({
      req(nrow(filters_summ_vars())>0)
      #ggplotly(dynamicTicks = T,
          ggplot(filters_summ_vars(), aes(x=Filter,y=Count,fill=Class,label=round(Count,3))) + geom_bar(stat="identity") + 
            geom_label(data=filters_summ_vars() %>% filter(Class=="PASS")) + 
            labs(y="% variants") + theme(axis.text.x = element_text(angle=45, hjust=1))
      #)
  })
  
  output$summary_genes_filters <- renderPlot({
      req(nrow(filters_summ_genes())>0)
      #ggplotly(dynamicTicks = T,
          ggplot(filters_summ_genes(), aes(x=Filter,y=Count,fill=Class, label=round(Count,3))) + geom_bar(stat="identity") + 
            geom_label(data=filters_summ_genes() %>% filter(Class=="PASS")) +
            labs(y="% genes") + theme(axis.text.x = element_text(angle=45, hjust=1))
      #)
  })
  
  callModule(plotModule, "genes_scatter", plot_data = RV$data$genes_scores, missingValues = c(99,-99), plotType = "scatter", plotOptions = list("size" = 1), variables = list("color"="Class", "label"="gene"))
  callModule(plotModule, "variants_scatter", plot_data = RV$data$variants_df, missingValues = c(99,-99), plotType = "bigdata", variables = list("color"="Class", "size" = 1))
  callModule(plotModule, "variants_barplot", plot_data = RV$data$variants_df, plotType = "barplot", variables = list("fill"="Class"), additionalOptions = list(format1, scale_y_sqrt()))
  
  ## Results ----------------------------
    ## Results GENES tab --------------------------
    output$genes_results <- renderUI({ 
      shiny::validate(need(RV$filtered_applied, "No filters applied, please set filters first"))
      tagList(
        fluidRow(
          column(8, withSpinner(plotlyOutput("GADO_rank"))),
          column(4, withSpinner(plotlyOutput("GADO_distribution")))
        ),
        DT::dataTableOutput("genesTable"),
        h3("Custom list genes"),
        DT::dataTableOutput("customGenesTable"),
        hr(),
        fluidRow(column(8), 
                 column(3, offset=1,downloadObjUI(id = "save_results", label = "Download results"))),
        br(),
        box(title = "Applied filters", id = "effective_filters_box", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
            DT::dataTableOutput("effective_filters_table")
        )
      )
    }) 
  
      output$GADO_rank <- renderPlotly({
          genedetail_tab <- as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
          if (length(input$genesTable_rows_selected) > 0) {
              Zscore <- genedetail_tab[input$genesTable_rows_selected, "gado_zscore"]
              GADO_plot <- ggplot(RV$data$genes_scores, aes(x=gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = Zscore, color="red") + geom_vline(xintercept = RV$data$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
          } else {
              GADO_plot <- ggplot(RV$data$genes_scores, aes(x=gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = RV$data$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
          }
          ggplotly( GADO_plot )
      })
      
      output$GADO_distribution <- renderPlotly({
          shiny::validate(need(input$genesTable_rows_selected, "Select a gene"))
          shiny::validate(need(gado_distribution, "Cohort GADO distribution not available"))
          genedetail_tab <- as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
          gene_name <- genedetail_tab[input$genesTable_rows_selected, "gene"]
          Zscore <- genedetail_tab[input$genesTable_rows_selected, "gado_zscore"]
          
          min_value <- min(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name])
          max_value <- max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name])
          req(min_value, max_value)
          GADO_dist <- ggplot(gado_distribution[gado_distribution$Hgnc == gene_name,], aes(x=Zscore)) + 
              geom_histogram(bins = max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name])*10) + 
              geom_vline(xintercept = Zscore, color = "red") +
              scale_x_continuous(breaks=round(seq(min_value, max_value,1)))
          ggplotly( GADO_dist )
       })
      
      candidate_genes_df <- reactive({
          as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% mutate(row_idx = row_number()) %>% select(-Class))
      })
      
      output$genesTable <- DT::renderDataTable({
        na_values <- base::unique(unlist(app_settings$fill_na$fill_na_genes))
        #Check if there are saved vars, if yes change row backgrounds accordingly
        if (is.null(RV$saved_vars)) {
          datatable(candidate_genes_df(), selection="single", options = list(scrollX = TRUE)) %>% formatStyle(names(candidate_genes_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        } else {
          datatable(candidate_genes_df(), selection="single", options = list(scrollX = TRUE)) %>% formatStyle(names(candidate_genes_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values)))) %>%
            formatStyle(  'gene',
                          target = 'row',
                          backgroundColor = styleEqual(saved_genes()$all, rep('yellow', length(saved_genes()$all)))) %>%
            formatStyle(  'gene',
                        target = 'row',
                        backgroundColor = styleEqual(saved_genes()$thisCase, rep('red', length(saved_genes()$thisCase))))
        }
          
      })
      
      genesTable_proxy <- DT::dataTableProxy("genesTable")
      
      customGenesTable_df <- reactive({
          candidate_genes_df() %>% filter(gene %in% RV$custom_genes$gene)
      })
      
      output$customGenesTable <- DT::renderDataTable({
        na_values <- base::unique(unlist(app_settings$fill_na$fill_na_genes))
        if (is.null(RV$saved_vars)) {
          datatable(customGenesTable_df(), selection="single") %>% 
            formatStyle(names(customGenesTable_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        } else {
          datatable(customGenesTable_df(), selection="single") %>% 
            formatStyle(names(customGenesTable_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values)))) %>%
            formatStyle(  'gene',
                          target = 'row',
                          backgroundColor = styleEqual(saved_genes()$all, rep('yellow', length(saved_genes()$all)))) %>%
            formatStyle(  'gene',
                          target = 'row',
                          backgroundColor = styleEqual(saved_genes()$thisCase, rep('red', length(saved_genes()$thisCase))))
        }
      })
      
      observeEvent(input$customGenesTable_rows_selected, {
          row_idx <- customGenesTable_df()[input$customGenesTable_rows_selected,"row_idx"]
          genesTable_proxy %>% selectRows(row_idx) 
      })
      
      output$effective_filters_table <- DT::renderDataTable(selection="none", {
        RV$effective_filters
      })
      
      #Constantly update results object ready for download
      toListenResults <- reactive({
        req(inherits(RV$data, "list"))
        list(RV$data$genes_scores, RV$data$variants_df, RV$data$comphet_df, RV$custom_genes)
      })
      
      observeEvent(toListenResults(), {
        if (nrow(RV$data$comphet_df) > 0) {
          files_to_save = list(
            "genes.tsv" = as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% select(-Class)),
            "customGenes.tsv" = as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% filter(gene %in% RV$custom_genes) %>% select(-Class)),
            "variants.tsv" = as.data.frame(RV$data$variants_df %>% filter(Class == "PASS", gene %in% RV$filtered_genes_list)),
            "comphet.tsv" = as.data.frame(RV$data$comphet_df %>% filter(Class == "PASS", gene %in% RV$filtered_genes_list) %>% 
                                            gather(key="Variant",value = "varID", v1:v2) %>% 
                                            inner_join(., RV$data$variants_df[RV$data$variants_df$Class == "PASS",], by=c("varID"="rec_id")) %>%
                                            select(-Class.x,-Class.y)), 
            "applied_genelist.tsv" = as.data.frame(RV$custom_genes),
            "applied_filters.json" = RV$filters_json)
        } else if (nrow(RV$data$comphet_df) == 0) {
          files_to_save = list(
            "genes.tsv" = as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% select(-Class)),
            "customGenes.tsv" = as.data.frame(RV$data$genes_scores %>% filter(Class == "PASS") %>% filter(gene %in% RV$custom_genes) %>% select(-Class)),
            "variants.tsv" = as.data.frame(RV$data$variants_df %>% filter(Class == "PASS", gene %in% RV$filtered_genes_list)),
            "applied_genelist.tsv" = as.data.frame(RV$custom_genes),
            "applied_filters.json" = RV$filters_json) 
        }
        callModule(downloadObj, id = "save_results",
        output_prefix=input$CaseCode,
        output_data = files_to_save,
        zip_archive = paste0(input$CaseCode, ".results.zip") )
      })
      
    ## Results VARIANTS tab --------------------------
      output$variants_results <- renderUI({ 
        shiny::validate(need(RV$filtered_applied, "No filters applied, please set filters first"))
        tagList(
          fluidRow(
            column(6, selectInput("vars_results_genes", "Show variants for:", choices = c("ALL GENES", "CUSTOM GENE LIST"), selected = "ALL GENES", multiple = FALSE)),
            column(6, actionButton("save_vars_varresult", label = "Save to preferred vars"))
          ),
          br(),
          box(title = "Single variants", id = "vars_results_box", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              DT::dataTableOutput("vars_results_table") ),
          box(title = "Compound het variants", id = "comphet_results_box", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              DT::dataTableOutput("comphet_results_table") )
        )
      })
      
      variants_pass_df <- reactive ({
        if (input$vars_results_genes == "ALL GENES") {
          RV$data$variants_df %>% filter(Class == "PASS")
        } else {
          RV$data$variants_df %>% filter(Class == "PASS", gene %in% RV$custom_genes$gene)
        }
      })
      
      comphet_pass_df <- reactive ({
          comphet_details <- RV$data$comphet_df %>% filter(Class == "PASS") %>% gather(key="variant",value = "varID", v1:v2) %>% select(-variant, -gene) %>% distinct()
          comphet_details <- merge(comphet_details,RV$data$variants_df[RV$data$variants_df$Class == "PASS",], by.x="varID",by.y="rec_id")
          as.data.frame(comphet_details %>% select(-Class.x,-Class.y,) %>% arrange(rec_id))   
      })
      
      output$vars_results_table <- DT::renderDataTable({
         na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
         na_values <- na_values[na_values != 0]
         if (is.null(RV$saved_vars)) {
           datatable(variants_pass_df(), 
                     selection="multiple",
                     options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
             formatStyle(names(variants_pass_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
         } else {
          datatable(variants_pass_df(), 
                   selection="multiple",
                   options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
             formatStyle(  'internal_id',
                         target = 'row',
                         backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(base::unique(saved_vars()$all))))) %>%
             formatStyle(  'rec_id',
                           target = 'row',
                           backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
             formatStyle(names(variants_pass_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
             
          }
      })
      
      output$comphet_results_table <- DT::renderDataTable({
        na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
        if (is.null(RV$saved_vars)) {
          datatable(comphet_pass_df(), 
                    selection="multiple",
                    options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>% 
            formatStyle(names(comphet_pass_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))  
        } else {
          datatable(comphet_pass_df(), 
                  selection="multiple",
                  options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>% 
            formatStyle(  'internal_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(saved_vars()$all)))) %>%
            formatStyle(  'rec_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
            formatStyle(names(comphet_pass_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        }
      })
      
      observeEvent(input$comphet_results_table_rows_selected, {
          gene_name <- base::unique(comphet_pass_df()[input$comphet_results_table_rows_selected, "gene.x"])[1]
          row_idx <- candidate_genes_df()[candidate_genes_df()$gene == gene_name,"row_idx"]
          #message(gene_name, "\t", row_idx)
          genesTable_proxy %>% selectRows(row_idx) 
      })
      
      observeEvent(input$vars_results_table_rows_selected, {
          gene_name <- base::unique(variants_pass_df()[input$vars_results_table_rows_selected, "gene"])[1]
          row_idx <- candidate_genes_df()[candidate_genes_df()$gene == gene_name,"row_idx"]
          #message(gene_name, "\t", row_idx)
          genesTable_proxy %>% selectRows(row_idx)
      })
      
      observeEvent(input$save_vars_varresult, {
        col_names <- colnames(variants_pass_df() %>% select(-Class))
        vars_to_save_1 <- variants_pass_df()[input$vars_results_table_rows_selected,] %>% select(-Class)
        vars_to_save_2 <- comphet_pass_df()[input$comphet_results_table_rows_selected,col_names]
        vars_to_save <- rbind(vars_to_save_1, vars_to_save_2)
        vars_to_save <- cbind(caseID=RV$data$pedigree, vars_to_save)
        
        showModal(
          modalDialog(
            textInput("save_var_note", "Variant note",
                      placeholder = 'Special note for these variants'
            ),
            span('You can input here any short note about the saved variants',
                 'Note that commas are not allowed and will be converted to underscores'),
            footer = tagList(
              actionButton("ok_note", "OK")
            )
          )
        )
      })
      
    ## PanelApp and genes lists tab ------------------------------
      output$gene_lists_results <- renderUI({ 
        shiny::validate(need(RV$filtered_applied, "No filters applied, please set filters first"))
        tagList(
          box(title = "PanelApp panels", id = "PanelApp_panels", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = FALSE,
              DT::dataTableOutput("PanelApp_panels_table"),
              DT::dataTableOutput("PanelApp_genes_table") ),
          box(title = "ClinVar pathogenic/likely pathogenic", id = "ClinVar_path", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              DT::dataTableOutput("ClinVar_table") ),
          box(title = "Other genes lists", id = "Other_genes_lists", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              DT::dataTableOutput("geneLists_table"),
              DT::dataTableOutput("geneLists_genes_table") )
        )
      })
      
      output$PanelApp_panels_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), rownames=F, {
          PanelApp_panels_df()
      })
      
      PanelApp_genes_df <- reactive ({
          shiny::req(input$PanelApp_panels_table_rows_selected)
          panelID <- PanelApp_panels_df()[input$PanelApp_panels_table_rows_selected, "id"]
          genes_df <- as.data.frame(PanelApp_genes %>% filter(panel_idx == panelID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., RV$data$genes_scores, by = c("entity_name" = "gene")))
          genes_df <- genes_df[order(genes_df$gado_zscore, decreasing = TRUE),]
      })
      
      output$PanelApp_genes_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
          shiny::validate(need(input$PanelApp_panels_table_rows_selected, "Select one panel"))
          PanelApp_genes_df()       
      })
      
      output$ClinVar_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
          ClinVar_df()
      })
      
      output$geneLists_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
          geneLists_df()
      })
      
      geneList_genes_df <- reactive ({
          shiny::req(input$geneLists_table_rows_selected)
          listID <- geneLists_df()[input$geneLists_table_rows_selected, "id"]
          genes_df <- as.data.frame(geneLists_genes %>% filter(genelist_idx == listID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., RV$data$genes_scores, by = c("entity_name" = "gene")))
          genes_df <- genes_df[order(genes_df$gado_zscore, decreasing = TRUE),]
      })
      
      output$geneLists_genes_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
          shiny::validate(need(input$geneLists_table_rows_selected, "Select one gene list"))
          geneList_genes_df()
      })
      
      observeEvent(input$ClinVar_table_rows_selected, {
          gene_name <- ClinVar_df()[input$ClinVar_table_rows_selected, "gene"]
          row_idx <- candidate_genes_df()[candidate_genes_df()$gene == gene_name,"row_idx"]
          genesTable_proxy %>% selectRows(row_idx) 
      })
      
      observeEvent(input$PanelApp_genes_table_rows_selected, {
          gene_name <- PanelApp_genes_df()[input$PanelApp_genes_table_rows_selected, "entity_name"]
          row_idx <- candidate_genes_df()[candidate_genes_df()$gene == gene_name,"row_idx"]
          genesTable_proxy %>% selectRows(row_idx) 
      })
      
      observeEvent(input$geneLists_genes_table_rows_selected, {
          gene_name <- geneList_genes_df()[input$geneLists_genes_table_rows_selected, "entity_name"]
          row_idx <- candidate_genes_df()[candidate_genes_df()$gene == gene_name,"row_idx"]
          genesTable_proxy %>% selectRows(row_idx) 
      })
      
  ## Gene detail tab -----------------------------
    # Get gene name of the selected candidate gene
    gene_name <- reactive ({
        symbol = candidate_genes_df()[input$genesTable_rows_selected, "gene"]
        updateSelectInput(session = session, inputId = "chr_coverage",selected = genes_bed$V1[genes_bed$V4 == symbol])
        message("Selected gene: ", symbol)
        return(symbol)
      })
      
    output$gene_detail_ui <- renderUI({
        shiny::validate(need(gene_name(), "No gene selected, please select a gene from gene results first"))
        tagList(
          fluidRow(
            column(2, h4("Gene symbol:")) , column(9, offset = 1, textOutput("Gene_symbol"))
          ),
          fluidRow(
            column(2, h4("Gene name:")) , column(9, offset = 1, textOutput("Gene_name"))
          ),
          fluidRow(
            column(2, h4("Link to GTeX:")) , column(9, offset = 1, uiOutput("GTeX_link"))
          ),
          fluidRow(br()),
          h3("Scores for the selected gene"),
          DT::dataTableOutput("geneDetail"),
          h3("Filtered variants for the selected gene"),
          DT::dataTableOutput("variantsTable"),
          h3("Filtered compound hets for the selected gene"),
          DT::dataTableOutput("comphetTable"),
          hr(),
          fluidRow(
            column(4, actionButton("save_vars_genedetail",label = "Save to preferred vars")),
            column(4, downloadObjUI("get_igv_session", label = "Download IGV session")),
            column(4, textOutput("igv_region"))
          ),
          br(),
          fluidRow(column(4, uiOutput("go_to_venus"))),
          
          hr(),
          box(title = "Regulatory regions details", id = "reg_regions_box", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              uiOutput("reg_regions_details")
          ),
          box(title = "PanelApp and ClinVar", id = "panelapp_clinvar_details", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              h3("PanelApp panels"),
              DT::dataTableOutput("panelapp_detail_tab"),
              h3("ClinVar associated diseases"),
              DT::dataTableOutput("clinvar_detail_tab"),
              h3("Interesting gene lists"),
              DT::dataTableOutput("genelists_detail_tab") ),
          box(title = "Associated HPO terms", id = "HPO_terms_details", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              DT::dataTableOutput("hpo_detail_tab")),
          box(title = "Pathways and ontology", id = "path_go_details", status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              div(DT::dataTableOutput("geneInfo"),style="font-size: 75%") ),
          box(title = "GTeX median expression", id = "gtex_expression", status = "info", solidHeader = TRUE, width=12,
              collapsible = TRUE, collapsed = TRUE,
              plotOutput("gtex_plot") )
        )
      })
    
    venus_link <- reactive({
      selected_vars <- gene_comphet_vars_df()[input$comphetTable_rows_selected, ]
      aa_changes <- selected_vars$aa_change[selected_vars$consequence == "missense_variant"]
      selected_vars <- gene_vars_df()[input$variantsTable_rows_selected, ]
      aa_changes <- base::unique(c(aa_changes, selected_vars$aa_change[selected_vars$consequence == "missense_variant"]))
      links <- list()
      for (aa_change in aa_changes)
        if (!is.na(aa_change)) {
          aa_change <- unlist(strsplit(aa_change,","))[1]
          id <- aa_change
          transcript_id <- unlist(strsplit(aa_change,":"))[1]
          aa_change <- unlist(strsplit(aa_change,":"))[2]
          aa_change <- gsub("p\\.","",aa_change)
          aa <- str_extract_all(aa_change, "[A-Za-z]{3}")[[1]]
          pos <- str_extract(aa_change, "\\d+")
          aa_change <- paste0(aa_codes[[aa[1]]],pos,aa_codes[[aa[2]]])
          url_link <- paste0("https://michelanglo.sgc.ox.ac.uk/venus_transcript?enst=",transcript_id, "&mutation=", aa_change, "&redirect")
          links[[id]] <- url_link
        }
      return(links)
    })
    
    output$Gene_symbol <- renderText({
        shiny::validate(need(gene_name() != "", 'No gene selected'))
        gene_name()
    })
    
    output$Gene_name <- renderText({
        shiny::validate(need(gene_name() != "", 'No gene selected'))
        genes_info[genes_info$symbol == gene_name(), "name"]
    })
    
    output$GTeX_link <- renderUI({
        shiny::validate(need(gene_name() != "", 'No gene selected'))
        ensg_id <- genes_info[genes_info$symbol == gene_name(), "ensembl_gene_id"]
        tags$a(href=paste0("https://gtexportal.org/home/gene/", gene_name()), paste0(gene_name(), "(", ensg_id, ")"), target="_blank")
    })
    
    output$geneDetail <- DT::renderDataTable(selection="none", {
        candidate_genes_df()[input$genesTable_rows_selected,]
    })
    
    output$geneInfo <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        shiny::validate(need(gene_name() != "", 'No gene selected'))
        selected <- list()
        
        for (n in names(gene_anno)) {
            id_list <- grep(gene_name(), gene_anno[[n]])
            selected[[n]] <- names(gene_anno[[n]])[id_list]
        }
        Additional_info <- data.frame(lapply(selected, "length<-", max(lengths(selected))), stringsAsFactors = F)
    })
    
    output$gtex_plot <- renderPlot ({
        gene_exp <- GTeX_data[GTeX_data$Description == gene_name(),]
        shiny::validate(need(nrow(gene_exp)>0, "Gene not found in GTeX"))
        gene_exp <- gather(gene_exp, key="tissue", value="median_TPM", 3:ncol(gene_exp))
        ggplot(gene_exp, aes(x=tissue, y=median_TPM)) + geom_bar(stat="identity") + format1
    })
    
    gene_comphet_vars_df <- reactive({
        comphet_details <- RV$data$comphet_df %>% filter(Class == "PASS", gene == gene_name()) %>% gather(key="variant",value = "varID", v1:v2) %>% select(-gene, -variant) %>% distinct()
        comphet_details <- merge(comphet_details,RV$data$variants_df[RV$data$variants_df$Class == "PASS",], by.x="varID",by.y="rec_id")
        as.data.frame(comphet_details %>% select(-Class.x,-Class.y,) %>% arrange(rec_id))
    })
    
    gene_vars_df <- reactive ({
        as.data.frame(RV$data$variants_df %>% filter(Class == "PASS", gene == gene_name(), rec_id %nin% gene_comphet_vars_df()$varID)) 
    })
    
    output$variantsTable <- DT::renderDataTable({
        #shiny::validate(need(gene_name != "", 'No gene selected'))
        na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
        na_values <- na_values[na_values != 0]
        if (is.null(RV$saved_vars)) {
          datatable(gene_vars_df(), 
                    selection="multiple",
                    options = list(scrollX = TRUE)) %>%
            formatStyle(names(gene_vars_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))   
        } else {
          datatable(gene_vars_df(), 
                selection="multiple",
                options = list(scrollX = TRUE)) %>%
            formatStyle(  'internal_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(saved_vars()$all)))) %>%
            formatStyle(  'rec_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
            formatStyle(names(gene_vars_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        }
    })
    
    output$comphetTable <- DT::renderDataTable({
        #shiny::validate(need(gene_name != "", 'No gene selected'))
        shiny::validate(need(nrow(gene_comphet_vars_df())>0, 'No compound het variants'))
        na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
        na_values <- na_values[na_values != 0]
        
        if (is.null(RV$saved_vars)) {
          datatable(gene_comphet_vars_df(), 
                    selection="multiple",
                    options = list(scrollX = TRUE)) %>%
            formatStyle(names(gene_comphet_vars_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values)))) 
        } else {
          datatable(gene_comphet_vars_df(), 
                selection="multiple",
                options = list(scrollX = TRUE)) %>%
            formatStyle(  'internal_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(saved_vars()$all)))) %>%
            formatStyle(  'rec_id',
                          target = 'row',
                          backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
            formatStyle(names(gene_comphet_vars_df()), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        }
    })
    
    observeEvent(input$save_vars_genedetail, {
      col_names <- colnames(RV$data$variants_df %>% select(-Class))
      vars_to_save_1 <- gene_vars_df()[input$variantsTable_rows_selected,] %>% select(-Class)
      vars_to_save_2 <- gene_comphet_vars_df()[input$comphetTable_rows_selected,col_names]
      vars_to_save <- rbind(vars_to_save_1, vars_to_save_2)
      
      #Remove samples GQ and GT cols since they are different for each case
      #GTs and samples order are collapsed and saved in new cols
      GT_cols <- grep("GT_",colnames(vars_to_save))
      GQ_cols <- grep("GQ_",colnames(vars_to_save))
      vars_to_save$GTs <- apply( as.data.frame(vars_to_save[ , GT_cols ]) , 1 , paste, collapse=",")
      vars_to_save$samples_order <- paste(gsub("GT_","",colnames(vars_to_save)[GT_cols]), collapse=",")
      vars_to_save <- vars_to_save %>% select(-GT_cols, -GQ_cols)
      
      RV$vars_to_save <- cbind(caseID=RV$data$pedigree, vars_to_save)
      
      showModal(
        modalDialog(
        textInput("save_var_note", "Variant note",
                  placeholder = 'Special note for these variants'
        ),
        span('You can input here any short note about the saved variants',
             'Note that commas are not allowed and will be converted to underscores'),
        footer = tagList(
          actionButton("ok_note", "OK")
        )
        )
      )
    })
    
    observeEvent(input$ok_note, {
      note_text <- gsub(",","_",input$save_var_note)
      RV$vars_to_save <- RV$vars_to_save %>% tibble::add_column(note=note_text, .after = "caseID")
      RV$saved_vars <- rbind(RV$saved_vars, RV$vars_to_save) %>% distinct()
      removeModal()
    })

    output$panelapp_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE),  {
      req(gene_name() != "" & !is.null(gene_name()))  
      gene_df <- PanelApp_genes[PanelApp_genes$entity_name == gene_name(), c("entity_name","panel_idx","confidence_level")]
        shiny::validate(need(nrow(gene_df)>0, "No PanelApp panels for this gene"))
        
        df <- merge(gene_df, PanelApp_panels_df(), by.x="panel_idx", by.y="id")
        df[,c("entity_name","confidence_level","panel_idx","name","disease_group","version","relevant_disorders")]
    })
    
    output$clinvar_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        ClinVar_df()[ClinVar_df()$gene == gene_name(),]
    })
    
    output$genelists_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE),  {
        listID <- geneLists_genes[geneLists_genes$entity_name == gene_name(), "genelist_idx"]
        geneLists_df()[geneLists_df()$id %in% listID,]
    })
    
    output$hpo_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
      HPO_ids <- base::unique(HPO_genes$HPO_id[HPO_genes$gene == gene_name()])
      HPO_table <- HPO_obo %>% filter(HPO_id %in% HPO_ids) %>% distinct()
    })
    
    IGV_session <- reactive ({
        shiny::validate(need(!is.null(input$variantsTable_rows_selected) | !is.null(input$comphetTable_rows_selected), "No rows selected"))
        chromosome <- base::unique(c(gene_comphet_vars_df()$chr[input$comphetTable_rows_selected], gene_vars_df()$chr[input$variantsTable_rows_selected]))
        start_pos <- base::unique(c(gene_comphet_vars_df()$start[input$comphetTable_rows_selected], gene_vars_df()$start[input$variantsTable_rows_selected]))
        end_pos <- base::unique(c(gene_comphet_vars_df()$end[input$comphetTable_rows_selected], gene_vars_df()$end[input$variantsTable_rows_selected]))
        
        shiny::validate(need(length(chromosome) == 1, "Please select vars on the same chromosome"))
        shiny::validate(need(!is.null(RV$data$bam_files), "No information on BAM files loaded"))
        
        start_pos <- min(start_pos) - 100
        end_pos <- max(end_pos) + 100
        region <- paste0(chromosome,":",start_pos,"-",end_pos)
        out_file <- paste0(chromosome,"_",start_pos,"-",end_pos,".xml")
        RV$selected_vars_region <- region
        
        igv_xml <- makeIGVxml(region, 
                              RV$data$affected_samples, 
                              RV$data$unaffected_samples,
                              VCF_files = RV$data$vcf_files,
                              SV_files = RV$data$sv_files,
                              BAM_files = RV$data$bam_files )
        
        #jigv_command <- paste("jigv --region", region,
        #                      paste('"',BAM_dir, RV$data$affected_samples, '.bam#', RV$data$affected_samples, '_affected"', sep="", collapse=" " ),
        #                      paste('"',BAM_dir, RV$data$unaffected_samples, '.bam#', RV$data$affected_samples, '_unaffected"', sep="", collapse=" " ),
        #                      paste0(VCF_dir, RV$data$pedigree, ".PASS.NORM.vcf.gz"), 
        #                      SV_file, collapse=" ")
        callModule(downloadObj, id="get_igv_session", output_prefix= out_file, output_data=igv_xml, col_names=FALSE)
        
        list(outfile=out_file, session_xml=igv_xml)
    })
    
    output$go_to_venus <- renderUI ({
      req(length(venus_link())>0)
      #message(str(venus_link()))
      mybutton <- list()
      for (i in 1:length(venus_link())) {
        mybutton[[paste0("venus_",i)]] <- actionButton(paste0("venus_",i), 
                                              paste0("VENUS prediction for ", names(venus_link())[i]), 
                                              onclick =paste0("window.open('", venus_link()[[i]], "', '_blank')"))
      }
      col_dimension <- ceiling(12 / length(venus_link()))
      buttons_row <- NULL
      for (n in names(mybutton)) {
        buttons_row <- tagList(buttons_row,tagList(column(col_dimension, mybutton[[n]], align="center")))
      }
      fluidRow(buttons_row)
    })
    
    selected_var <- reactive({
      var1 <- NULL
      var2 <- NULL
      if(nrow(gene_comphet_vars_df())>0 & length(input$comphetTable_rows_selected) > 0) {
        var1 <- gene_comphet_vars_df()[input$comphetTable_rows_selected, c("start","reg_id")]
      }
      if(nrow(gene_vars_df()) & length(input$variantsTable_rows_selected) > 0) {
        var2 <- gene_vars_df()[input$variantsTable_rows_selected, c("start","reg_id")]
      }
      vars <- rbind(var1, var2)
      
      if (is.null(vars)) {
        selected_var <- list(status=0, var_pos=0, reg_ids=0)
      } else if (nrow(vars)==0) {
        selected_var <- list(status=0, var_pos=0, reg_ids=0)
      } else if (nrow(vars)==1) {
        if (is.na(vars$reg_id)) {
          selected_var <- list(status=0, var_pos=0, reg_ids=0)
        } else {
        reg_ids <- unlist(strsplit(vars$reg_id, ","))
        selected_var <- list(status=1, var_pos=vars$start, reg_ids=reg_ids)
        callModule(SQlite_query,"GREENDB_query",db=GREENDB_file, var_position=selected_var$var_pos, regions=selected_var$reg_ids)
        }
      } else if (nrow(vars) > 1) {
        selected_var <- list(status=2, var_pos=0, reg_ids=0)
      }
  
      message("status:", selected_var$status, "; var_pos:", selected_var$var_pos, "; reg_ids:", selected_var$reg_ids)
      return(selected_var)
    })
    
    output$reg_regions_details <- renderUI({
      if(selected_var()$status == 1) {
        query_results_UI("GREENDB_query", boxes=TRUE)
      } else {
        verbatimTextOutput("reg_region_detail_error")
      }
    })
    
    output$reg_region_detail_error <- renderText({
      text <- NULL
      if (selected_var()$status == 2) {
        text <- "You can only select a single variant and it must have regulatory region annotations"
      } else if(selected_var()$status == 0) {
        text <- "No variant selected or the variant has no regulatory region annotations"
      }
      return(text)
    })
    
    output$igv_region <- renderText({
        paste0("Generated session: ", IGV_session()$outfile)
    })
    
  ## Expansion Hunter tab -----------------------
    output$exphunter_selection <- renderUI({
        selectInput("exphunter_loci", h3("Gene:"), choices = sort(base::unique(RV$data$ExpHunter$PEDIGREE$LocusId)), multiple = FALSE)
    })
    
    output$exphunter_variants <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        RV$data$ExpHunter$PEDIGREE %>% filter(LocusId == input$exphunter_loci) %>% select(VariantId,RepeatUnit,ReferenceRegion) %>% distinct()
    })
    
    output$exphunter_plot <- renderPlot({ 
        shiny::validate(need(input$exphunter_loci != "", "No gene selected"))
        ggplot(RV$data$ExpHunter$PEDIGREE[RV$data$ExpHunter$PEDIGREE$LocusId == input$exphunter_loci,], aes(y=sampleID, x=Genotype, fill=Coverage)) + 
            geom_tile() + 
            facet_grid(group~VariantId, scales="free") +
            format1
    })
    
  ## Explore coverage tab ----------------------
    output$gene_cov_select <- renderUI({
        genes_list <- sort(base::unique(genes_bed$V4[genes_bed$V1 == input$chr_coverage]))
        if (!is.null(gene_name())) {
            selected_gene = gene_name()
        } else {
            selected_gene = "NO_GENE"
        }
        selectInput("gene_coverage", "Gene: ", choices = c("NO_GENE",genes_list), selected = selected_gene)
    })
    
    observeEvent(input$plot_coverage, {
        filename <- paste0(Coverage_dir,"/indexcov-", input$chr_coverage, ".bed.gz")
        shiny::validate(need(file_test("-f", filename), "No coverage data"))
        
        RV$cov_file <- filename
    })
    
    output$smooth_dimension <- renderText({
        paste0("Window size: ", as.numeric(input$smooth_factor) * 16.4, " kb")
    })
    
    output$coverage_plot <- renderPlotly({
        req(RV$cov_file)
        cov_file <- gzfile(RV$cov_file)
        cov_data <- read.table(cov_file, sep="\t", header=T, comment.char = "@", stringsAsFactors = F)
        cov_df <- cov_data %>% select(all_of(c("X.chrom","start","end", RV$data$all_samples))) %>% gather(.,key="sample",value="norm_cov", 4:(3+length(RV$data$all_samples)))
        cov_df$middle_point <- round(cov_df$start + ((cov_df$end - cov_df$start) / 2))
        
        #if a gene is selected restrict data to just this gene
        if (input$gene_coverage != "NO_GENE") {
            plot_chr <- genes_bed$V1[genes_bed$V4 == input$gene_coverage]
            plot_start <- genes_bed$V2[genes_bed$V4 == input$gene_coverage] - as.numeric(input$cov_region_pad)*1000
            plot_end <- genes_bed$V3[genes_bed$V4 == input$gene_coverage] + as.numeric(input$cov_region_pad)*1000
            cov_df <- cov_df %>% filter(X.chrom == plot_chr, middle_point >= plot_start & middle_point <= plot_end)   
        }
        
        #apply smooth factor to coverage data
        d <- as.data.frame(cov_df %>% group_by(sample,G=trunc(0:(n()-1)/input$smooth_factor)) %>% summarise(mean=mean(norm_cov), pos=median(middle_point), .groups="keep"))
        
        cov_plot <- ggplot(d, aes(x=pos/100000, y=mean, color=sample)) + 
            geom_line(size=0.5, alpha=0.5) + scale_y_continuous(limits=c(0,3), breaks=0:3) + 
            theme(axis.text.x=element_text(size=8, angle=45, hjust=1)) + 
            labs(x="position (x 100kb)", title=paste0("Chromosome: ", input$chr_coverage), subtitle = paste0("Gene: ", input$gene_coverage))
        
        ggplotly(cov_plot, tooltip = c("color","y"), dynamicTicks = T)
    })
    
    outputOptions(output, "vars_filters_UI", suspendWhenHidden = FALSE)
    outputOptions(output, "segregation_controls", suspendWhenHidden = FALSE)
    outputOptions(output, "GQfilter_controls", suspendWhenHidden = FALSE)
    outputOptions(output, "genes_filters_UI", suspendWhenHidden = FALSE)
    outputOptions(output, "custom_bed_check", suspendWhenHidden = FALSE)
    outputOptions(output, "ROH_filters_UI", suspendWhenHidden = FALSE)    
    
  ## Known variants tab --------------------------
    output$known_clinvar_tab <- DT::renderDataTable({
      ok_cols <- c(vis_cols,
                   colnames(RV$data$known_clinvar)[grep("GT_",colnames(RV$data$known_clinvar))])
      hide_cols <- setdiff(colnames(RV$data$known_clinvar), ok_cols)
      hide_cols_idx <- which(colnames(RV$data$known_clinvar) %in% hide_cols)
      datatable(
        RV$data$known_clinvar, extensions = 'Buttons', selection="single", filter = "top",
        options = list(
          scrollX = TRUE,
          dom = 'lBfrtip',
          buttons = c('copy', 'csv', 'excel','colvis'),
          columnDefs = list(
            list(targets = hide_cols_idx, visible = FALSE)
          )
        )
      )
    })
    
    output$known_cosmic_tab <- DT::renderDataTable({
      ok_cols <- c(vis_cols,
                   colnames(RV$data$known_cosmic)[grep("GT_",colnames(RV$data$known_cosmic))])
      hide_cols <- setdiff(colnames(RV$data$known_cosmic), ok_cols)
      hide_cols_idx <- which(colnames(RV$data$known_clinvar) %in% hide_cols)
      datatable(
        RV$data$known_cosmic, extensions = 'Buttons', selection="none", filter = "top",
        options = list(
          scrollX = TRUE,
          dom = 'lBfrtip',
          buttons = c('copy', 'csv', 'excel','colvis'),
          columnDefs = list(
            list(targets = hide_cols, visible = FALSE)
          )
        )
      )
    })
    
    clinvar_var_info <- reactive({
      req(input$known_clinvar_tab_rows_selected)
      infos <- getClinvarInfo(RV$data$known_clinvar$known_ids[input$known_clinvar_tab_rows_selected])
      return(infos)
    })
    
    output$clinvar_var_sig <- renderText({
      if (is.null(clinvar_var_info()$significance) | is.na(clinvar_var_info()$significance) | clinvar_var_info()$significance == "") {
        return ("NONE")
      } else {
        return(clinvar_var_info()$significance)
      }
    })
    
    output$clinvar_var_pheno <- renderText({
      if (is.null(clinvar_var_info()$phenotype) | is.na(clinvar_var_info()$phenotype) | clinvar_var_info()$phenotype == "") {
        return ("NONE")
      } else {
        return(clinvar_var_info()$phenotype)
      }
    }) 
  ## Preferred vars -------------------------
    observeEvent(input$preferred_vars_remove, {
      shiny::validate(need(input$preferred_vars_tab_rows_selected, "Select a variant"))
      RV$saved_vars <- RV$saved_vars[-input$preferred_vars_tab_rows_selected,]
      if(nrow(RV$saved_vars) == 0) {RV$saved_vars <- NULL}
    })
    
    observeEvent(input$preferred_vars_reset, {
      RV$saved_vars <- NULL
    })
    
    observeEvent(input$preferred_vars_file, {
      req(input$preferred_vars_file)
      tryCatch({
        preferred_vars_df <- read.table(input$preferred_vars_file$datapath,header=T,sep=",",stringsAsFactors = F, as.is=T)
        RV$saved_vars <- preferred_vars_df %>% select(-X)
        RV$notifications[["preferred_vars"]] <- notificationItem(
          text = paste0("Loaded ", nrow(RV$saved_vars), " preferred variants"),
          icon = icon("check-circle"),
          status = "success")
      }, error=function(cond) {
        RV$notifications[["preferred_vars"]] <- notificationItem(
          text = paste0("Failed loading preferred variants"),
          icon = icon("exclamation-circle"),
          status = "danger")
      })
      
    })
    
    output$preferred_vars_tab <- DT::renderDataTable({
      req(RV$saved_vars)
      ok_cols <- c("caseID", "note", vis_cols,
                   colnames(RV$saved_vars)[grep("GT_",colnames(RV$saved_vars))])
      hide_cols <- setdiff(colnames(RV$saved_vars), ok_cols)
      hide_cols_idx <- which(colnames(RV$saved_vars) %in% hide_cols)
      datatable(
        RV$saved_vars, extensions = 'Buttons', selection="multiple",
        editable = list(target = 'column', disable = list(columns = c(1,3:ncol(RV$saved_vars)))),
        options = list(
          pageLength = 20, lengthMenu = c(10, 20, 50, 100, 200),
          scrollX = TRUE,
          dom = 'lBfrtip',
          buttons = c('copy', 'csv', 'excel','colvis'),
          autoWidth = TRUE,
          columnDefs = list(
            list(width = '400px', targets = c(2)),
            list(targets = hide_cols_idx, visible = FALSE)
          )
        )
      )
    })
    
    observeEvent(input$preferred_vars_tab_cell_edit,{
      RV$saved_vars[input$preferred_vars_tab_cell_edit$row,input$preferred_vars_tab_cell_edit$col] <<- input$preferred_vars_tab_cell_edit$value
    })
    
    saved_genes <- reactive ({
      thisCase <- "HACKTHIS"
      all <- "NOTHING_SAVED"
      if(!is.null(RV$saved_vars)) {
        thisCase <- c(thisCase, base::unique(RV$saved_vars$gene[RV$saved_vars$caseID == RV$data$pedigree]))
        all <- base::unique(RV$saved_vars$gene)
      }
      return(list(thisCase=thisCase,all=all))
    })
    
    saved_vars <- reactive ({
      thisCase <- "HACKTHIS"
      all <- "NOTHING_SAVED"
      if(!is.null(RV$saved_vars)) {
        thisCase <- c(thisCase, base::unique(RV$saved_vars$rec_id[RV$saved_vars$caseID == RV$data$pedigree]))
        all <- base::unique(RV$saved_vars$internal_id)
      }
      return(list(thisCase=thisCase,all=all))
    })
  ## Cohort mode ------------------
    cohort_results <- observeEvent(input$apply_cohort, {
      req(RV$vars_filters_df, RV$segregation_filters_df, RV$genes_filters_df, RV$regions_filters_df)
      n <- length(input$cohort_samples)
      message("Started cohort analysis for ", n, " samples")
      #updateTabItems(session, "tabs", "cohort_analysis")
      withProgress(message = 'Cohort analysis', value = 0, {
                   
        cohort <- list()
        RV$cohort_files <- list(applied_genelist.txt = as.data.frame(RV$custom_genes),
                             applied_filters.json = RV$filters_json)
        registerDoParallel(cores=input$cohort_threads)
        sample_results <- foreach (sampleID=input$cohort_samples, .packages = c("dplyr","tidyr"), .verbose = T, .combine = combine_df, .multicombine = TRUE) %dopar% {
          #sampleID <- gsub("\\.RData[.enc]*","",file, perl = T)
          incProgress(1/n, detail = paste("Doing", sampleID))
          
          cohort$vars_pass_GQ <- NULL
          cohort$vars_pass_BED <- NULL
          cohort$vars_pass_ROH <- NULL 
          cohort$vars_pass_segregation <- NULL
          cohort$vars_pass_filters <- NULL
          cohort$vars_filters_df <- NULL
          cohort$segregation_filters_df <- NULL
          cohort$genes_filters_df <- NULL
          cohort$regions_filters_df <- NULL
          cohort$filters_json <- NULL
          
          #Load data from plain or encrypted object
          if (file_test("-f", paste0(data_dir,"/",sampleID,".RData.enc"))) {
            cohort$data <- decrypt_datafile(paste0(data_dir,"/",sampleID,".RData.enc"), pwd = input$pwd)
          } else if (file_test("-f", paste0(data_dir,"/",sampleID,".RData"))) {
            cohort$data <- readRDS(paste0(data_dir,"/",sampleID,".RData"))
          }
          
          shiny::validate(need(inherits(cohort$data, "list"), FALSE))
          cohort$data$variants_df <- cohort$data$variants_df %>% replace_na(app_settings$fill_na$fill_na_vars)
          cohort$data$segregation_df$sup_dnm[cohort$data$segregation_df$sup_dnm < 0] <- 0
          cohort$data$genes_scores <- cohort$data$genes_scores %>% replace_na(app_settings$fill_na$fill_na_genes)
          #cohort$data$genes_scores$cohort_norm_Z <- apply(cohort$data$genes_scores,1,function(x) computeNormZ(x["gene"],x["gado_zscore"]))
          
          cohort$GQ_cols_all <- which(colnames(cohort$data$variants_df) %in% paste("GQ", cohort$data$all_samples, sep="_"))
          cohort$GQ_cols_affected <- which(colnames(cohort$data$variants_df) %in% paste("GQ", cohort$data$affected_samples, sep="_"))
          
          #Apply same filters
          cohort$vars_pass_filters <- callModule(getPASSVars_filters, "variants_filters", filters_settings$VARIANTS, cohort$data$variants_df, cohort$data$comphet_df)
          #message("PASS vars ", length(cohort$vars_pass_filters$vars))
          cohort$genes_pass_filters <- callModule(getPASSGenes_filters, "genes_filters", filters_settings$GENES, cohort$data$genes_scores)
          #message("PASS genes ", length(cohort$genes_pass_filters))
          cohort$vars_pass_segregation <- callModule(segregationModule, "segregation", segregation_df = cohort$data$segregation_df, cols_names = segregation_cols)
          #message("PASS segregation ", length(cohort$vars_pass_segregation))
          cohort$vars_pass_GQ <- callModule(GQfilterModule, "GQ_filter", 
                                        variants_df = cohort$data$variants_df, 
                                        GQ_cols = cohort$GQ_cols_all, 
                                        affected_cols = cohort$GQ_cols_affected,
                                        exclude_var_type = sv_vars)
          #message("PASS GQ ", length(cohort$vars_pass_GQ))
          
          if (inherits(RV$customBed_ranges, "GRanges")) {
            cohort$vars_pass_BED <- callModule(bedfilterModule,"bed_filter",variants_ranges=cohort$data$variants_ranges, bed_ranges=RV$customBed_ranges)
          } else {
            cohort$vars_pass_BED <- cohort$data$variants_df$rec_id
          }
          #message("PASS BED ", length(cohort$vars_pass_BED))
          
          #Variants results
          ## FILTER VARS
          #Get the final list of vars in accepted comphet intersecting filters and segregation
          cohort$comphet_final_vars <- cohort$data$comphet_df %>% filter(
            rec_id %in% intersect(cohort$vars_pass_segregation,cohort$vars_pass_filters$comphet)) %>% 
            gather(key="Variant",value="VarID",v1:v2) %>% select(VarID)
          
          #GET THE FINAL LIST OF ACCEPTED VARS
          #the final list of segregating accepted vars is now equal to
          #single vars passing filters and segregation
          #vars part of a comphet passing all filters + segregation
          cohort$accepted_vars_list <- base::unique(intersect(intersect(intersect(
            cohort$vars_pass_GQ, 
            cohort$vars_pass_BED),
            cohort$vars_pass_filters$vars),
            cohort$vars_pass_segregation))
          
          cohort$accepted_vars_list <- base::unique(c(cohort$accepted_vars_list, cohort$comphet_final_vars$VarID))
          #message("Final PASS count ", length(cohort$accepted_vars_list))
          
          #return variants dataframe with updated Class column (PASS/FILTER)
          cohort$data$variants_df <- as.data.frame(cohort$data$variants_df %>% 
                                                    filter(rec_id %in% cohort$accepted_vars_list) )
          
          #comphet results
          cohort$data$comphet_df <- as.data.frame(cohort$data$comphet_df %>% 
                                  filter(v1 %in% cohort$data$variants_df$rec_id,
                                          v2 %in% cohort$data$variants_df$rec_id,
                                          rec_id %in% intersect(cohort$vars_pass_segregation,cohort$vars_pass_filters$comphet)))
    
          #genes results
          cohort$filtered_vars_list <- base::unique(c(
            cohort$data$variants_df$rec_id, 
            cohort$data$comphet_df$rec_id)) 
          
          #Filters on gene scores are applied
          cohort$data$genes_df <- as.data.frame(cohort$data$genes_df %>% 
                                                  filter(variants %in% cohort$filtered_vars_list,
                                                         gene %in% cohort$genes_pass_filters))    
          
          #genes scores
          cohort$data$genes_scores <- as.data.frame(cohort$data$genes_scores %>% 
                                                      filter(gene %in% cohort$data$genes_df$gene))
          
          #Create results table
          output <- list()
          cohort$data$genes_scores$caseID <- sampleID
          cohort$data$variants_df$caseID <- sampleID
          cols_to_remove <- c(paste("GQ", cohort$data$all_samples, sep="_"),
                              paste("GT", cohort$data$all_samples, sep="_"))
          cohort$data$variants_df <- cohort$data$variants_df %>% select(-all_of(cols_to_remove))
          
          output$genes.tsv = rbind(output$cohort_files$genes.tsv, 
                                         cohort$data$genes_scores)
          output$customGenes.tsv = rbind(output$cohort_files$customGenes.tsv,
                                               cohort$data$genes_scores %>% filter(gene %in% cohort$custom_genes))
          output$variants.tsv = rbind(output$cohort_files$variants.tsv,
                                            cohort$data$variants_df %>% filter(gene %in% cohort$data$genes_df$gene))
  
          comphet_PASS_vars <- cohort$data$comphet_df %>% filter(gene %in% cohort$data$genes_df$gene)
          if (nrow(comphet_PASS_vars) > 0) {
            #message("There are comphet to save")
            comphet_to_save <- comphet_PASS_vars %>% 
              gather(key="Variant",value = "varID", v1:v2) %>% 
              inner_join(., cohort$data$variants_df, by=c("varID"="rec_id"))
            
            comphet_to_save$caseID <- sampleID
            
            output$comphet.tsv = rbind(output$cohort_files$comphet.tsv, comphet_to_save)
          }
          return(output)
        }
        stopImplicitCluster()
        RV$cohort_files <- c(RV$cohort_files, sample_results)
        RV$cohort_exit_status <- paste0("Finished with\n\t", nrow(RV$cohort_files$genes.tsv), " genes\n\t", nrow(RV$cohort_files$variants.tsv), " variants")
        callModule(downloadObj, id = "save_cohort",
                   output_prefix= "Cohort",
                   output_data = RV$cohort_files,
                   zip_archive ="Cohort.results.zip")
      })
    })
    
    output$cohort_exit_status <- renderText(RV$cohort_exit_status)
    
    output$cohort_genes_df <- DT::renderDataTable({
      na_values <- base::unique(unlist(app_settings$fill_na$fill_na_genes))
      #Check if there are saved vars, if yes change row backgrounds accordingly
      if (is.null(RV$saved_vars)) {
        datatable(RV$cohort_files$genes.tsv, selection="single") %>% formatStyle(names(RV$cohort_files$genes.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
      } else {
        datatable(RV$cohort_files$genes.tsv, selection="single") %>% formatStyle(names(RV$cohort_files$genes.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values)))) %>%
          formatStyle(  'gene',
                        target = 'row',
                        backgroundColor = styleEqual(saved_genes()$all, rep('yellow', length(saved_genes()$all)))) %>%
          formatStyle(  'gene',
                        target = 'row',
                        backgroundColor = styleEqual(saved_genes()$thisCase, rep('red', length(saved_genes()$thisCase))))
      }
    })
    
    output$cohort_vars_df <- DT::renderDataTable({
      na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
      na_values <- na_values[na_values != 0]
      if (is.null(RV$saved_vars)) {
        datatable(RV$cohort_files$variants.tsv, 
                  selection="multiple",
                  options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
          formatStyle(names(RV$cohort_files$variants.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
      } else {
        datatable(RV$cohort_files$variants.tsv, 
                  selection="multiple",
                  options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
          formatStyle(  'internal_id',
                        target = 'row',
                        backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(base::unique(saved_vars()$all))))) %>%
          formatStyle(  'rec_id',
                        target = 'row',
                        backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
          formatStyle(names(RV$cohort_files$variants.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        
      }
    })
    
    output$cohort_comphet_df <- DT::renderDataTable({
      na_values <- base::unique(unlist(app_settings$fill_na$fill_na_vars))
      na_values <- na_values[na_values != 0]
      if (is.null(RV$saved_vars)) {
        datatable(RV$cohort_files$comphet.tsv, 
                  selection="multiple",
                  options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
          formatStyle(names(RV$cohort_files$comphet.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
      } else {
        datatable(RV$cohort_files$comphet.tsv, 
                  selection="multiple",
                  options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200))) %>%
          formatStyle(  'internal_id',
                        target = 'row',
                        backgroundColor = styleEqual(saved_vars()$all, rep('yellow', length(base::unique(saved_vars()$all))))) %>%
          formatStyle(  'rec_id',
                        target = 'row',
                        backgroundColor = styleEqual(saved_vars()$thisCase, rep('red', length(saved_vars()$thisCase)))) %>%
          formatStyle(names(RV$cohort_files$comphet.tsv), backgroundColor = styleEqual(na_values, rep('gray', length(na_values))))
        
      }
    })

  ## CLOSING OPERATIONS ---------------- 
    onStop(function() {
      df <- isolate(RV$saved_vars)
      if (!is.null(df)) {
        write.table(df, file=paste0(user_dir,"/Saved_variants.csv"), sep=",", row.names=F)
      }
    })   
}

## Run the application ------------------
shinyApp(ui = ui, server = server)

# VARIANT EXPLORER
# Author: Edoardo Giacopuzzi
# Explore and filter annotated variants from V2

# Input are encrypted RData objects created with Prepare_data.R
# Each object contains data from VARAN V2 and var2reg, ROH data and Exp Hunter data

# TODO Set up configurable filters

#Install needed packages if missing
list.of.packages <- c("cyphr","shiny", "DT", "dplyr", "plotly", "kinship2", "tidyr", "shinydashboard", "gridExtra", "ggplot2", "jsonlite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {install.packages(new.packages)}

#Packages from Bioconductor
if(!("GenomicRanges" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install("GenomicRanges")
}

#Packages from github
if(!("scattermore" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
    devtools::install_github('exaexa/scattermore', upgrade = FALSE)
}
if(!("scattermore" %in% installed.packages()[,"Package"])) {
    if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
    devtools::install_github('daattali/shinycssloaders',upgrade = FALSE)
}

#options(repos = BiocManager::repositories())
#getOption("repos")
library(cyphr)
library(shiny)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(kinship2)
library(tidyr)
library(shinydashboard)
library(GenomicRanges)
library(gridExtra)
library(jsonlite)
library(scattermore)
library(shinycssloaders)
source("plotModule.R")
source("downloadModule.R")
source("segregationModule.R")
source("intersectBedModule.R")
source("GQModule.R")

APP_VERSION <- "1.0.5"
resource_dir <- "Resources"
PanelApp_dir <- paste0(resource_dir, "/PanelApp")
GeneLists_dir <- paste0(resource_dir, "/geneLists")
Coverage_dir <- paste0(resource_dir, "/coverage")
BAM_dir <- "/well/gel/HICF2/HICF2_hg38_remap/RareDisease_data/BAM/"
VCF_dir <- "/well/gel/HICF2/HICF2_hg38_remap/RareDisease_data/VCF/"
SV_file <- "/well/gel/HICF2/HICF2_hg38_remap/RareDisease_data/CNV/HICF2_RareDisease_SV.PASS.vcf.gz"

#Set data dir containing variants tables
data_dir <- "encrypted_data"

#################
### FUNCTIONS ###
#################

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

makeBoolean <- function(myvalues,myvector, logic="include") {
    if (logic=="include") { 
        true_value <- 1
    }
    if (logic=="exclude") { 
        true_value <- 0
    }
    
    outlist <- list()
    for (v in myvalues) {
        if (v %in% myvector) {
            outlist[[v]] = true_value
        } else {
            outlist[[v]] = c(0,1)
        }
    }
    return(outlist)
}

`%nin%` = Negate(`%in%`)

#############################
### Load supporting files ###
#############################

#Temporary ped and disease data are loaded from an Robject
#load("Peds_and_HPOs.RData")

##Load hg38 chrom sizes
chr_sizes_file <- getResource("hg38_chromSizes.txt")
chr_sizes <- read.table(chr_sizes_file,header=F, sep="\t", stringsAsFactors = F)
genome_size <- sum(as.numeric(chr_sizes$V2))

##Load gene names and IDs
genes_info_file <- getResource("/hgnc_complete_set.txt")
genes_info <- read.table(genes_info_file, header=T, sep="\t", stringsAsFactors = F)

##Load genes coordinates
genes_bed_file <- gzfile(getResource("gencode.v33.basic.genes.bed.gz"))
genes_bed <- read.table(genes_bed_file, header=F, sep="\t",stringsAsFactors = F)

##Load pathways and GO groups (suppose .gmt file from MSigDB)
gene_anno <- list()
gene_anno[["pathways"]] <- readGMT(getResource("c2.cp.v7.0.symbols.gmt"))
gene_anno[["GO_BP"]] <- readGMT(getResource("c5.bp.v7.0.symbols.gmt"))
gene_anno[["GO_MF"]] <- readGMT(getResource("c5.mf.v7.0.symbols.gmt"))
gene_anno[["GO_CC"]] <- readGMT(getResource("c5.cc.v7.0.symbols.gmt"))

## Load GTeX median expression
GTeX_file <- gzfile(getResource("GTEx_v8_median_TPM.gct.gz"))
GTeX_data <- read.table(GTeX_file,sep="\t",header=T, stringsAsFactors = F)

## Load PanelApp index table
#PanelApp index table contains an id column with panel id number
#Each panel file is then names GRCh38_Panel_[id].bed
PanelApp_file <- getResource("PanelApp_index.tsv")
PanelApp_data <- read.table(PanelApp_file, sep="\t", header = T, stringsAsFactors = F)

PanelApp_genes <- data.frame(entity_name=character(), confidence_level=character(), panel_idx=character(), stringsAsFactors = F)

for (id in PanelApp_data$id) {
    filename <- paste0(PanelApp_dir,"/GRCh38_Panel_", id, ".bed")
    mydf <- read.table(filename, sep="\t", header=T, stringsAsFactors = F, comment.char = "@")
    mydf$panel_idx <- id
    mydf <- mydf[,c("entity_name","confidence_level","panel_idx")]
    PanelApp_genes <- rbind(PanelApp_genes, mydf)
}

## Load additional genes list
#Expect a GeneLists_index.tsv file in the genelist folder
#This file must contain an id column and this ids mus match lists file name
#Example: list1.tsv, id=list1
geneLists_file <- getResource("GeneLists_index.tsv")
geneLists_data <- read.table(geneLists_file, sep="\t", header = T, stringsAsFactors = F)

geneLists_genes <- data.frame(entity_name=character(), genelist_idx=character(), stringsAsFactors = F)

for (id in geneLists_data$id) {
    filename <- paste0(GeneLists_dir,"/", id, ".tsv")
    mydf <- data.frame(entity_name=scan(filename, what="", sep="\n"), genelist_idx = id, stringsAsFactors = F)
    geneLists_genes <- rbind(geneLists_genes, mydf)
}

## Load ClinVar pathogenic genes
# Tab-separated file with gene symbol and list of diseases
ClinVar_file <- getResource("clinvar_path_likelypath_20190704.tsv")
ClinVar_genes <- read.table(ClinVar_file, sep="\t", header = T, stringsAsFactors = F)

##Load GADO distribution
GADO_file <- gzfile(getResource("GADO_distribution.tsv.gz"))
gado_distribution <- read.table(GADO_file, sep="\t", header=T, stringsAsFactors = F)

#######################
### Set environment ###
#######################

##Set axes options for plots
genes_axes_options <- c(
    "GADO Zscore"= "gado_zscore", 
    "Exomiser Pheno score"= "exomiser_gene_pheno_score",
    "gnomAD pLI"= "pLI_gnomad",
    "GDI phred"= "GDI_phred",
    "RVIS intolerance"= "RVIS",
    "EDS reg space score"= "EDS")
variants_axes_options <- c(
    "Maximum population AF"= "max_pop_af",
    "d score"= "d_score",
    "CADD phred"= "CADD_PhredScore", 
    "DANN score"= "DANN_score",
    "ReMM score" = "ReMM_score",
    "SpliceAI score" = "SpliceAI_SNP_SpliceAI_max",
    "dbscSNV splice score" = "dbscSNV_ada",
    "phyloP100 conservation" = "PhyloP100",
    "REVEL score" = "REVEL_score",
    "MCAP score" = "MCAP_score",
    "LinSight score" = "LinSight"
    )
variants_bar_options <- c(
    "Region type" = "reg_type",
    "Variant consequence" = "consequence",
    "Variant type" = "var_type",
    "PanelApp" = "PanelApp",
    "Chromosome" = "chr")

##Plots formatting styles
format1 <- theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
                 axis.text.y = element_text(size=12))

##Names of columns where filter are applied
filter_cols_gene <- c("pLI_gnomAD","GDI_phred")
filter_cols_vars <- c("d_score","SpliceAI_SpliceAI_max","max_pop_af","consequence","cohort_af")

##consequence groups for various variant categories
reg_vars <- c("enhancer_variant","promoter_variant","bivalent_variant","silencer_variant","insulator_variant") 

##var_type groups for various variant categories
sv_vars <- c("DEL","DUP","INV","DEL:ME")

##classification of reg regions sources
reg_sources <- list(
    "database" = c(
        "ENCODE cCREs" = "BENGI", 
        "FOCS" = "FOCS", 
        "HACER" = "HACER", 
        "FANTOM5" = "FANTOM5", 
        "Ensembl Regulatory Build" = "EnsemblRegBuild",
        "RefSeq Regulatory Build" = "RefSeqRegBuild",
        "VISTA enhancers" = "VISTA",
        "EPD6 curated promoters" = "EPD6"),
    "computational" = c(
        "ENCODE HMM profile" = "ENCODE-HMM",
        "DeepLearning DECRES" = "DECRES",
        "SegWey Encyclopedia" = "SegWey"),
    "experimental" = c(
        "CRISPRi-FlowFISH" = "FulcoEtAl2019",
        "CRISPR-Perturb" = "GasperiniEtAl2019",
        "Hi-C screening" = "JungEtAl2019")
)

##Set segregation columns names
segregation_cols <- c(
    "het_affected" = "het_aff", 
    "het_unaffected" = "het_unaff", 
    "hom_affected" = "hom_aff", 
    "hom_unaffected" = "hom_unaff", 
    "comphet_affected" = "comphet_aff" )

##Set reactive objects
RV <- reactiveValues(
    notifications = list(),
    data = 0,
    custom_genes = character(),
    customBed_ranges = FALSE,
    filters_summ_genes = data.frame(),
    filters_summ_vars = data.frame(),
    selected_vars_region = "NONE",
    accepted_reg_db = NULL,
    cov_plot = NULL,
    selected_gene = FALSE,
    messages = list (jigv = messageItem(from="Variant Explorer",
                                    message = "HOW TO USE JIGV SCRIPT
                                    Copy the downloaded script to humbug and launch it.
                                    bash downloaded_script.sh
                                    Leaving the humbug session live, follow the instruction in the link to set up putty
                                    Then you should be able to see IGV on your browser at localhost:5001
                                    Configuration for putty
                                    Source port: 5001
                                    Destination: localhost:5001", 
                                    href="https://www.skyverge.com/blog/how-to-set-up-an-ssh-tunnel-with-putty/")
                      )
    )

############################################
### Read variants data and filter config ###
############################################

#and look for files in data_dir
files <- list.files(data_dir, pattern = ".RData")
samplesID <- gsub("\\.RData[.enc]*","",files, perl = T)

filter_definitions <- read_json("Filters_definitions.json")

######################
### USER INTERFACE ###
######################

ui <- dashboardPage(
    dashboardHeader(
        title = paste0("VarExplorer v", APP_VERSION),
        dropdownMenuOutput("MessageMenu"),
        dropdownMenuOutput("NotificationMenu")
    ),

    dashboardSidebar(
        #drop-down list of cases
        selectInput("CaseCode", h3("Case code:"), choices = sort(unique(samplesID))),
        #h3("Password:"),
        textInput("pwd","Password:",placeholder = "Enter decryption password"),
        actionButton("decrypt_button",label = "Load data"),
        
        #Custom gene list
        fileInput(inputId = "custom_file",label = "Custom gene list:", multiple=FALSE, accept="text/plain", placeholder = "gene list txt file"),
        #Custom genomic regions
        fileInput(inputId = "custom_bed",label = "Region BED:", multiple=FALSE, accept=".bed", placeholder = "region bed file"),
        
        sidebarMenu(
            menuItem("Variants overview", tabName = "overview", icon = icon("th")),
            menuItem("Filters settings", tabName = "filters", icon = icon("th")),
            menuItem("Filter explorer", tabName = "filter_explorer", icon = icon("th")),
            menuItem("Filter results", tabName = "filter_results", icon = icon("th")),
            menuItem("PanelApp and gene lists", tabName = "gene_lists", icon = icon("th")),
            menuItem("Gene details", tabName = "gene_details", icon = icon("th")),
            menuItem("Expansion Hunter", tabName= "expansion_hunter", icon = icon("th")),
            menuItem("Explore coverage", tabName= "coverage_explorer", icon = icon("th"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "overview",
                    fluidRow(
                        box(title = "Pedigree", width=6, status = "primary", solidHeader = TRUE,
                            withSpinner(plotOutput("ped"))
                        ),
                        box(title = "Variant types", width=6, status = "primary", solidHeader = TRUE,
                            withSpinner(plotOutput("Var_type_plot"))
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
            tabItem(tabName = "filters",
                    column(12, align="center", actionButton(inputId = "Apply_filters", label = "Apply filters")),
                    fluidRow(hr()),
                    
                    fluidRow(
                        box(title = "Variant filters", width = 12, status = "primary", solidHeader = TRUE,
                        fluidRow(
                            column(4,
                                h4("General filters"),
                                uiOutput("d_score"),
                                uiOutput("Max_pop_AF"),
                                uiOutput("Cohort_AF"),
                                uiOutput("var_consequence"),
                                uiOutput("custom_bed_check"),
                                checkboxGroupInput("vars_anno_regions","Exclude from:", 
                                                   choices = c("Segmental duplications" = "SegDup", "Low complexity regions" = "LowComplexity", "Highly variable genes" = "TopVariableGenes"),
                                                   inline = FALSE) ),
                            column(4,
                                h4("Missense filters"),
                                uiOutput("CADD"),
                                uiOutput("DANN"),
                                uiOutput("MCAP"),
                                uiOutput("REVEL") ),
                            column(4,
                                h4("Splice filters"),   
                                uiOutput("spliceAI")) ),
                        fluidRow(
                            column(4, 
                                h4("Regulatory filters"),
                                uiOutput("DB_sources"),
                                selectInput("reg_connected_gene", "Gene connection:", choices = c("ALL" = "ALL", "Closest gene" = "closest_gene", "From database" = "reg_db"), multiple = TRUE, selected="ALL"),
                                uiOutput("LinSight"),
                                uiOutput("ReMM"),
                                uiOutput("PhyloP100"),
                                uiOutput("LoF_tolerance"),
                                checkboxGroupInput("NC_anno_regions","Included in:", 
                                                   choices = c("TFBS" = "TFBS", "DNase peak" = "DNase", "Ultra-conserved element" = "UCNE"),
                                                   inline = TRUE) ),
                            column(4,
                                h4("Compound het filter"),
                                "At least one of the variant must have the selected consequence",
                                uiOutput("comphet_consequence") ) )
                        )
                    ),

                    fluidRow(
                        box(title = "Gene filters", id = "gene_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        uiOutput("pLI"),
                        uiOutput("GDI"),
                        uiOutput("RVIS"),
                        uiOutput("EDS")),
            
                        box(title = "Segregation filters", id = "segregation_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        fluidRow(
                            column(6, textOutput("Total_affected")),
                            column(6, textOutput("Total_unaffected")) ),
                        uiOutput("segregation_controls"),
                        "Number of affected / unaffected carrying the variant with the given genotype",
                        "Homozygous and heterozygous conditions evaluated as AND, comphet as OR",
                        uiOutput("GQfilter_controls")
                        ) 
                    ),
                    fluidRow(box(title = "ROH regions filter", id = "ROH_filters_box", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE,
                        uiOutput("ROH_sample_select"),
                        uiOutput("ROH_filters_UI") )
                    )
            ),
            tabItem(tabName = "filter_explorer",
                    h3("Summary of filters effect"),
                    fluidRow(column(6,withSpinner(plotlyOutput("summary_variants_filters", width = "100%"))), column(6,withSpinner(plotlyOutput("summary_genes_filters", width = "100%")))),
                    
                    h3("Genes filter evaluation"),
                    withSpinner(plotSelectedUI("genes_scatter", variables=list("x"=genes_axes_options,"y"=genes_axes_options), plotly=TRUE)),
                    br(),
                    
                    h3("Variants filter evaluation"),
                    withSpinner(plotSelectedUI("variants_scatter", variables=list("x"=variants_axes_options,"y"=variants_axes_options), set_limits=c("x","y"), plotly=FALSE)),
                    br(),
                    withSpinner(plotSelectedUI("variants_barplot", variables=list("x"=variants_bar_options), plotly=FALSE))    
            ),
            tabItem(tabName = "filter_results",
                    fluidRow(
                        column(8, withSpinner(plotlyOutput("GADO_rank"))),
                        column(4, withSpinner(plotlyOutput("GADO_distribution")))
                    ),
                    DT::dataTableOutput("genesTable"),
                    h3("Custom list genes"),
                    DT::dataTableOutput("customGenesTable"),
                    hr(),
                    fluidRow(column(8), 
                             column(3, offset=1,downloadObjUI(id = "save_results", label = "Download results")))    
            ),
            tabItem(tabName = "gene_lists",
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

            ),
            tabItem(tabName = "gene_details",
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
                    fluidRow(br()),
                    fluidRow(
                        column(4, downloadObjUI("get_jigv_script", label = "Download JIGV script")),
                        column(4, textOutput("igv_region"))
                    ),
                    fluidRow(hr()),
                    box(title = "PanelApp and ClinVar", id = "panelapp_clinvar_details", status = "info", solidHeader = TRUE, width = 12,
                        collapsible = TRUE, collapsed = TRUE,
                        h3("PanelApp panels"),
                        DT::dataTableOutput("panelapp_detail_tab"),
                        h3("ClinVar associated diseases"),
                        DT::dataTableOutput("clinvar_detail_tab"),
                        h3("Interesting gene lists"),
                        DT::dataTableOutput("genelists_detail_tab") ),
                    box(title = "Pathways and ontology", id = "path_go_details", status = "info", solidHeader = TRUE, width = 12,
                        collapsible = TRUE, collapsed = TRUE,
                        div(DT::dataTableOutput("geneInfo"),style="font-size: 75%") ),
                    box(title = "GTeX median expression", id = "gtex_expression", status = "info", solidHeader = TRUE, width=12,
                        collapsible = TRUE, collapsed = TRUE,
                        plotOutput("gtex_plot") )
                    
            ),
            tabItem(tabName = "expansion_hunter",
                    uiOutput("exphunter_selection"),
                    DT::dataTableOutput("exphunter_variants"),
                    plotOutput("exphunter_plot", width = "100%", height = "800px")
            ),
            tabItem(tabName = "coverage_explorer",
                    fluidRow(
                        column(3,selectInput("chr_coverage", "Chromosome: ", choices = paste("chr", c(1:22,"X","Y","M"), sep=""))),
                        column(3,uiOutput("gene_cov_select")),
                        column(3,textInput("cov_region_pad", label = "Flanking region (kb): ", value = "500")),
                        column(3,actionButton("plot_coverage","Plot coverage")) ),
                    withSpinner(plotlyOutput("coverage_plot"))
            )
        )
    )
)

########################
### SERVER FUNCTIONS ###
########################

server <- function(input, output, session) {
    
    ################
    ### Load Files 
    ################
    
    observeEvent(input$decrypt_button, {
        #Rememeber that df in pre-processed object are generated with fread so slicing works differently [,..indexes]
        
        #Reset notifications, messages, custom bed and custom genes
        RV$messages <- list()
        RV$custom_genes <- character()
        RV$customBed_ranges <- FALSE
        RV$notifications <- list()
        RV$cov_plot <- NULL
        RV$selected_gene <- FALSE
        
        
        #Load data from plain or encrypted object
        if (file_test("-f", paste0(data_dir,"/",input$CaseCode,".RData.enc"))) {
            RV$data <- decrypt_datafile(paste0(data_dir,"/",input$CaseCode,".RData.enc"), pwd = input$pwd)
        } else if (file_test("-f", paste0(data_dir,"/",input$CaseCode,".RData"))) {
            RV$data <- readRDS(paste0(data_dir,"/",input$CaseCode,".RData"))
        }
        if (inherits(RV$data, "list")) {
            RV$data$variants_df <- RV$data$variants_df %>% replace_na(filter_definitions$fill_na_vars)
            RV$data$genes_scores <- RV$data$genes_scores %>% replace_na(filter_definitions$fill_na_genes)
            RV$data$ROH_data$ROHClass <- cut(RV$data$ROH_data$Length_bp, 
                                             breaks = c(0,500000,2000000,max(RV$data$ROH_data$Length_bp)), 
                                             labels = c("small (< 500kb)","medium (500kb-2Mb)","large (>= 2Mb)"))
            RV$notifications[["decrypt"]] <- notificationItem(
                text = paste0("Loaded data for ", input$CaseCode),
                icon = icon("check-circle"),
                status = "success")
            RV$notifications[["pedigree"]] <- notificationItem(
                text = paste0("Individuals: ", RV$data$
                                  n_all_samples),
                icon = icon("check-circle"),
                status = "success")
            RV$notifications[["variants"]] <- notificationItem(
                text = paste0("variants loaded: ", length(unique(RV$data$variants_df$var_id))),
                icon = icon("check-circle"),
                status = "success")
            RV$notifications[["genes"]] <- notificationItem(
                text = paste0("distinct genes: ", length(unique(RV$data$genes_scores$gene))),
                icon = icon("check-circle"),
                status = "success")
            
            GQ_cols_all <- which(colnames(RV$data$variants_df) %in% paste("GQ", RV$data$all_samples, sep="_"))
            RV$GQ_cols_all <- GQ_cols_all
            RV$GQ_cols_affected <- which(colnames(RV$data$variants_df) %in% paste("GQ", RV$data$affected_samples, sep="_"))
            RV$maxGQ <- max(RV$data$variants_df[,..GQ_cols_all], na.rm = T)
            
        } else {
            RV$notifications[["decrypt"]] <- notificationItem(
                text = "Error loading data! Missing file or wrong password!",
                icon = icon("exclamation-circle"),
                status = "danger")
        }
    })
    
    observeEvent(input$custom_file, {
        req(input$custom_file)
        tryCatch({
            RV$custom_genes <- scan(input$custom_file$datapath,what="",sep="\n")
            RV$notifications[["custom_file"]] <- notificationItem(
                text = paste0(length(RV$custom_genes), " genes loaded"),
                icon = icon("check-circle"),
                status = "success")
        }, error=function(cond) {
            RV$notifications[["custom_file"]] <- notificationItem(
                text = paste0("Failed loading custom genes list"),
                icon = icon("exclamation-circle"),
                status = "danger")
        })
    })

    observeEvent(input$custom_bed, {
        req(input$custom_bed)
        message(input$custom_bed$datapath)
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
    
    output$Disease <- renderText({
        paste0(HPOs$Disease[HPOs$X.CaseID == input$CaseCode])
    })
    
    output$NotificationMenu <- renderMenu({
        dropdownMenu(type = "notifications", .list = RV$notifications)
    })
    
    output$MessageMenu <- renderMenu({
        dropdownMenu(type = "messages", .list = RV$messages)
    })
    
    #######################################
    ### Apply filters on button click
    #######################################
        
    observeEvent(input$Apply_filters, {
        if ("ALL" %in% input$consequence) { 
            RV$accepted_consequence <- sort(unique(RV$data$variants_df$consequence))    
        } else {
            RV$accepted_consequence <- input$consequence
        }
        
        if ("ALL" %in% input$comphet_consequence_filter) { 
            RV$comphet_consequence <- sort(unique(RV$data$variants_df$consequence))    
        } else {
            RV$comphet_consequence <- input$comphet_consequence_filter
        }
        
        if ("ALL" %in% input$reg_db_select) { 
            RV$accepted_reg_db <- c(RV$accepted_reg_db, reg_sources$database)    
        } else {
            RV$accepted_reg_db <- c(RV$accepted_reg_db, input$reg_db_select)
        }
        
        if ("ALL" %in% input$reg_comp_select) { 
            RV$accepted_reg_db <- c(RV$accepted_reg_db, reg_sources$computational)     
        } else {
            RV$accepted_reg_db <- c(RV$accepted_reg_db, input$reg_comp_select)
        }
        
        if ("ALL" %in% input$reg_exp_select) { 
            RV$accepted_reg_db <- c(RV$accepted_reg_db, reg_sources$experimental)     
        } else {
            RV$accepted_reg_db <- c(RV$accepted_reg_db, input$reg_exp_select)
        }
        
        if ("ALL" %in% input$reg_connected_gene) { 
            RV$accepted_connected_gene <- c("closest_gene", "reg_db")     
        } else {
            RV$accepted_connected_gene <- input$reg_connected_gene
        }
    })    

    #######################################
    ### Data reactive tables configuration
    #######################################
    
    segregating_vars <- reactive({
        input$Apply_filters
        shiny::validate(need(inherits(RV$data, "list"), FALSE))
        isolate(
            callModule(segregationModule, "segregation", segregation_df = RV$data$segregation_df, cols_names = segregation_cols)
        )
    })
    
    variants_df <- reactive({
        input$Apply_filters
        shiny::validate(need(inherits(RV$data, "list"), FALSE))
        
        isolate({
            
            #Get rec_id for reg variants in accepted db sources
            accepted_reg_recid <- RV$data$variants_df$rec_id[grep(paste(RV$accepted_reg_db, collapse="|"), RV$data$variants_df$db_source)]
           
            #Read which regions are selected from checkboxes and convert to 0/1 values
            NC_reg_anno <- makeBoolean(c("TFBS","DNase","UCNE"), input$NC_anno_regions, logic = "include")
            var_reg_anno <- makeBoolean(c("SegDup","LowComplexity","TopVariableGenes"), input$vars_anno_regions, logic = "exclude")
            
            ##ROH FILTER, BED FILTER and GQ FILTER
            #These filters are high-level and applied before all others
            #call the bed regions module and get rec_id for variants in the ROH regions
            if (!is.null(RV$data$ROH_ranges)) {
                #message("call bed module")
                vars_pass_ROH <- callModule(bedfilterModule,"ROH_filter",variants_ranges=RV$data$variants_ranges, bed_ranges=RV$data$ROH_ranges[[input$ROH_sample]], multiplier=1000)
            } else {
                vars_pass_ROH <- RV$data$variants_df$rec_id
            }
            
            #call the bed regions module and get rec_id for variants in the custom BED regions
            if (inherits(RV$customBed_ranges, "GRanges")) {
                vars_pass_BED <- callModule(bedfilterModule,"bed_filter",variants_ranges=RV$data$variants_ranges, bed_ranges=RV$customBed_ranges)
            } else {
                vars_pass_BED <- RV$data$variants_df$rec_id
            }
            
            #Get rec_id for variants passing the GQ filter
            vars_pass_GQ <- callModule(GQfilterModule, "GQ_filter", 
                                       variants_df = RV$data$variants_df, 
                                       GQ_cols = RV$GQ_cols_all, 
                                       affected_cols = RV$GQ_cols_affected,
                                       exclude_var_type = sv_vars)
            
            #Pre-filter variants_df and retain only vars passing ROH and GQ filters
            vars_pass <- unique(intersect(intersect(vars_pass_ROH, vars_pass_GQ), vars_pass_BED))
            prefilter_vars <- as.data.frame(RV$data$variants_df %>% filter(rec_id %in% vars_pass))
            
            ## FILTER VARS  
            #first select variants that pass the filters
            pass_vars <- as.data.frame(prefilter_vars %>% 
                filter(! (
                    d_score < input$d_score_filter | 
                    consequence %nin% RV$accepted_consequence |
                    max_pop_af > input$MaxPopAF_filter |
                    cohort_af > input$CohortAF_filter |
                    SegDup %nin% var_reg_anno$SegDup |
                    LowComplexity %nin% var_reg_anno$LowComplexity |
                    TopVariableGenes %nin% var_reg_anno$TopVariableGenes |
                    (reg_type == "splicing" & 
                         (SpliceAI_SNP_SpliceAI_max < input$spliceAI_filter & 
                              SpliceAI_INDEL_SpliceAI_max < input$spliceAI_filter) ) |
                    (consequence == "missense_variant" &
                         (CADD_PhredScore < input$CADD_filter & 
                              REVEL_score < input$REVEL_filter &
                              MCAP_score < input$MCAP_filter &
                              DANN_score < input$DANN_filter)) |
                    (consequence %in% reg_vars &
                         (ReMM_score < input$ReMM_filter |
                              LinSight < input$LinSight_filter |
                              PhyloP100 < input$PhyloP100_filter |
                              LoF_tolerance > input$LoFtolerance_filter) &
                         rec_id %nin% accepted_reg_recid &
                         reg_type %nin% RV$accepted_connected_gene &
                         TFBS %nin% NC_reg_anno$TFBS &
                         DNase %nin% NC_reg_anno$DNase &
                         UCNE %nin% NC_reg_anno$UCNE)
                )))
            
            #COMPHET FILTERING
            #Select comphet where both vars pass variants filters
            comphet_accepted <- RV$data$comphet_df %>% filter(
                    v1 %in% pass_vars$rec_id &
                    v2 %in% pass_vars$rec_id)
            
            #COMPHET specific consequence
            #This filter allow to select comphet where at least 1 var as the given consequence
            #From pass vars select only the ones with accepted consequence
            #based on the comphet consequence filter
            comphet_accepted_vars <- as.data.frame(pass_vars %>%
                filter(consequence %in% RV$comphet_consequence) )
            
            #Get the final list of vars in accepted comphet
            #from comphet select combo when: 
            #both vars in the combo passed the general variants filters (comphet_accepted df)
            #at least one var in the combo pass the comphet consequence filter
            #the combo passed segregation filters
            comphet_single_vars <- comphet_accepted %>% filter(
                rec_id %in% segregating_vars() &
                (v1 %in% comphet_accepted_vars$rec_id |
                v2 %in% comphet_accepted_vars$rec_id) ) %>% gather(key="Variant",value="VarID",v1:v2) %>% select(VarID)
            
            #GET THE FINAL LIST OF ACCEPTED VARS
            #the final list of segregating accepted vars is now equal to
            #single vars passing filters and segregation
            #vars part of a comphet passing all filters + segregation
            accepted_vars_list <- intersect(segregating_vars(), pass_vars$rec_id)
            accepted_vars_list <- unique(c(accepted_vars_list, comphet_single_vars$VarID))
            
            #return variants dataframe with updated Class column (PASS/FILTER)
            as.data.frame(RV$data$variants_df %>% mutate(Class = ifelse(
                rec_id %in% accepted_vars_list,
                "PASS","FILTER") ) )
        })
    })
    
    comphet_df <- reactive({
        filtered_vars_list <- variants_df()$rec_id[variants_df()$Class == "PASS"]
        
        as.data.frame(RV$data$comphet_df %>% mutate(Class = ifelse(
            v1 %in% filtered_vars_list & 
            v2 %in% filtered_vars_list &
            rec_id %in% segregating_vars(), "PASS", "FILTER")))
    })
    
    genes_df <- reactive({
        RV$filtered_vars_list <- unique(c(
            variants_df()$rec_id[variants_df()$Class == "PASS"], 
            comphet_df()$rec_id[comphet_df()$Class == "PASS"])) 
        
        #Filters on gene scores are applied
        genes_above_score <- as.data.frame(RV$data$genes_scores %>% filter(
            pLI_gnomad >= input$pLI_filter & 
            GDI_phred <= input$GDI_filter &
            RVIS <= input$RVIS_filter &
            EDS >= input$EDS_filter))
            
        as.data.frame(RV$data$genes_df %>% mutate(Class = ifelse(
            variants %in% RV$filtered_vars_list & 
            gene %in% genes_above_score$gene, "PASS", "FILTER")))    
    })
    
    genes_scores <- reactive({
        RV$filtered_genes_list <- unique(genes_df()$gene[genes_df()$Class == "PASS"]) 
        as.data.frame(RV$data$genes_scores %>% mutate(Class = ifelse(
            gene %in% RV$filtered_genes_list, "PASS", "FILTER")))
    })
    
    filters_summ_genes <- reactive ({
        tot_genes <- nrow(genes_scores())
        
        PASS_counts <- c(
            genes_scores() %>% filter(pLI_gnomad >= input$pLI_filter) %>% nrow(),
            genes_scores() %>% filter(GDI_phred <= input$GDI_filter) %>% nrow(),
            genes_scores() %>% filter(RVIS <= input$RVIS_filter) %>% nrow(),
            genes_scores() %>% filter(EDS >= input$EDS_filter) %>% nrow() )
    
        filters_summ_genes <- data.frame(
            Filter=c("pLI gnomAD", "GDI phred", "RVIS", "EDS"), 
            PASS=PASS_counts, 
            FILTERED=(tot_genes-PASS_counts) )
        filters_summ_genes <- gather(filters_summ_genes, key="Class", value="Count", PASS:FILTERED)   
    })
    
    filters_summ_vars <- reactive ({
        tot_vars <- RV$data$variants_df %>% select(var_id) %>% distinct() %>% nrow()
        PASS_counts <- c(
            variants_df() %>% filter(d_score >= input$d_score_filter) %>% select(var_id) %>% distinct() %>% nrow(),
            variants_df() %>% filter(max_pop_af <= input$MaxPopAF_filter) %>% select(var_id) %>% distinct() %>% nrow(),
            variants_df() %>% filter(consequence %in% RV$accepted_consequence) %>% select(var_id) %>% distinct() %>% nrow(),
            variants_df() %>% filter(cohort_af <= input$CohortAF_filter) %>% select(var_id) %>% distinct() %>% nrow()
        )
        filters_summ_vars <- data.frame(Filter=c("d score","MaxPop AF","consequence","cohort AF"),
            PASS=PASS_counts, 
            FILTERED=(tot_vars-PASS_counts))
        filters_summ_vars <- gather(filters_summ_vars, key="Class", value="Count", PASS:FILTERED)
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

    
    ##################
    ### Overview Tab
    ##################
    
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
        Total_ROH_bychr <- aggregate(RV$data$ROH_data$Length_bp, list(RV$data$ROH_data$Sample, RV$data$ROH_data$ROHClass, RV$data$ROH_data$Chromosome), FUN=sum)
        Total_ROH_bychr <- merge(Total_ROH_bychr, chr_sizes, by.x="Group.3", by.y="V1")
        Total_ROH_bychr$pctROH <- Total_ROH_bychr$x / Total_ROH_bychr$V2
        ggplot(Total_ROH_bychr[Total_ROH_bychr$Group.1 == input$ROH_plot_sample,], aes(x=Group.3,y=pctROH,fill=Group.2)) + 
            geom_bar(stat="identity", position = position_dodge(0.9)) + 
            geom_hline(yintercept = 0.1, linetype="dashed") + 
            labs(x="chromosome", y="fraction of chromosome within ROH", fill="ROH size", title=paste0("ROH distribution by chromosome for ", input$ROH_plot_sample)) + 
            format1
    })
    
    #################
    ### Filters tab
    #################
    
    output$d_score <- renderUI({
        sliderInput("d_score_filter", "Min d_score:",
                    min = 0,
                    max = max(RV$data$variants_df$d_score, na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    output$custom_bed_check <- renderUI({
        shiny::validate(need(inherits(RV$customBed_ranges, "GRanges"), "No custom BED"))
        bedcontrolUI("bed_filter", label = "Select only vars in custom regions")
    })
    
    output$spliceAI <- renderUI({
        sliderInput("spliceAI_filter", "Min spliceAI score:",
                    min = 0,
                    max = max(RV$data$variants_df$SpliceAI_SNP_SpliceAI_max, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$CADD <- renderUI({
        sliderInput("CADD_filter", "Min CADD phred:",
                    min = 0,
                    max = max(RV$data$variants_df$CADD_PhredScore, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$DB_sources <- renderUI({
        values <- sort(unique(RV$data$variants_df$db_source))
        values_db <- reg_sources[["database"]][reg_sources[["database"]] %in% values]
        values_db <- c("ALL" = "ALL", values_db)
        values_computational <- reg_sources[["computational"]][reg_sources[["computational"]] %in% values]
        values_computational <- c("ALL" = "ALL", values_computational)
        values_experimental <- reg_sources[["experimental"]][reg_sources[["experimental"]] %in% values]
        values_experimental <- c("ALL" = "ALL", values_experimental)
        
        tagList(
            selectInput("reg_db_select", "Database sources:", choices = values_db, multiple=TRUE, selected="ALL"),
            selectInput("reg_comp_select", "Computational sources:", choices = values_computational, multiple=TRUE, selected="ALL"),
            selectInput("reg_exp_select", "Experimental sources:", choices = values_experimental, multiple=TRUE, selected="ALL")
        )
        
    })
    
    output$LinSight <- renderUI({
        sliderInput("LinSight_filter", "Min LinSight score:",
                    min = 0,
                    max = max(RV$data$variants_df$LinSight, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$ReMM <- renderUI({
        sliderInput("ReMM_filter", "Min ReMM score:",
                    min = 0,
                    max = max(RV$data$variants_df$ReMM_score, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$PhyloP100 <- renderUI({
        sliderInput("PhyloP100_filter", "Min PhyloP100 score:",
                    min = min(RV$data$variants_df$PhyloP100, na.rm = T),
                    max = max(RV$data$variants_df$PhyloP100, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$LoF_tolerance <- renderUI({
        sliderInput("LoFtolerance_filter", "Max LoF tolerance score:",
                    min = 0,
                    max = max(RV$data$variants_df$LoF_tolerance, na.rm = T),
                    value = 1,
                    step = 0.02)        
    })
    
    output$REVEL <- renderUI({
        sliderInput("REVEL_filter", "Min REVEL score:",
                    min = 0,
                    max = max(RV$data$variants_df$REVEL_score, na.rm = T),
                    value = 0,
                    step = 0.01)
    })
    
    output$DANN <- renderUI({
        sliderInput("DANN_filter", "Min DANN score:",
                    min = 0,
                    max = max(RV$data$variants_df$DANN_score, na.rm = T),
                    value = 0,
                    step = 0.01)
    })
    
    output$MCAP <- renderUI({
        sliderInput("MCAP_filter", "Min M-CAP score:",
                    min = 0,
                    max = max(RV$data$variants_df$MCAP_score, na.rm = T),
                    value = 0,
                    step = 0.01)
    })
    
    output$Max_pop_AF <- renderUI({
        sliderInput("MaxPopAF_filter", "Max population AF:",
                    min = 0,
                    max = max(RV$data$variants_df$max_pop_af, na.rm = T),
                    value = max(RV$data$variants_df$max_pop_af, na.rm = T),
                    step = 0.001)
    })
    
    output$Cohort_AF <- renderUI({
        sliderInput("CohortAF_filter", "Max cohort AF:",
                    min = 0,
                    max = max(RV$data$variants_df$cohort_af, na.rm = T),
                    value = max(RV$data$variants_df$cohort_af, na.rm = T),
                    step = 0.001)
    })
    
    output$pLI <- renderUI({
        sliderInput("pLI_filter", "Min gnomad pLI:",
                    min = 0,
                    max = max(RV$data$genes_scores$pLI_gnomad[RV$data$genes_scores$pLI_gnomad != 99], na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    output$GDI <- renderUI({
        startvalue = 0
        sliderInput("GDI_filter", "Max phred GDI:",
                    min = 0,
                    max = max(RV$data$genes_scores$GDI_phred, na.rm = T),
                    value = max(RV$data$genes_scores$GDI_phred, na.rm = T),
                    step = 0.05)
    })
    
    output$RVIS <- renderUI({
        startvalue = 0
        sliderInput("RVIS_filter", "Max RVIS score:",
                    min = min(RV$data$genes_scores$RVIS[RV$data$genes_scores$RVIS != -99], na.rm = T),
                    max = max(RV$data$genes_scores$RVIS, na.rm = T),
                    value = max(RV$data$genes_scores$RVIS, na.rm = T),
                    step = 0.05)
    })
    
    output$EDS <- renderUI({
        startvalue = 0
        sliderInput("EDS_filter", "Min EDS score:",
                    min = 0,
                    max = max(RV$data$genes_scores$EDS[RV$data$genes_scores$EDS != 99], na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    output$segregation_controls <- renderUI({
        segregationUI("segregation", choices_affected = RV$data$values_affected, choices_unaffected=RV$data$values_unaffected)
    })
    
    output$GQfilter_controls <- renderUI({
        message(RV$maxGQ)
        GQfilterUI("GQ_filter",maxGQ = RV$maxGQ, defaultGQ=10)    
    })
    
    output$var_consequence <- renderUI({
        var_types <- sort(unique(RV$data$variants_df$consequence))
        names(var_types) <- sort(unique(RV$data$variants_df$consequence))
        var_types <- c("ALL" = "ALL", var_types)
        selectInput("consequence", "Variant consequence:", choices = var_types, multiple=TRUE, selected="ALL")
    })
    
    output$comphet_consequence <- renderUI({
        var_types <- sort(unique(RV$data$variants_df$consequence))
        names(var_types) <- sort(unique(RV$data$variants_df$consequence))
        var_types <- c("ALL" = "ALL", var_types)
        selectInput("comphet_consequence_filter", "Variant consequence:", choices = var_types, multiple=TRUE, selected="ALL")
    })
    
    output$ROH_sample_select <- renderUI ({
        selectInput("ROH_sample", label = "Select sample:", choices = names(RV$data$ROH_ranges), multiple = FALSE, selected = "AFFECTED_SHARED")
    })
    
    output$ROH_filters_UI <- renderUI({ 
        shiny::validate(need(!is.null(RV$data$ROH_ranges), "No ROH regions loaded"))
    
        bedcontrolUI("ROH_filter", label = "Select only vars in ROH regions", 
            slider_config = c(
                "label" = "ROH min dimension (kb)",
                "min" = 0,
                "max" = max(RV$data$ROH_ranges[[input$ROH_sample]]$value) / 1000,
                "step" = 10,
                "value" = 250
            ))
    })

    ########################
    ### Filter Explorer tab
    ########################
    
    output$summary_variants_filters <- renderPlotly({
        ggplotly(
            ggplot(filters_summ_vars(), aes(x=Filter,y=Count,fill=Class)) + geom_bar(stat="identity") + labs(y="N variants") + theme(axis.text.x = element_text(angle=45, hjust=1))
        )
    })
    
    output$summary_genes_filters <- renderPlotly({
        ggplotly(
            ggplot(filters_summ_genes(), aes(x=Filter,y=Count,fill=Class)) + geom_bar(stat="identity") + labs(y="N genes") + theme(axis.text.x = element_text(angle=45, hjust=1))
        )
    })
    
    #output$genes_scatter <- renderPlotly({
    #    ggplotly(
    #        ggplot(genes_scores(), aes_string(x=input$genes_X_axis, y=input$genes_Y_axis, color="Class", label="gene")) + geom_point(size=1) 
    #    )
    #})
    
    callModule(plotModule, "genes_scatter", plot_data = genes_scores(), missingValues = c(99,-99), plotType = "scatter", plotOptions = list("size" = 1), variables = list("color"="Class", "label"="gene"))
    callModule(plotModule, "variants_scatter", plot_data = variants_df(), missingValues = c(99,-99), plotType = "bigdata", variables = list("color"="Class", "size" = 1))
    callModule(plotModule, "variants_barplot", plot_data = variants_df(), plotType = "barplot", variables = list("fill"="Class"), additionalOptions = list(format1, scale_y_sqrt()))
    
    ########################
    ### Filter Results tab
    ########################
    
    output$GADO_rank <- renderPlotly({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        if (length(input$genesTable_rows_selected) > 0) {
            Zscore <- genedetail_tab[input$genesTable_rows_selected, "gado_zscore"]
            GADO_plot <- ggplot(genes_scores(), aes(x=gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = Zscore, color="red") + geom_vline(xintercept = RV$data$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        } else {
            GADO_plot <- ggplot(genes_scores(), aes(x=gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = RV$data$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        }
        ggplotly( GADO_plot )
    })
    
    output$GADO_distribution <- renderPlotly({
        shiny::validate(need(input$genesTable_rows_selected, "Select a gene"))
        
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "gene"]
        Zscore <- genedetail_tab[input$genesTable_rows_selected, "gado_zscore"]
        
        GADO_dist <- ggplot(gado_distribution[gado_distribution$Hgnc == gene_name,], aes(x=Zscore)) + 
            geom_histogram(bins = max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name])*10) + 
            geom_vline(xintercept = Zscore, color = "red") +
            scale_x_continuous(breaks=round(seq(min(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name]),max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name]),1)))
        ggplotly( GADO_dist )
     })
    
    candidate_genes_df <- reactive({
        as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% mutate(row_idx = row_number()) %>% select(-Class))
    })
    
    output$genesTable <- DT::renderDataTable(selection="single", {
        candidate_genes_df()
    })
    
    genesTable_proxy <- DT::dataTableProxy("genesTable")
    
    customGenesTable_df <- reactive({
        candidate_genes_df() %>% filter(gene %in% RV$custom_genes)
    })
    
    output$customGenesTable <- DT::renderDataTable(selection="single", {
        customGenesTable_df()
    })
    
    observeEvent(input$customGenesTable_rows_selected, {
        row_idx <- customGenesTable_df()[input$customGenesTable_rows_selected,"row_idx"]
        genesTable_proxy %>% selectRows(row_idx) 
    })
    
    callModule(downloadObj, id = "save_results",
    output_prefix=input$CaseCode,
    output_data = list(
        "genes.tsv" = as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class)),
        "customGenes.tsv" = as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% filter(gene %in% RV$custom_genes) %>% select(-Class)),
        "variants.tsv" = as.data.frame(variants_df() %>% filter(Class == "PASS", gene %in% RV$filtered_genes_list)),
        "comphet.tsv" = as.data.frame(comphet_df() %>% filter(Class == "PASS", gene %in% RV$filtered_genes_list) %>% 
                            gather(key="Variant",value = "varID", v1:v2) %>% 
                            inner_join(., variants_df()[variants_df()$Class == "PASS",], by=c("varID"="rec_id")) %>%
                            select(-Class.x,-Class.y))),
    zip_archive = paste0(input$CaseCode, ".results.zip") )
        
    ################################
    ### PanelApp and genes lists tab
    ################################
    
    output$PanelApp_panels_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
        PanelApp_panels_df()
    })
    
    PanelApp_genes_df <- reactive ({
        shiny::req(input$PanelApp_panels_table_rows_selected)
        panelID <- PanelApp_panels_df()[input$PanelApp_panels_table_rows_selected, "id"]
        genes_df <- as.data.frame(PanelApp_genes %>% filter(panel_idx == panelID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., genes_scores(), by = c("entity_name" = "gene")))
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
        genes_df <- as.data.frame(geneLists_genes %>% filter(genelist_idx == listID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., genes_scores(), by = c("entity_name" = "gene")))
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
    
    ####################
    ### Gene detail tab
    ####################
    
    #genedetail_tab <- reactive ({
    #    candidate_genes_df()
    #})
    
    gene_name <- reactive ({
        symbol = candidate_genes_df()[input$genesTable_rows_selected, "gene"]
        updateSelectInput(session = session, inputId = "chr_coverage",selected = genes_bed$V1[genes_bed$V4 == symbol])
        RV$selected_gene <- symbol
        candidate_genes_df()[input$genesTable_rows_selected, "gene"]
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
        #shiny::validate(gene_name != "")
        
        comphet_details <- comphet_df() %>% filter(Class == "PASS", gene == gene_name()) %>% gather(key="variant",value = "varID", v1:v2) %>% select(-gene, -variant) %>% distinct()
        #shiny::validate(need(nrow(comphet_details)>0, 'No compound het variants'))
        comphet_details <- merge(comphet_details,variants_df()[variants_df()$Class == "PASS",], by.x="varID",by.y="rec_id")
        as.data.frame(comphet_details %>% select(-Class.x,-Class.y,) %>% arrange(rec_id))
    })
    
    gene_vars_df <- reactive ({
        #shiny::validate(gene_name != "")
        #comphet_details <- comphet_df() %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", V1:V2)
        
        #as.data.frame(variants_df() %>% filter(Class == "PASS", Gene == gene_name))
        as.data.frame(variants_df() %>% filter(Class == "PASS", gene == gene_name(), rec_id %nin% gene_comphet_vars_df()$varID)) 
    })
    
    output$variantsTable <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
        #shiny::validate(need(gene_name != "", 'No gene selected'))
        
        as.data.frame(gene_vars_df())
    })
    
    output$comphetTable <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
        #shiny::validate(need(gene_name != "", 'No gene selected'))
        shiny::validate(need(nrow(gene_comphet_vars_df())>0, 'No compound het variants'))

        as.data.frame(gene_comphet_vars_df())
    })
    
    output$panelapp_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE),  {
        panelID <- PanelApp_genes[PanelApp_genes$entity_name == gene_name(), "panel_idx"]
        conf_level <- PanelApp_genes[PanelApp_genes$entity_name == gene_name(), "confidence_level"]
        df <- PanelApp_panels_df()[PanelApp_panels_df()$id %in% panelID,]
        df$gene <- gene_name()
        df$confidence_level <- conf_level
        df[,c("gene","confidence_level","id","name","disease_group","version","relevant_disorders")]
    })
    
    output$clinvar_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        ClinVar_df()[ClinVar_df()$gene == gene_name(),]
    })
    
    output$genelists_detail_tab <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE),  {
        listID <- geneLists_genes[geneLists_genes$entity_name == gene_name(), "genelist_idx"]
        geneLists_df()[geneLists_df()$id == listID,]
    })
    
    
    JIGV_script <- reactive ({
        shiny::validate(need(!is.null(input$variantsTable_rows_selected) | !is.null(input$comphetTable_rows_selected), "No rows selected"))
        chromosome <- unique(c(gene_comphet_vars_df()$chr[input$comphetTable_rows_selected], gene_vars_df()$chr[input$variantsTable_rows_selected]))
        start_pos <- unique(c(gene_comphet_vars_df()$start[input$comphetTable_rows_selected], gene_vars_df()$start[input$variantsTable_rows_selected]))
        end_pos <- unique(c(gene_comphet_vars_df()$end[input$comphetTable_rows_selected], gene_vars_df()$end[input$variantsTable_rows_selected]))
        
        shiny::validate(need(length(chromosome) == 1, "Please select vars on the same chromosome"))
        
        start_pos <- min(start_pos)
        end_pos <- max(end_pos)
        region <- paste0(chromosome,":",start_pos,"-",end_pos)
        out_file <- paste0(chromosome,"_",start_pos,"-",end_pos,".sh")
        RV$selected_vars_region <- region
        
        
        jigv_command <- paste("jigv --region", region,
                              paste('"',BAM_dir, RV$data$affected_samples, '.bam#', RV$data$affected_samples, '_affected"', sep="", collapse=" " ),
                              paste('"',BAM_dir, RV$data$unaffected_samples, '.bam#', RV$data$affected_samples, '_unaffected"', sep="", collapse=" " ),
                              paste0(VCF_dir, RV$data$pedigree, ".PASS.NORM.vcf.gz"), 
                              SV_file, collapse=" ")
        message(out_file)
        message(jigv_command)
        list(outfile=out_file, command=jigv_command)
    })
    
    callModule(downloadObj, id="get_jigv_script", output_prefix= JIGV_script()$outfile, output_data=JIGV_script()$command, col_names=FALSE)

    output$igv_region <- renderText({
        paste0("Generated script: ", JIGV_script()$outfile)
    })
    
    ########################
    ### Expansion Hunter tab
    ########################
    
    output$exphunter_selection <- renderUI({
        selectInput("exphunter_loci", h3("Gene:"), choices = sort(unique(RV$data$ExpHunter$PEDIGREE$LocusId)), multiple = FALSE)
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
    
    ########################
    ### Explore coverage tab
    ########################
    
    output$gene_cov_select <- renderUI({
        genes_list <- unique(genes_bed$V4[genes_bed$V1 == input$chr_coverage])
        if (RV$selected_gene != FALSE) {
            selected_gene = RV$selected_gene
        } else {
            selected_gene = "NO_GENE"
        }
        selectInput("gene_coverage", "Gene: ", choices = c("NO_GENE",genes_list), selected = selected_gene)
    })
    
    observeEvent(input$plot_coverage, {
        filename <- paste0(Coverage_dir,"/indexcov-", input$chr_coverage, ".bed.gz")
        shiny::validate(need(file_test("-f", filename), "No coverage data"))
        
        cov_file <- gzfile(filename)
        cov_data <- read.table(cov_file, sep="\t", header=T, comment.char = "@", stringsAsFactors = F)
        cov_df <- cov_data %>% select(all_of(c("X.chrom","start","end", RV$data$all_samples))) %>% gather(.,key="sample",value="norm_cov", 4:(3+length(RV$data$all_samples)))
        cov_df$middle_point <- round(cov_df$start + ((cov_df$end - cov_df$start) / 2))
        
        if (input$gene_coverage != "NO_GENE") {
            plot_start <- genes_bed$V2[genes_bed$V4 == input$gene_coverage] - as.numeric(input$cov_region_pad)*1000
            plot_end <- genes_bed$V3[genes_bed$V4 == input$gene_coverage] + as.numeric(input$cov_region_pad)*1000
            cov_df <- cov_df %>% filter(middle_point >= plot_start & middle_point <= plot_end)   
        }
        RV$cov_plot <- ggplot(cov_df, aes(x=middle_point/100000, y=norm_cov, color=sample)) + 
            geom_line(size=0.5, alpha=0.5) + scale_y_continuous(limits=c(0,3), breaks=0:3) + 
            theme(axis.text.x=element_text(size=8, angle=45, hjust=1)) + 
            labs(x="position (x 100kb)", title=paste0("Chromosome: ", input$chr_coverage), subtitle = paste0("Gene: ", input$gene_coverage))
        
    })
    
    output$coverage_plot <- renderPlotly({
        req(RV$cov_plot)
        ggplotly(RV$cov_plot, tooltip = c("color","y"), dynamicTicks = T)
    })
    

    
}

# Run the application 
shinyApp(ui = ui, server = server)

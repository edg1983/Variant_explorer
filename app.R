# VARIANT EXPLORER
# Author: Edoardo Giacopuzzi
# Explore and filter annotated variants from V2

# Input are tables from the VARAN V2 and var2reg

# TODO Add separate table reporting all Clinvar pathogenic vars
# TODO Set up configurable filters

library(shiny)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(kinship2)
library(tidyr)
library(shinydashboard)
source("plotModule.R")
source("downloadModule.R")
source("segregationModule.R")

resource_dir <- "Resources"
PanelApp_dir <- paste0(resource_dir, "/PanelApp")
GeneLists_dir <- paste0(resource_dir, "/geneLists")
segregation_cols <- c(
    "het_affected" = "het_aff", 
    "het_unaffected" = "het_unaff", 
    "hom_affected" = "hom_aff", 
    "hom_unaffected" = "hom_unaff", 
    "comphet_affected" = "comphet_aff" )

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

`%nin%` = Negate(`%in%`)

#############################
### Load supporting files ###
#############################

#Temporary ped and disease data are loaded from an Robject
load("Peds_and_HPOs.RData")

##Read ped file
#ped_file <- "~/Servers/well/gel/HICF2/HICF2_hg38_remap/data/PED_sexfix/All_samples.ped"
#ped_df <- read.table(ped_file, header=F, sep="\t", stringsAsFactors = F)
#ped_df$V5[ped_df$V5==1] <- "male"
#ped_df$V5[ped_df$V5==2] <- "female"
#ped_df$V5[ped_df$V5==0] <- "unknown"
#ped_df$V6[ped_df$V6 == 0] <- 1
#all_peds <- with(ped_df, pedigree(dadid = V3, momid = V4, affected = V6, sex = V5, famid =V1, missid = 0,id = V2))

##Load gene names and IDs
genes_info_file <- getResource("/hgnc_complete_set.txt")
genes_info <- read.table(genes_info_file, header=T, sep="\t", stringsAsFactors = F)

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
GADO_file <- gzfile(getResource("GADO_distribution.RData"))
gado_distribution <- read.table(GADO_file, sep="\t", header=T, stringsAsFactors = F)

#######################
### Set environment ###
#######################

##Set axes options for plots
genes_axes_options <- c(
    "GADO Zscore"= "Gado_zscore", 
    "Exomiser Pheno Rank"= "Exomiser_GenePhenoScore",
    "gnomAD pLI"= "pLI_gnomad", 
    "GDI score"= "GDI_phred")
variants_axes_options <- c("Maximum population AF"= "MaxPopAF",
    "d score"= "d_score",
    "CADD phred"= "CADD_PhredScore", 
    "DANN score"= "DANN_DANN",
    "ReMM score" = "ReMM_score",
    "SpliceAI score" = "SpliceAI_SpliceAI_max",
    "dbscSNV splice score" = "dbscSNV_ada")
variants_bar_options <- c(
    "Region type" = "Reg_type",
    "Variant consequence" = "Consequence",
    "Variant type" = "VarType",
    "PanelApp" = "PanelApp",
    "Chromosome" = "Chr")

##Plots formatting styles
format1 <- theme(axis.text.x = element_text(angle=45, hjust=1))

##Names of columns where filter are applied
filter_cols_gene <- c("pLI_gnomAD","GDI_phred")
filter_cols_vars <- c("d_score","SpliceAI_SpliceAI_max","MaxPopAF","Consequence","cohortAF")

##Consequence for regulatory vars
reg_vars <- c("enhancer_variant","promoter_variant","bivalent_variant","silencer_variant") 

##Set reactive objects
RV <- reactiveValues(
    variant_df = data.frame(), 
    genes_scores = data.frame(), 
    genes_df = data.frame(),
    comphet_df = data.frame(),
    load_status = character(),
    custom_genes = character(),
    filters_summ_genes = data.frame(),
    filters_summ_vars = data.frame(),
    values_affected = 0,
    values_unaffected = 0)

##########################
### Read variants data ###
##########################

#Set data dir containing variants tables and look for files
data_dir <- "example_data/new_tables"

files <- list.files(data_dir, pattern = "*.vars.tsv.gz")
samplesID <- gsub("V2\\.|\\.var2reg.vars.tsv.gz","",files)

#We can switch to read from a input table containing sampleID and ped locations...
#samples_df <- data.frame(ID=samplesID,file=files)

######################
### USER INTERFACE ###
######################

ui <- dashboardPage(
    dashboardHeader(
        title = "Variant Explorer",
        dropdownMenu(type = "notifications",
                     notificationItem(
                         text = "Beta version, functionality limited and bugs expected",
                         icon = icon("exclamation-triangle"),
                         status = "warning"
                     )
        )
    ),

    dashboardSidebar(
        #drop-down list of cases
        selectInput("CaseCode", h3("Case code:"), choices = sort(unique(samplesID))),
        h3("Disease:"),
        textOutput("Disease"),
        
        h3("Custom gene list"),
        fileInput(inputId = "custom_file",label = "Custom gene list:", multiple=FALSE, accept="text/plain", placeholder = "gene list txt file"),
        #actionButton(inputId = "Load_file", label = "Load file"),
        textOutput("Loading_result"),
        sidebarMenu(
            menuItem("Variants overview", tabName = "overview", icon = icon("th")),
            menuItem("Filters settings", tabName = "filters", icon = icon("th")),
            menuItem("Filter explorer", tabName = "filter_explorer", icon = icon("th")),
            menuItem("Filter results", tabName = "filter_results", icon = icon("th")),
            menuItem("PanelApp and gene lists", tabName = "gene_lists", icon = icon("th")),
            menuItem("Gene details", tabName = "gene_details", icon = icon("th"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "overview",
                    plotOutput("ped"),
                    plotOutput("Var_consequence_plot"),
                    plotOutput("Var_type_plot"),
                    plotOutput("Var_PopAF_plot"),
                    plotOutput("Gene_segregation_plot")    
            ),
            tabItem(tabName = "filters",
                    column(12, align="center", actionButton(inputId = "Apply_filters", label = "Apply filters")),
                    fluidRow(hr()),
                    
                    fluidRow(
                        box(title = "Variant filters", width = 12, status = "primary", solidHeader = TRUE,
                        #tabPanel("General",
                        fluidRow(
                            column(4,
                                h4("General filters"),
                                uiOutput("d_score"),
                                uiOutput("Max_pop_AF"),
                                uiOutput("Cohort_AF"),
                                uiOutput("var_consequence") ),
                            column(4,
                                h4("Missense filters"),
                                uiOutput("CADD"),
                                uiOutput("REVEL") ),
                            column(4,
                                h4("Splice filters"),   
                                uiOutput("spliceAI")) ),
                        fluidRow(
                            column(4, 
                                h4("Regulatory filters"),
                                uiOutput("DANN"),
                                uiOutput("ReMM"),
                                uiOutput("PhyloP100") ),
                            column(4,
                                h4("Compound het filter"),
                                "At least one of the variant must have the selected consequence",
                                uiOutput("comphet_consequence") ) )
                        )
                    ),
                        #tabPanel("Splicing", uiOutput("spliceAI"))
                    #) ),
                    
                    fluidRow(
                        box(title = "Gene filters", id = "gene_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        uiOutput("pLI") ),
            
                        box(title = "Segregation filters", id = "segregation_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        fluidRow(
                            column(6, textOutput("Total_affected")),
                            column(6, textOutput("Total_unaffected")) ),
                        uiOutput("segregation_controls"),
                                #column(3, uiOutput("recessive")),
                                #column(3, uiOutput("dominant")),
                                #column(3, uiOutput("denovo")),
                                #column(3, uiOutput("comphet")) ),
                        "Number of affected / unaffected carrying the variant with the given genotype",
                        "Homozygous and heterozygous conditions evaluated as AND, comphet as OR"
                        ),
                    )
            ),
            tabItem(tabName = "filter_explorer",
                    h3("Summary of filters effect"),
                    fluidRow(column(6,plotlyOutput("summary_variants_filters", width = "100%")), column(6,plotlyOutput("summary_genes_filters", width = "100%"))),
                    
                    h3("Genes filter evaluation"),
                    fluidRow(column(6,selectInput("genes_X_axis", h4("X axis:"), choices = genes_axes_options, selected = "Gado_zscore")),
                             column(6,selectInput("genes_Y_axis", h4("Y axis:"), choices = genes_axes_options, selected = "pLI_gnomAD"))
                    ),
                    plotlyOutput("genes_scatter",width = "100%"),
                    br(),
                    
                    h3("Variants filter evaluation"),
                    plotSelectedUI("variants_scatter", variables=list("x"=variants_axes_options,"y"=variants_axes_options), plotly=TRUE),
                    br(),
                    plotSelectedUI("variants_barplot", variables=list("x"=variants_bar_options), plotly=TRUE)    
            ),
            tabItem(tabName = "filter_results",
                    fluidRow(
                        column(8, plotlyOutput("GADO_rank")),
                        column(4, plotlyOutput("GADO_distribution"))
                    ),
                    DT::dataTableOutput("genesTable"),
                    h3("Custom list genes"),
                    DT::dataTableOutput("customGenesTable"),
                    hr(),
                    fluidRow(column(8), 
                             column(3, offset=1,downloadObjUI(id = "save_results")))    
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
                    br(),
                    box(title = "Pathways and ontology", id = "path_go_details", status = "info", solidHeader = TRUE, width = 12,
                        collapsible = TRUE, collapsed = TRUE,
                        div(DT::dataTableOutput("geneInfo"),style="font-size: 75%") ),
                    box(title = "GTeX median expression", id = "gtex_expression", status = "info", solidHeader = TRUE, width=12,
                        collapsible = TRUE, collapsed = TRUE,
                        plotOutput("gtex_plot") )
                    
            )
        )
    )
)

########################
### SERVER FUNCTIONS ###
########################

server <- function(input, output) {
    
    ################
    ### Load Files 
    ################
    
    observeEvent(input$CaseCode, {
        RV$load_status = ""
        RV$custom_genes = character()
        #reset("custom_file")
        
        gado_file <- paste0(data_dir, "/", paste0(input$CaseCode,".txt"))
        variants_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.vars.tsv.gz")))
        genes_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.genes.tsv.gz")))
        comphet_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.comphet.tsv.gz")))

        gado_score <- read.table(gado_file, sep="\t", header=T, stringsAsFactors = F)
        RV$gado90 <- quantile(gado_score$Zscore, .9)[[1]]
        
        RV$comphet_df <- read.table(comphet_file, header=T, sep="\t", stringsAsFactors=F)
        RV$comphet_df$ID <- paste0("CompHet_", RV$comphet_df$ID)
        RV$comphet_df$Class <- "PASS"
        
        RV$variants_df <- read.table(variants_file, header=T, sep="\t", stringsAsFactors = F)
        
        RV$variants_df$MaxPopAF[is.na(RV$variants_df$MaxPopAF)] <- 0
        RV$variants_df$Consequence[is.na(RV$variants_df$Consequence)] <- "non_coding_region"
        RV$variants_df[is.na(RV$variants_df)] <- 0
        RV$variants_df$Class <- "PASS"

        seg_df1 <- as.data.frame(
            RV$comphet_df %>% select(ID,num_aff) %>% 
                mutate(hom_aff=0, hom_unaff=0, het_aff=0, het_unaff=0) %>%
                rename(comphet_aff = num_aff) )
        seg_df2 <- as.data.frame(
            RV$variants_df %>% select(ID,hom_aff,hom_unaff,het_aff,het_unaff) %>%
                mutate(comphet_aff = 0) )
        RV$segregation_df <- rbind(seg_df1, seg_df2)
        
        RV$genes_df <- read.table(genes_file, header=T, sep="\t", stringsAsFactors = F)
        RV$genes_df$Class <- "PASS"
        RV$genes_scores <- as.data.frame(RV$genes_df %>% select(Gene,Gado_zscore,Exomiser_GenePhenoScore,pLI_exac,pLI_gnomad,Class) %>% distinct() %>% arrange(desc(Gado_zscore)))
        RV$genes_df <- as.data.frame(RV$genes_df %>% select(-ID,-Gado_zscore,-Exomiser_GenePhenoScore,-pLI_exac,-pLI_gnomad))
        
        genes_with_comphet <- as.data.frame(
            RV$comphet_df %>% select(Gene,ID) %>% group_by(Gene) %>% 
                mutate(Variants=paste(ID,collapse = ","), Variants_n=n(), Inh_model="comphet", Class="PASS") %>% 
                select(-ID) %>% 
                distinct()
        )
        RV$genes_df <- rbind(RV$genes_df, genes_with_comphet)
        RV$genes_df <- as.data.frame(RV$genes_df %>% separate_rows(Variants, sep=","))
        RV$genes_scores$pLI_gnomad[is.na(RV$genes_scores$pLI_gnomad)] <- 0
        RV$genes_scores$pLI_exac[is.na(RV$genes_scores$pLI_exac)] <- 0

        ###TEMP - remove vars not following any seg model
        #RV$genes_df <- as.data.frame(RV$genes_df %>% filter(Inh_model != "other"))
        #RV$variants_df <- as.data.frame(RV$variants_df %>% filter(ID %in% RV$genes_df$Variants))
        ###
        
        RV$total_affected <- sum(all_peds[input$CaseCode]$affected)
        RV$total_unaffected <- length(all_peds[input$CaseCode]$id) - RV$total_affected 
        RV$values_affected <- seq(0,RV$total_affected)
        names(RV$values_affected) <- seq(0,RV$total_affected)
        RV$values_affected <- c("NOT_ALLOWED" = (RV$total_affected + 1), RV$values_affected)
        RV$values_unaffected <- seq(0,RV$total_unaffected)
    })
    
    observeEvent(input$custom_file, {
        req(input$custom_file)
        RV$custom_genes <- scan(input$custom_file$datapath,what="",sep="\n")
        RV$load_status <- paste0(length(RV$custom_genes), " genes loaded")
    })
    
    output$Loading_result <- renderText({
        RV$load_status
    })
    
    output$Disease <- renderText({
        paste0(HPOs$Disease[HPOs$X.CaseID == input$CaseCode])
    })
    
    #######################################
    ### Apply filters on button click
    #######################################
        
    observeEvent(input$Apply_filters, {
        if ("ALL" %in% input$Consequence) { 
            RV$accepted_consequence <- sort(unique(RV$variants_df$Consequence))    
        } else {
            RV$accepted_consequence <- input$Consequence
        }
        
        if ("ALL" %in% input$comphet_consequence_filter) { 
            RV$comphet_consequence <- sort(unique(RV$variants_df$Consequence))    
        } else {
            RV$comphet_consequence <- input$comphet_consequence_filter
        }
    })    

    #######################################
    ### Data reactive tables configuration
    #######################################
    
    segregating_vars <- reactive({
        input$Apply_filters
        
        isolate(
            callModule(segregationModule, "segregation", segregation_df = RV$segregation_df, cols_names = segregation_cols)
        )
    })
    
    variants_df <- reactive({
        input$Apply_filters
        
        isolate({
            #first select variants that pass the filters
            pass_vars <- as.data.frame(RV$variants_df %>% 
                filter(! (
                d_score < input$d_score_filter | 
                    Consequence %nin% RV$accepted_consequence |
                    MaxPopAF > input$MaxPopAF_filter |
                    cohortAF > input$CohortAF_filter |
                    (Reg_type == "splicing" & 
                         (SpliceAI_SNP_SpliceAI_max < input$spliceAI_filter & 
                              SpliceAI_DEL_SpliceAI_max < input$spliceAI_filter) ) |
                    (Consequence == "missense_variant" &
                         (CADD_PhredScore < input$CADD_filter & 
                              REVEL_score < input$REVEL_filter)) |
                    (Consequence %in% reg_vars &
                         (ReMM_score < input$ReMM_filter |
                              DANN_DANN < input$DANN_filter |
                              PhyloP100 < input$PhyloP100_filter))
                )))
            
            #COMPHET FILTERING
            #Select comphet where both vars pass variants filters
            comphet_accepted <- RV$comphet_df %>% filter(
                    V1 %in% pass_vars$ID &
                    V2 %in% pass_vars$ID)
            
            #COMPHET specific consequence
            #This filter allow to select comphet where at least 1 var as the given consequence
            #From pass vars select only the ones with accepted consequence
            #based on the comphet consequence filter
            comphet_accepted_vars <- as.data.frame(pass_vars %>%
                filter(Consequence %in% RV$comphet_consequence) )
            
            #Get the final list of vars in accepted comphet
            #from comphet select combo when: 
            #both vars in the combo passed the general variants filters (comphet_accepted df)
            #at least one var in the combo pass the comphet consequence filter
            #the combo passed segregation filters
            comphet_single_vars <- comphet_accepted %>% filter(
                ID %in% segregating_vars() &
                (V1 %in% comphet_accepted_vars$ID |
                V2 %in% comphet_accepted_vars$ID) ) %>% gather(key="Variant",value="VarID",V1:V2) %>% select(VarID)
            
            #GET THE FINAL LIST OF ACCEPTED VARS
            #the final list of segregating accepted vars is now equal to
            #single vars passing filters and segregation
            #vars part of a comphet passing all filters + segregation
            accepted_vars_list <- intersect(segregating_vars(), pass_vars$ID)
            accepted_vars_list <- unique(c(accepted_vars_list, comphet_single_vars$VarID))
            
            #return variants dataframe with updated Class column (PASS/FILTER)
            as.data.frame(RV$variants_df %>% mutate(Class = ifelse(
                ID %in% accepted_vars_list,
                "PASS","FILTER") ) )
        })
    })
    
    comphet_df <- reactive({
        filtered_vars_list <- variants_df()$ID[variants_df()$Class == "PASS"]
        
        as.data.frame(RV$comphet_df %>% mutate(Class = ifelse(
            V1 %in% filtered_vars_list & 
            V2 %in% filtered_vars_list &
            ID %in% segregating_vars(), "PASS", "FILTER")))
    })
    
    genes_df <- reactive({
        RV$filtered_vars_list <- unique(c(
            variants_df()$ID[variants_df()$Class == "PASS"], 
            comphet_df()$ID[comphet_df()$Class == "PASS"])) 
        
        genes_above_score <- as.data.frame(RV$genes_scores %>% mutate(Class = ifelse(
            pLI_gnomad >= input$pLI_filter, "PASS", "FILTER")))
            
        as.data.frame(RV$genes_df %>% mutate(Class = ifelse(
            Variants %in% RV$filtered_vars_list & 
            Gene %in% genes_above_score$Gene, "PASS", "FILTER")))    
    })
    
    genes_scores <- reactive({
        RV$filtered_genes_list <- unique(genes_df()$Gene[genes_df()$Class == "PASS"]) 
        as.data.frame(RV$genes_scores %>% mutate(Class = ifelse(
            Gene %in% RV$filtered_genes_list, "PASS", "FILTER")))
    })
    
    filters_summ_genes <- reactive ({
        PASS_counts <- c(
        #genes_df() %>% filter(GDI_phred <= input$GDI_filter) %>% nrow(),
        genes_scores() %>% filter(pLI_gnomad >= input$pLI_filter) %>% nrow() )
        
        filters_summ_genes <- data.frame(Filter=c("pLI gnomAD"), 
            PASS=PASS_counts, 
            FILTERED=(nrow(RV$genes_df)-PASS_counts) )
        filters_summ_genes <- gather(filters_summ_genes, key="Class", value="Count", PASS:FILTERED)   
    })
    
    filters_summ_vars <- reactive ({
        tot_vars <- RV$variants_df %>% select(Chr,Pos,Ref,Alt) %>% nrow()
        PASS_counts <- c(
            variants_df() %>% filter(d_score >= input$d_score_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter((Reg_type == "splicing" & SpliceAI_SNP_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing") %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(MaxPopAF <= input$MaxPopAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(Consequence %in% RV$accepted_consequence) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(cohortAF <= input$CohortAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow()
        )
        filters_summ_vars <- data.frame(Filter=c("d score","SpliceAI","MaxPop AF","Consequence","cohort AF"),
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
        paste0("\tTotal number of affected individuals in this pedigree: ", RV$total_affected)
    })
    
    output$Total_unaffected <- renderText({
        paste0("\tTotal number of unaffected individuals in this pedigree: ", RV$total_unaffected)
    })
    
    output$ped <- renderPlot({
        plot.pedigree(all_peds[input$CaseCode], mar = c(5, 3, 5, 3))
    })
    
    output$Var_consequence_plot <- renderPlot({
        ggplot(RV$variants_df, aes(x=Consequence)) + geom_bar()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    output$Var_type_plot <- renderPlot({
        ggplot(RV$variants_df, aes(x=VarType)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    output$Var_PopAF_plot <- renderPlot({
        ggplot(RV$variants_df, aes(x=MaxPopAF, y=cohortAF, color=VarType)) + geom_point()
    })
    
    output$Gene_segregation_plot <- renderPlot({
        gene_segregation_df <- as.data.frame(RV$genes_df %>% select(Gene,Inh_model) %>% distinct())
        ggplot(gene_segregation_df, aes(x=Inh_model)) + geom_bar() + labs(y="genes count")
    })
    
    #################
    ### Filters tab
    #################
    
    output$d_score <- renderUI({
        startvalue = 0
        sliderInput("d_score_filter", "Min d_score:",
                    min = 0,
                    max = max(RV$variants_df$d_score, na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    output$spliceAI <- renderUI({
        startvalue = 0
        sliderInput("spliceAI_filter", "Min spliceAI score:",
                    min = 0,
                    max = max(RV$variants_df$SpliceAI_SNP_SpliceAI_max, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$CADD <- renderUI({
        startvalue = 0
        sliderInput("CADD_filter", "Min CADD phred:",
                    min = 0,
                    max = max(RV$variants_df$CADD_PhredScore, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$DANN <- renderUI({
        startvalue = 0
        sliderInput("DANN_filter", "Min DANN score:",
                    min = 0,
                    max = max(RV$variants_df$DANN_DANN, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$ReMM <- renderUI({
        startvalue = 0
        sliderInput("ReMM_filter", "Min ReMM score:",
                    min = 0,
                    max = max(RV$variants_df$ReMM_score, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$PhyloP100 <- renderUI({
        startvalue = 0
        sliderInput("PhyloP100_filter", "Min PhyloP100 score:",
                    min = 0,
                    max = max(RV$variants_df$PhyloP100, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$REVEL <- renderUI({
        startvalue = 0
        sliderInput("REVEL_filter", "Min REVEL score:",
                    min = 0,
                    max = max(RV$variants_df$REVEL_score, na.rm = T),
                    value = 0,
                    step = 0.02)
    })
    
    output$Max_pop_AF <- renderUI({
        startvalue = 0
        sliderInput("MaxPopAF_filter", "Max gnomAD AF:",
                    min = 0,
                    max = max(RV$variants_df$MaxPopAF, na.rm = T),
                    value = max(RV$variants_df$MaxPopAF, na.rm = T),
                    step = 0.001)
    })
    
    output$Cohort_AF <- renderUI({
        startvalue = 0
        sliderInput("CohortAF_filter", "Max cohort AF:",
                    min = 0,
                    max = max(RV$variants_df$cohortAF, na.rm = T),
                    value = max(RV$variants_df$cohortAF, na.rm = T),
                    step = 0.001)
    })
    
    output$pLI <- renderUI({
        startvalue = 0
        sliderInput("pLI_filter", "Min gnomad pLI:",
                    min = 0,
                    max = max(RV$genes_scores$pLI_gnomad, na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    #output$GDI <- renderUI({
    #    startvalue = 0
    #    sliderInput("GDI_filter", "Max phred GDI:",
    #                min = 0,
    #                max = max(RV$genes_df$GDI_phred, na.rm = T),
    #                value = max(RV$genes_df$GDI_phred, na.rm = T),
    #                step = 0.1)
    #})
    
    output$segregation_controls <- renderUI({
        segregationUI("segregation", choices_affected = RV$values_affected, choices_unaffected=RV$values_unaffected)
    })
    
    output$var_consequence <- renderUI({
        var_types <- sort(unique(RV$variants_df$Consequence))
        names(var_types) <- sort(unique(RV$variants_df$Consequence))
        var_types <- c("ALL" = "ALL", var_types)
        selectInput("Consequence", "Variant consequence:", choices = var_types, multiple=TRUE, selected="ALL")
    })
    
    output$comphet_consequence <- renderUI({
        var_types <- sort(unique(RV$variants_df$Consequence))
        names(var_types) <- sort(unique(RV$variants_df$Consequence))
        var_types <- c("ALL" = "ALL", var_types)
        selectInput("comphet_consequence_filter", "Variant consequence:", choices = var_types, multiple=TRUE, selected="ALL")
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
    
    output$genes_scatter <- renderPlotly({
        ggplotly(
            ggplot(genes_scores(), aes_string(x=input$genes_X_axis, y=input$genes_Y_axis, color="Class", label="Gene")) + geom_point(size=1) 
        )
    })
    
    callModule(plotSelectedModule, "variants_scatter", plot_data = variants_df(), plotType = "scatter", variables = list("color"="Class"))
    callModule(plotSelectedModule, "variants_barplot", plot_data = variants_df(), plotType = "barplot", variables = list("fill"="Class"), additionalOptions = list(format1, scale_y_sqrt()))
    
    ########################
    ### Filter Results tab
    ########################
    
    output$GADO_rank <- renderPlotly({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        if (length(input$genesTable_rows_selected) > 0) {
            Zscore <- genedetail_tab[input$genesTable_rows_selected, "Gado_zscore"]
            GADO_plot <- ggplot(genes_scores(), aes(x=Gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = Zscore, color="red") + geom_vline(xintercept = RV$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        } else {
            GADO_plot <- ggplot(genes_scores(), aes(x=Gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = RV$gado90, linetype = "dashed") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        }
        ggplotly( GADO_plot )
    })
    
    output$GADO_distribution <- renderPlotly({
        validate(need(input$genesTable_rows_selected, "Select a gene"))
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        Zscore <- genedetail_tab[input$genesTable_rows_selected, "Gado_zscore"]
        GADO_dist <- ggplot(gado_distribution[gado_distribution$Hgnc == gene_name,], aes(x=Zscore)) + 
            geom_histogram(bins = max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name])*10) + 
            geom_vline(xintercept = Zscore, color = "red") +
            scale_x_continuous(breaks=round(seq(min(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name]),max(gado_distribution$Zscore[gado_distribution$Hgnc == gene_name]),1)))
        ggplotly( GADO_dist )
     })
    
    output$genesTable <- DT::renderDataTable(selection="single", {
        as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class)) 
    })
    
    genesTable_proxy <- DT::dataTableProxy("genesTable")
    
    customGenesTable_df <- reactive({
        as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% mutate(row_idx = row_number()) %>% filter(Gene %in% RV$custom_genes) %>% select(-Class))
    })
    
    output$customGenesTable <- DT::renderDataTable(selection="single", {
        customGenesTable_df()
    })
    
    observeEvent(input$customGenesTable_rows_selected, {
        mydf <- customGenesTable_df()
        row_idx <- mydf[input$customGenesTable_rows_selected,"row_idx"]
        RV$load_stats <- row_idx
        genesTable_proxy %>% selectRows(row_idx) 
    })
    
    callModule(downloadObj, id = "save_results",
        output_prefix=input$CaseCode,
        output_data = list(
            "genes.tsv" = as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class)),
            "customGenes.tsv" = as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% filter(Gene %in% RV$custom_genes) %>% select(-Class)),
            "variants.tsv" = as.data.frame(variants_df() %>% filter(Class == "PASS", Gene %in% RV$filtered_genes_list)),
            "comphet.tsv" = as.data.frame(comphet_df() %>% filter(Class == "PASS", Gene %in% RV$filtered_genes_list) %>% 
                                gather(key="Variant",value = "varID", V1:V2) %>% 
                                inner_join(., variants_df()[variants_df()$Class == "PASS",], by=c("varID"="ID")) %>%
                                select(-Class.x,-Class.y))),
        zip_archive = paste0(input$CaseCode, ".results.zip")
    )

    ################################
    ### PanelApp and genes lists tab
    ################################
    
    output$PanelApp_panels_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
        PanelApp_panels_df()
    })
    
    output$PanelApp_genes_table <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        validate(need(input$PanelApp_panels_table_rows_selected, "Select one panel"))
        
        panelID <- PanelApp_panels_df()[input$PanelApp_panels_table_rows_selected, "id"]

        PanelApp_genes_df <- as.data.frame(PanelApp_genes %>% filter(panel_idx == panelID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., genes_scores(), by = c("entity_name" = "Gene")))
        PanelApp_genes_df <- PanelApp_genes_df[order(PanelApp_genes_df$Gado_zscore, decreasing = TRUE),]        
    })
    
    output$ClinVar_table <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        ClinVar_df()
    })
    
    output$geneLists_table <- DT::renderDataTable(selection="single", options = list(scrollX = TRUE), {
        geneLists_df()
    })
    
    output$geneLists_genes_table <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        validate(need(input$PanelApp_panels_table_rows_selected, "Select one gene list"))
        
        listID <- geneLists_df()[input$geneLists_table_rows_selected, "id"]
        
        geneList_genes_df <- as.data.frame(geneLists_genes %>% filter(genelist_idx == listID, entity_name %in% RV$filtered_genes_list) %>%  left_join(., genes_scores(), by = c("entity_name" = "Gene")))
        geneList_genes_df <- geneList_genes_df[order(geneList_genes_df$Gado_zscore, decreasing = TRUE),]        
    })
    
    ####################
    ### Gene detail tab
    ####################
    
    output$Gene_symbol <- renderText({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        gene_name
    })
    
    output$Gene_name <- renderText({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        genes_info[genes_info$symbol == gene_name, "name"]
    })
    
    output$GTeX_link <- renderUI({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        ensg_id <- genes_info[genes_info$symbol == gene_name, "ensembl_gene_id"]
        tags$a(href=paste0("https://gtexportal.org/home/gene/", gene_name), paste0(gene_name, "(", ensg_id, ")"), target="_blank")
    })
    
    output$geneDetail <- DT::renderDataTable(selection="none", {
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        genedetail_tab[input$genesTable_rows_selected,]
    })
    
    output$geneInfo <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        selected <- list()
        
        for (n in names(gene_anno)) {
            id_list <- grep(gene_name, gene_anno[[n]])
            selected[[n]] <- names(gene_anno[[n]])[id_list]
        }
        Additional_info <- data.frame(lapply(selected, "length<-", max(lengths(selected))), stringsAsFactors = F)
    })
    
    output$gtex_plot <- renderPlot ({
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        gene_exp <- GTeX_data[GTeX_data$Description == gene_name,]
        validate(need(nrow(gene_exp)>0, "Gene not found in GTeX"))
        gene_exp <- gather(gene_exp, key="tissue", value="median_TPM", 3:ncol(gene_exp))
        ggplot(gene_exp, aes(x=tissue, y=median_TPM)) + geom_bar(stat="identity") + format1
    })
    
    output$variantsTable <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS"))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        
        comphet_details <- comphet_df() %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", V1:V2)
        
        #as.data.frame(variants_df() %>% filter(Class == "PASS", Gene == gene_name))
        as.data.frame(variants_df() %>% filter(Class == "PASS", Gene == gene_name, ID %nin% comphet_details$varID))
    })
    
    output$comphetTable <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS"))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        
        comphet_details <- comphet_df() %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", V1:V2) %>% select(-Gene, -Variant) %>% distinct()
        validate(need(nrow(comphet_details)>0, 'No compound het variants'))
        comphet_details <- merge(comphet_details,variants_df()[variants_df()$Class == "PASS",], by.x="varID",by.y="ID")
        as.data.frame(comphet_details %>% select(-Class.x,-Class.y,) %>% arrange(ID))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

# VARIANT EXPLORER
# Author: Edoardo Giacopuzzi
# Explore and filter annotated variants from V2

# Input are tables from the VARAN V2 filtering and annotation

# TODO Add ClinVar pathogenic gene filter
# TODO Add separate table reporting all Clinvar pathogenic vars

library(shiny)
library(DT)
library(plotly)
library(ggplot2)
library(kinship2)
library(tidyr)
library(shinydashboard)
source("plotModule.R")
source("downloadModule.R")

#Set data dir containing variants tables and look for files
data_dir <- "example_data"

files <- list.files(data_dir, pattern = "*.vars.tsv.gz")
samplesID <- gsub("V2\\.|\\.var2reg.vars.tsv.gz","",files)

#We can switch to read from a input table containing sampleID and ped locations...
#samples_df <- data.frame(ID=samplesID,file=files)

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
genes_info_file <- "hgnc_complete_set.txt"
genes_info <- read.table(genes_info_file, header=T, sep="\t", stringsAsFactors = F)

##Load pathways and GO groups (suppose .gmt file from MSigDB)
gene_anno <- list()
gene_anno[["pathways"]] <- readGMT("c2.cp.v7.0.symbols.gmt")
gene_anno[["GO_BP"]] <- readGMT("c5.bp.v7.0.symbols.gmt")
gene_anno[["GO_MF"]] <- readGMT("c5.mf.v7.0.symbols.gmt")
gene_anno[["GO_CC"]] <- readGMT("c5.cc.v7.0.symbols.gmt")

## Load GTeX median expression
GTeX_file <- gzfile("GTEx_v8_median_TPM.gct.gz")
GTeX_data <- read.table(GTeX_file,sep="\t",header=T, stringsAsFactors = F)

#######################
### Set environment ###
#######################

`%nin%` = Negate(`%in%`)

##Set axes options for plots
genes_axes_options <- c(
    "GADO Zscore"= "Gado_zscore", 
    "Exomiser Pheno Rank"= "Exomiser_GenePhenoScore",
    "gnomAD pLI"= "pLI_gnomAD", 
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

##Set reactive objects
RV <- reactiveValues(
    variant_df = data.frame(), 
    genes_scores = data.frame(), 
    genes_df = data.frame(),
    comphet_df = data.frame(),
    load_status = character(),
    custom_genes = character(),
    filters_summ_genes = data.frame(),
    filters_summ_vars = data.frame())

######################
### USER INTERFACE ###
######################

ui <- dashboardPage(
    dashboardHeader(
        title = "Variant Explorer",
        dropdownMenu(type = "notifications",
                     notificationItem(
                         text = "Beta version now working only on 004Int001",
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
        textInput(inputId = "custom_file",label = "File:",value = "", placeholder = "gene list file"),
        actionButton(inputId = "Load_file", label = "Load file"),
        textOutput("Loading_result"),
        sidebarMenu(
            menuItem("Variants overview", tabName = "overview", icon = icon("th")),
            menuItem("Filters settings", tabName = "filters", icon = icon("th")),
            menuItem("Filter explorer", tabName = "filter_explorer", icon = icon("th")),
            menuItem("Filter results", tabName = "filter_results", icon = icon("th")),
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
                            column(6,
                                h4("General filters"),
                                uiOutput("d_score"),
                                uiOutput("Max_pop_AF"),
                                uiOutput("Cohort_AF"),
                                uiOutput("var_consequence") ),
                            column(6,
                                h4("Splice filters"),   
                                uiOutput("spliceAI")) 
                        ) 
                    ),
                        #tabPanel("Splicing", uiOutput("spliceAI"))
                    #) ),
                    
                    fluidRow(
                        box(title = "Gene filters", id = "gene_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        uiOutput("pLI"),
                        uiOutput("GDI") ),
            
                        box(title = "Segregation filters", id = "segregation_filters", status = "primary", solidHeader = TRUE,
                        collapsible = TRUE,
                        textOutput("Total_affected"),
                        fluidRow(
                                column(3, uiOutput("recessive")),
                                column(3, uiOutput("dominant")),
                                column(3, uiOutput("denovo")),
                                column(3, uiOutput("comphet")) ),
                        "Number of affected individuals in which a variants segregates",
                        "Works with OR logic"
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
                    plotlyOutput("GADO_rank"),
                    DT::dataTableOutput("genesTable"),
                    h3("Custom list genes"),
                    DT::dataTableOutput("customGenesTable"),
                    hr(),
                    fluidRow(column(8), 
                             column(3, offset=1,downloadObjUI(id = "save_results")))    
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
        variants_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.vars.tsv.gz")))
        genes_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.genes.tsv.gz")))
        comphet_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.comphet.tsv.gz")))
        segregation_file <- gzfile(paste0(data_dir, "/", paste0("V2.",input$CaseCode,".var2reg.segregation.tsv.gz")))
        
        RV$comphet_df <- read.table(comphet_file, header=T, sep="\t", stringsAsFactors=F)
        RV$comphet_df$Class <- "PASS"
        
        RV$segregation_df <- read.table(segregation_file, header=T, sep="\t", stringsAsFactors=F)
        
        RV$variants_df <- read.table(variants_file, header=T, sep="\t", stringsAsFactors = F)
        RV$variants_df$MaxPopAF[is.na(RV$variants_df$MaxPopAF)] <- 0
        RV$variants_df$cohortAF[is.na(RV$variants_df$cohortAF)] <- 0
        RV$variants_df$SpliceAI_SpliceAI_max[is.na(RV$variants_df$SpliceAI_SpliceAI_max)] <- 0
        RV$variants_df$Consequence[is.na(RV$variants_df$Consequence)] <- "non_coding_region"
        RV$variants_df$Class <- "PASS"

        RV$genes_df <- read.table(genes_file, header=T, sep="\t", stringsAsFactors = F)
        RV$genes_df <- as.data.frame(RV$genes_df %>% separate_rows(Variants, sep=","))
        RV$genes_df$GDI_phred[is.na(RV$genes_df$GDI_phred)] <- max(RV$genes_df$GDI_phred, na.rm = T)
        RV$genes_df$pLI_gnomAD[is.na(RV$genes_df$pLI_gnomAD)] <- 0
        RV$genes_df$Class <- "PASS"
        
        ###TEMP - remove vars not following any seg model
        #RV$genes_df <- as.data.frame(RV$genes_df %>% filter(Inh_model != "other"))
        #RV$variants_df <- as.data.frame(RV$variants_df %>% filter(ID %in% RV$genes_df$Variants))
        ###

        RV$genes_scores <- as.data.frame(RV$genes_df %>% select(Gene,Gado_zscore,Exomiser_GenePhenoScore,GDI_phred,pLI_gnomAD,Class) %>% distinct() %>% arrange(desc(Gado_zscore)))
        
        RV$total_affected <- sum(all_peds[input$CaseCode]$affected)
        RV$values_segregation <- c(RV$total_affected+1,seq(1:RV$total_affected))
        names(RV$values_segregation) <- c("NOT_ACCEPTED", seq(1:RV$total_affected))
    })
    
    observeEvent(input$Load_file, {
        if (input$custom_file == "") {
            RV$load_status = "Insert a file name"
        } else if (!file.exists(input$custom_file)) {
            RV$load_status = "File do not exists"
        } else {
            RV$custom_genes <- scan(input$custom_file,what="",sep="\n")
            RV$load_status <- paste0(length(RV$custom_genes), " genes loaded")
        }
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
    })    
    #    RV$segregating_vars_df <- as.data.frame(RV$segregation_df %>% 
    #        filter(
    #            recessive >= input$recessive_filter | 
    #            dominant >= input$dominant_filter |
    #            deNovo >= input$denovo_filter |
    #            comphet >= input$comphet_filter) )
    #    cat(file=stderr(), "segregating df: ", length(RV$segregating_vars_df$ID),"\n")
        
    #    segregating_vars_list <- RV$segregating_vars_df$ID
    #    comphet_single_vars <- RV$comphet_df %>% filter(ID %in% segregating_vars_list) %>% gather(key="Variant",value="VarID",var1:var2) %>% select(VarID)
    #    segregating_vars_list <- unique(c(segregating_vars_list, comphet_single_vars$VarID))
        
    #    RV$variants_df <- as.data.frame(RV$variants_df %>% mutate(Class = ifelse(
    #        ID %in% segregating_vars_list &
    #            d_score >= input$d_score_filter & 
    #            Consequence %in% RV$accepted_consequence &
    #            MaxPopAF <= input$MaxPopAF_filter &
    #            cohortAF <= input$CohortAF_filter &
    #            ((Reg_type == "splicing" & SpliceAI_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing"),"PASS","FILTER")))
    #    RV$filtered_vars_list <- RV$variants_df$ID[RV$variants_df$Class == "PASS"]
        
    #    RV$comphet_df <- as.data.frame(RV$comphet_df %>% mutate(Class = ifelse(
    #        var1 %in% RV$filtered_vars_list & 
    #            var2 %in% RV$filtered_vars_list &
    #            ID %in% RV$segregating_vars_df$ID, "PASS", "FILTER")))
    #    
    #    RV$filtered_vars_list <- unique(c(
    #        RV$variants_df$ID[RV$variants_df$Class == "PASS"], 
    #        RV$comphet_df$ID[RV$comphet_df$Class == "PASS"])) 

    #    RV$genes_df <- as.data.frame(RV$genes_df %>% mutate(Class = ifelse(
    #        Variants %in% RV$filtered_vars_list & 
    #            GDI_phred <= input$GDI_filter &
    #            pLI_gnomAD >= input$pLI_filter, "PASS", "FILTER")))
    #    
    #    RV$filtered_genes_list <- unique(RV$genes_df$Gene[RV$genes_df$Class == "PASS"])
        
    #    RV$genes_scores <- as.data.frame(RV$genes_scores %>% mutate(Class = ifelse(
    #        Gene %in% RV$filtered_genes_list, "PASS", "FILTER")))
    
    #    PASS_counts <- c(
    #        RV$genes_df %>% filter(GDI_phred <= input$GDI_filter) %>% nrow(),
    #        RV$genes_df %>% filter(pLI_gnomAD >= input$pLI_filter) %>% nrow() )
        
    #    RV$filters_summ_genes <- data.frame(Filter=c("pLI gnomAD","GDI phred"), 
    #                                     PASS=PASS_counts, 
    #                                     FILTERED=(nrow(RV$genes_df)-PASS_counts) )
    #    RV$filters_summ_genes <- gather(RV$filters_summ_genes, key="Class", value="Count", PASS:FILTERED)
        
    #    tot_vars <- RV$variants_df %>% select(Chr,Pos,Ref,Alt) %>% nrow()
    #    PASS_counts <- c(
    #        RV$variants_df %>% filter(d_score >= input$d_score_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
    #        RV$variants_df %>% filter((Reg_type == "splicing" & SpliceAI_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing") %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
    #        RV$variants_df %>% filter(MaxPopAF <= input$MaxPopAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
    #        RV$variants_df %>% filter(Consequence %in% RV$accepted_consequence) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
    #        RV$variants_df %>% filter(cohortAF <= input$CohortAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow()
    #    )
    #    RV$filters_summ_vars <- data.frame(Filter=c("d score","SpliceAI","MaxPop AF","Consequence","cohort AF"),
    #                                    PASS=PASS_counts, 
    #                                    FILTERED=(tot_vars-PASS_counts))
    #    RV$filters_summ_vars <- gather(RV$filters_summ_vars, key="Class", value="Count", PASS:FILTERED)
    #})
    
    #######################################
    ### Data reactive tables configuration
    #######################################
    
    segregating_vars_df <- reactive({
        input$Apply_filters
        
        isolate(
            as.data.frame(RV$segregation_df %>% filter(
            recessive >= input$recessive_filter | 
            dominant >= input$dominant_filter |
            deNovo >= input$denovo_filter |
            comphet >= input$comphet_filter) )
        )
    })
    
    variants_df <- reactive({
        input$Apply_filters
        
        isolate({
            segregating_vars_list <- segregating_vars_df()$ID
            comphet_single_vars <- RV$comphet_df %>% filter(ID %in% segregating_vars_list) %>% gather(key="Variant",value="VarID",var1:var2) %>% select(VarID)
            segregating_vars_list <- unique(c(segregating_vars_list, comphet_single_vars$VarID))

            as.data.frame(RV$variants_df %>% mutate(Class = ifelse(
                ID %in% segregating_vars_list &
                d_score >= input$d_score_filter & 
                Consequence %in% RV$accepted_consequence &
                MaxPopAF <= input$MaxPopAF_filter &
                cohortAF <= input$CohortAF_filter &
                ((Reg_type == "splicing" & SpliceAI_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing"),"PASS","FILTER")))
        })
    })
    
    comphet_df <- reactive({
        filtered_vars_list <- variants_df()$ID[variants_df()$Class == "PASS"]
        
        as.data.frame(RV$comphet_df %>% mutate(Class = ifelse(
            var1 %in% filtered_vars_list & 
            var2 %in% filtered_vars_list &
            ID %in% segregating_vars_df()$ID, "PASS", "FILTER")))
    })
    
    genes_df <- reactive({
        RV$filtered_vars_list <- unique(c(
            variants_df()$ID[variants_df()$Class == "PASS"], 
            comphet_df()$ID[comphet_df()$Class == "PASS"])) 
            
        as.data.frame(RV$genes_df %>% mutate(Class = ifelse(
            Variants %in% RV$filtered_vars_list & 
            GDI_phred <= input$GDI_filter &
            pLI_gnomAD >= input$pLI_filter, "PASS", "FILTER")))    
    })
    
    genes_scores <- reactive({
        RV$filtered_genes_list <- unique(genes_df()$Gene[genes_df()$Class == "PASS"]) 
        as.data.frame(RV$genes_scores %>% mutate(Class = ifelse(
            Gene %in% RV$filtered_genes_list, "PASS", "FILTER")))
    })
    
    filters_summ_genes <- reactive ({
        PASS_counts <- c(
        genes_df() %>% filter(GDI_phred <= input$GDI_filter) %>% nrow(),
        genes_df() %>% filter(pLI_gnomAD >= input$pLI_filter) %>% nrow() )
        
        filters_summ_genes <- data.frame(Filter=c("pLI gnomAD","GDI phred"), 
            PASS=PASS_counts, 
            FILTERED=(nrow(RV$genes_df)-PASS_counts) )
        filters_summ_genes <- gather(filters_summ_genes, key="Class", value="Count", PASS:FILTERED)   
    })
    
    filters_summ_vars <- reactive ({
        tot_vars <- RV$variants_df %>% select(Chr,Pos,Ref,Alt) %>% nrow()
        PASS_counts <- c(
            variants_df() %>% filter(d_score >= input$d_score_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter((Reg_type == "splicing" & SpliceAI_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing") %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(MaxPopAF <= input$MaxPopAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(Consequence %in% RV$accepted_consequence) %>% select(Chr,Pos,Ref,Alt) %>% nrow(),
            variants_df() %>% filter(cohortAF <= input$CohortAF_filter) %>% select(Chr,Pos,Ref,Alt) %>% nrow()
        )
        filters_summ_vars <- data.frame(Filter=c("d score","SpliceAI","MaxPop AF","Consequence","cohort AF"),
            PASS=PASS_counts, 
            FILTERED=(tot_vars-PASS_counts))
        filters_summ_vars <- gather(filters_summ_vars, key="Class", value="Count", PASS:FILTERED)
    })
    
    ##################
    ### Overview Tab
    ##################
    
    output$Total_affected <- renderText({
        paste0("\tTotal number of affected individuals in this pedigree: ", RV$total_affected)
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
                    max = max(RV$variants_df$SpliceAI_SpliceAI_max, na.rm = T),
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
                    max = max(RV$genes_df$pLI_gnomAD, na.rm = T),
                    value = 0,
                    step = 0.05)
    })
    
    output$GDI <- renderUI({
        startvalue = 0
        sliderInput("GDI_filter", "Max phred GDI:",
                    min = 0,
                    max = max(RV$genes_df$GDI_phred, na.rm = T),
                    value = max(RV$genes_df$GDI_phred, na.rm = T),
                    step = 0.1)
    })
    
    output$recessive <- renderUI({
        selectInput("recessive_filter", "N recessive:", choices = RV$values_segregation, multiple=FALSE)
    })
    output$dominant <- renderUI({
        selectInput("dominant_filter", "N dominant:", choices = RV$values_segregation, multiple=FALSE)
    })
    output$denovo <- renderUI({
        selectInput("denovo_filter", "N denovo:", choices = RV$values_segregation, multiple=FALSE)
    })
    output$comphet <- renderUI({
        selectInput("comphet_filter", "N comphet:", choices = RV$values_segregation, multiple=FALSE)
    })
    
    output$var_consequence <- renderUI({
        var_types <- sort(unique(RV$variants_df$Consequence))
        names(var_types) <- sort(unique(RV$variants_df$Consequence))
        var_types <- c("ALL" = "ALL", var_types)
        selectInput("Consequence", "Variant consequence:", choices = var_types, multiple=TRUE, selected="ALL")
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
            GADO_plot <- ggplot(genes_scores(), aes(x=Gado_zscore, fill=Class)) + geom_histogram(bins=100) + geom_vline(xintercept = Zscore, color="red") + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        } else {
            GADO_plot <- ggplot(genes_scores(), aes(x=Gado_zscore, fill=Class)) + geom_histogram(bins=100) + scale_fill_brewer(palette="Set3") + scale_y_sqrt() 
        }
        ggplotly( GADO_plot )
    })
    
    output$genesTable <- DT::renderDataTable(selection="single", {
        as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% select(-Class)) 
    })
    
    genesTable_proxy <- DT::dataTableProxy("genesTable")
    
    output$customGenesTable <- DT::renderDataTable(selection="single", {
        as.data.frame(genes_scores() %>% filter(Class == "PASS") %>% mutate(row_idx = row_number()) %>% filter(Gene %in% RV$custom_genes) %>% select(-Class))
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
                                gather(key="Variant",value = "varID", var1:var2) %>% 
                                inner_join(., variants_df()[variants_df()$Class == "PASS",], by=c("varID"="ID")) %>%
                                select(-Class.x,-Class.y))),
        zip_archive = paste0(input$CaseCode, ".results.zip")
    )
    
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
    
    output$geneInfo <- DT::renderDataTable(selection="none", {
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
        
        comphet_details <- comphet_df() %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", var1:var2)
        
        #as.data.frame(variants_df() %>% filter(Class == "PASS", Gene == gene_name))
        as.data.frame(variants_df() %>% filter(Class == "PASS", Gene == gene_name, ID %nin% comphet_details$varID))
    })
    
    output$comphetTable <- DT::renderDataTable(selection="none", options = list(scrollX = TRUE), {
        genedetail_tab <- as.data.frame(genes_scores() %>% filter(Class == "PASS"))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        
        comphet_details <- comphet_df() %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", var1:var2) %>% select(-Gene, -Variant) %>% distinct()
        validate(need(nrow(comphet_details)>0, 'No compound het variants'))
        comphet_details <- merge(comphet_details,variants_df()[variants_df()$Class == "PASS",], by.x="varID",by.y="ID")
        as.data.frame(comphet_details %>% select(-Class.x,-Class.y))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

# VARIANT EXPLORER
# Author: Edoardo Giacopuzzi
# Explore and filter annotated variants from V2

# Input are tables from the VARAN V2 filtering and annotation

library(shiny)
library(DT)
library(plotly)
library(ggplot2)
library(kinship2)

#Set data dir containing variants tables and look for files
data_dir <- "example_data"

files <- list.files(data_dir, pattern = "*.vars.tsv.gz")
samplesID <- gsub("V2\\.|\\.var2reg.vars.tsv.gz","",files)

#We can switch to read from a input table containing sampleID and ped locations...
#samples_df <- data.frame(ID=samplesID,file=files)

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
pathway_file <- "c2.cp.v7.0.symbols.gmt"
pathways <- list()
con  <- file(pathway_file, open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myvector <- (strsplit(oneLine, "\t"))
    pathways[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
}
close(con)

GO_file <- "c5.bp.v7.0.symbols.gmt"
GO_bp <- list()
con  <- file(GO_file, open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myvector <- (strsplit(oneLine, "\t"))
    GO_bp[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
}
close(con)

GO_file <- "c5.cc.v7.0.symbols.gmt"
GO_cc <- list()
con  <- file(GO_file, open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myvector <- (strsplit(oneLine, "\t"))
    GO_cc[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
}
close(con)

GO_file <- "c5.mf.v7.0.symbols.gmt"
GO_mf <- list()
con  <- file(GO_file, open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myvector <- (strsplit(oneLine, "\t"))
    GO_mf[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
}
close(con)

#######################
### Set environment ###
#######################

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
variants_bar_options <- c("PanelApp" = "PanelApp",
    "Region type" = "Reg_type",
    "Variant consequence" = "Consequence",
    "Variant type" = "VarType",
    "Chromosome" = "Chr")

##Set reactive objects
RV <- reactiveValues(
    variant_df = data.frame(), 
    genes_scores = data.frame(), 
    genes_df = data.frame(),
    comphet_df = data.frame(),
    load_status = character(),
    custom_genes = character())

######################
### USER INTERFACE ###
######################

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel(fluidRow(h1("Variant Explorer"), h3("Beta version now working only on 004Int001"))),

    sidebarLayout(
        ################    
        # Side-bar Panel
        ################
        sidebarPanel(
            #drop-down list of cases
            selectInput("CaseCode", h3("Case code:"), choices = sort(unique(samplesID))),
            
            h3("Disease:"),
            textOutput("Disease"),
            
            h3("Custom gene list"),
            textInput(inputId = "custom_file",label = "File:",value = "", placeholder = "gene list file"),
            actionButton(inputId = "Load_file", label = "Load file"),
            textOutput("Loading_result")
        ),
        
        ############
        # Main-panel
        ############
        mainPanel(
            tabsetPanel(id = "mainarea",
                tabPanel("Overview",
                    plotOutput("ped"),
                    plotOutput("Var_consequence_plot"),
                    plotOutput("Var_type_plot"),
                    plotOutput("Var_PopAF_plot"),
                    plotOutput("Gene_segregation_plot")
                ),
                    tabPanel("Filters",
                    br(),
                    actionButton(inputId = "Apply_filters", label = "Apply filters"),
                         
                    h3("Variant filters"),
                    uiOutput("d_score"),
                    uiOutput("spliceAI"),
                    uiOutput("Max_pop_AF"),
                    uiOutput("Cohort_AF"),
                    uiOutput("var_consequence"),
                    
                    hr(), 
                    h3("Gene filters"),
                    uiOutput("pLI"),
                    uiOutput("GDI"),
                    
                    hr(),
                    h3("Segregation filters"),
                    h4("Filters based on the number of affected individuals in which the variants follow the desired model"),
                    h4("Works with OR logic"),
                    textOutput("Total_affected"),
                    br(),
                    fillRow(flex=4,
                            uiOutput("recessive"),
                            uiOutput("dominant"),
                            uiOutput("denovo"),
                            uiOutput("comphet"))
                ),
                tabPanel("Filter Explorer",
                    h3("Genes filter evaluation"),
                    fluidRow(selectInput("genes_X_axis", h4("X axis:"), choices = genes_axes_options, selected = "Gado_zscore"),
                        selectInput("genes_Y_axis", h4("Y axis:"), choices = genes_axes_options, selected = "pLI_gnomAD")
                    ),
                    plotlyOutput("genes_scatter",width = "100%"),
                    br(),
                    
                    h3("Variants filter evaluation"),
                    fluidRow(selectInput("variants_X_axis", h4("X axis:"), choices = variants_axes_options, selected = "MaxPopAF"),
                             selectInput("variants_Y_axis", h4("Y axis:"), choices = variants_axes_options, selected = "DANN_DANN")
                    ),
                    plotlyOutput("variants_scatter", width = "100%"),
                    br(),
                    selectInput("variants_bar_axis", h4("Category:"), choices = variants_bar_options, selected = "Reg_type"),
                    plotlyOutput("variants_bar",width = "100%")
                    
                ),
                tabPanel("Filter Results",
                    plotlyOutput("GADO_rank"),
                    #plotlyOutput("GADO_rank_filt"),
                    DT::dataTableOutput("genesTable"),
                    h3("Custom list genes"),
                    DT::dataTableOutput("customGenesTable")
                ),
                tabPanel("Gene details View",
                    fluidRow(
                        column(2, h4("Gene symbol:")) , column(9, offset = 1, textOutput("Gene_symbol"))
                    ),
                    fluidRow(
                        column(2, h4("Gene name:")) , column(9, offset = 1, textOutput("Gene_name"))
                    ),
                    fluidRow(
                        column(2, h4("Link to GTeX:")) , column(9, offset = 1, uiOutput("GTeX_link"))
                    ),
                    br(),
                    h3("Scores for the selected gene"),
                    DT::dataTableOutput("geneDetail"),
                    h3("Filtered variants for the selected gene"),
                    DT::dataTableOutput("variantsTable"),
                    h3("Filtered compound hets for the selected gene"),
                    DT::dataTableOutput("comphetTable"),
                    hr(),
                    h3("Further details on the gene"),
                    DT::dataTableOutput("geneInfo"),
                )
            )
        )
    )
)

########################
### SERVER FUNCTIONS ###
########################
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$CaseCode, {
        RV$load_status == ""
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
        RV$variants_df$Class <- "PASS"
        
        RV$genes_df <- read.table(genes_file, header=T, sep="\t", stringsAsFactors = F)
        RV$genes_df <- as.data.frame(RV$genes_df %>% separate_rows(Variants))
        RV$genes_df$Class <- "PASS"

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
    
    observeEvent(input$Apply_filters, {
        if (input$Consequence == "ALL") { 
            accepted_consequence <- sort(unique(RV$variants_df$Consequence))    
        } else {
            accepted_consequence <- input$Consequence
        }
        
        segregating_vars_list <- RV$segregation_df %>% filter(
            recessive >= input$recessive_filter | 
            dominant >= input$dominant_filter |
            deNovo >= input$denovo_filter |
            comphet >= input$comphet_filter) %>% select(ID)
        segregating_vars_list <- segregating_vars_list$ID
            
        RV$variants_df <- as.data.frame(RV$variants_df %>% mutate(Class = ifelse(
            ID %in% segregating_vars_list &
            d_score >= input$d_score_filter & 
            Consequence %in% accepted_consequence &
            MaxPopAF <= input$MaxPopAF_filter &
            cohortAF <= input$CohortAF_filter &
            ((Reg_type == "splicing" & SpliceAI_SpliceAI_max >= input$spliceAI_filter) | Reg_type != "splicing"),"PASS","FILTER")))
        filtered_vars_list <- RV$variants_df$ID[RV$variants_df$Class == "PASS"]
        
        RV$comphet_df <- as.data.frame(RV$comphet_df %>% mutate(Class = ifelse(
            var1 %in% filtered_vars_list & 
            var2 %in% filtered_vars_list &
            ID %in% segregating_vars_list, "PASS", "FILTER")))
        filtered_vars_list <- unique(c(filtered_vars_list, RV$comphet_df$ID[RV$comphet_df$Class == "PASS"]))
        
        RV$genes_df <- as.data.frame(RV$genes_df %>% mutate(Class = ifelse(
            Variants %in% filtered_vars_list & 
            GDI_phred <= input$GDI_filter &
            pLI_gnomAD >= input$pLI_filter, "PASS", "FILTER")))
        filtered_genes_list <- unique(RV$genes_df$Gene[RV$genes_df$Class == "PASS"])
        
        RV$genes_scores <- as.data.frame(RV$genes_scores %>% mutate(Class = ifelse(
            Gene %in% filtered_genes_list, "PASS", "FILTER")))
    })
    
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

    output$Disease <- renderText({
        paste0(HPOs$Disease[HPOs$X.CaseID == input$CaseCode])
    })
    
    output$Total_affected <- renderText({
        paste0("Total number of affected individuals in this pedigree: ", RV$total_affected)
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
        gene_segregation_df <- as.data.frame(RV$genes_df %>% select(Gene,inh_model) %>% distinct())
        ggplot(gene_segregation_df, aes(x=inh_model)) + geom_bar() + labs(y="genes count")
    })
    
    GADO_rank_plot <- eventReactive(c(input$Apply_filters, input$CaseCode),{
        ggplot(RV$genes_scores, aes(x=Gado_zscore, fill=Class)) + geom_histogram(bins=100) + scale_fill_brewer(palette="Set3") + scale_y_sqrt()
    })
    
    output$GADO_rank <- renderPlotly({
        ggplotly(GADO_rank_plot())
    })
    
    genesTable_df <- eventReactive(c(input$Apply_filters, input$CaseCode),{
        as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
    })
    
    output$genesTable <- DT::renderDataTable(selection="single", {
        genesTable_df() 
    })
    
    genesTable_proxy <- DT::dataTableProxy("genesTable")
    
    customGenesTable_df <- eventReactive(c(input$Apply_filters,input$Load_file),{
        as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% mutate(row_idx = row_number()) %>% filter(Gene %in% RV$custom_genes) %>% select(-Class))
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
    
    #output$GADO_rank_filt <- renderPlotly({
    #    validate(need(nrow(RV$filtered_genes_scoress) > 0, 'no filters selected'))
    #    ggplotly( 
    #        ggplot(RV$filtered_genes_plot_df, aes(x=Gado_zscore, fill=class)) + geom_histogram(bins=100) + scale_fill_brewer(palette="Set3") + scale_y_sqrt()
    #    )
    #})
    
    genes_scatter_plot <- eventReactive(c(input$Apply_filters,input$genes_X_axis,input$genes_Y_axis),{
        ggplot(RV$genes_scores, aes_string(x=input$genes_X_axis, y=input$genes_Y_axis, color="Class", label="Gene")) + geom_point() 
    })
    
    output$genes_scatter <- renderPlotly({
        ggplotly( genes_scatter_plot() )
    })
    
    variants_scatter_plot <- eventReactive(c(input$Apply_filters,input$variants_X_axis,input$variants_Y_axis),{
        ggplot(RV$variants_df, aes_string(x=input$variants_X_axis, y=input$variants_Y_axis, color="Class", label="Gene")) + geom_point() 
    })
    
    output$variants_scatter <- renderPlotly({
        ggplotly( variants_scatter_plot() )
    })
    
    variants_bar_plot <- eventReactive(c(input$Apply_filters, input$variants_bar_axis),{
        ggplot(RV$variants_df, aes_string(x=as.factor(input$variants_bar_axis), fill="Class")) + geom_bar() 
    })
    
    output$variants_bar <- renderPlotly({
        ggplotly( variants_bar_plot() )
    })
    
    output$Gene_symbol <- renderText({
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        gene_name
    })
    
    output$Gene_name <- renderText({
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        genes_info[genes_info$symbol == gene_name, "name"]
    })
    
    output$GTeX_link <- renderUI({
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        ensg_id <- genes_info[genes_info$symbol == gene_name, "ensembl_gene_id"]
        tags$a(href=paste0("https://gtexportal.org/home/gene/", gene_name), paste0(gene_name, "(", ensg_id, ")"), target="_blank")
    })
    
    output$geneDetail <- DT::renderDataTable(selection="none", {
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
        genedetail_tab[input$genesTable_rows_selected,]
    })
    
    output$geneInfo <- DT::renderDataTable(selection="none", {
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS") %>% select(-Class))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        selected <- list()
        pathways_list <- grep(gene_name, pathways)
        selected[["pathways"]] <- names(pathways)[pathways_list]
        
        GO_bp_list <- grep(gene_name, GO_bp)
        selected[["GO_BP"]] <- names(GO_bp)[GO_bp_list]
        
        GO_mf_list <- grep(gene_name, GO_mf)
        selected[["GO_MF"]] <- names(GO_mf)[GO_mf_list]
        
        GO_cc_list <- grep(gene_name, GO_cc)
        selected[["GO_CC"]] <- names(GO_cc)[GO_cc_list]
        
        Additional_info <- data.frame(lapply(selected, "length<-", max(lengths(selected))), stringsAsFactors = F)
    })
    
    output$variantsTable <- DT::renderDataTable(selection="none", {
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS"))
        #validate(need(nrow(genedetail_tab)>0, 'No variants passing the filters'))
        
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        
        as.data.frame(RV$variants_df %>% filter(Class == "PASS", Gene == gene_name))
    })
    
    output$comphetTable <- DT::renderDataTable(selection="none", {
        genedetail_tab <- as.data.frame(RV$genes_scores %>% filter(Class == "PASS"))
        gene_name <- genedetail_tab[input$genesTable_rows_selected, "Gene"]
        validate(need(gene_name != "", 'No gene selected'))
        comphet_details <- RV$comphet_df %>% filter(Class == "PASS", Gene == gene_name) %>% gather(key="Variant",value = "varID", var1:var2)
        validate(need(nrow(comphet_details)>0, 'No compound het variants'))
        comphet_details <- merge(comphet_details,RV$variants_df, by.x="varID",by.y="ID")
        as.data.frame(comphet_details %>% select(-Gene,-Variant,-varID))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

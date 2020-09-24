query_results_UI <- function(id,boxes=FALSE) {
  ns <- NS(id)
  region_select <- tagList(selectInput(ns("region_select"), 
                                       "Regulatory region",
                                       choices = c("global","by tissue"),
                                       selected = "global", 
                                       multiple = FALSE) )
  
  plot_output <- tagList(withSpinner(plotlyOutput(ns("region_plot"))))
  
  if (boxes == TRUE) {
    details <- tagList(
      box(title = "Reg regions details", width = 12, status = "primary", solidHeader = T, 
          collapsed = T, collapsible = T,
          DT::dataTableOutput(ns("regions"))),
      box(title = "Controlled genes details", width = 12, status = "primary", solidHeader = T, 
          collapsed = T, collapsible = T,
          DT::dataTableOutput(ns("genes_details"))),
      box(title = "Phenotypes details", width = 12, status = "primary", solidHeader = T, 
          collapsed = T, collapsible = T,
          DT::dataTableOutput(ns("pheno_details"))),
      box(title = "TFBS details", width = 12, status = "primary", solidHeader = T, 
          collapsed = T, collapsible = T,
          DT::dataTableOutput(ns("TFBS_details")))
      )
  } else {
    details <- tagList(
      DT::dataTableOutput(ns("regions")),
      br(),
      DT::dataTableOutput(ns("genes_details")),
      br(),
      DT::dataTableOutput(ns("pheno_details")),
      br(),
      DT::dataTableOutput(ns("TFBS_details"))
      )
  }
  
  output_UI <- tagList(
    fluidRow(
      column(4,region_select),
      column(8,plot_output) ),
    br(),
    details )

  return(output_UI)
} 

SQlite_query <- function(input, output, session, db, var_position, regions) {

  conn <- dbConnect(RSQLite::SQLite(), db)
  regions_query <- paste("\'", regions,"\'", collapse=",", sep="")

  regions_df <- reactive({
    query = paste0("SELECT r.regionID, r.chromosome, r.start, r.stop, \
          r.type, r.std_type, r.DB_source, r.PhyloP100_median, r.constrain_pct, \
          GROUP_CONCAT(DISTINCT g.gene_symbol) AS controlled_genes, \
          r.closestGene_symbol, r.closestGene_ensg, r.closestGene_dist, \
          GROUP_CONCAT(DISTINCT t.cell_or_tissue) AS cell_or_tissue, \
          GROUP_CONCAT(DISTINCT m.method) AS detection_method, \
          GROUP_CONCAT(DISTINCT p.phenotype) AS phenotypes \
          FROM GRCh38_Regions AS r \
          LEFT JOIN tissues AS t ON r.regionID = t.regionID \
          LEFT JOIN genes AS g ON r.regionID = g.regionID \
          LEFT JOIN methods AS m ON r.regionID = m.regionID \
          LEFT JOIN phenos AS p ON r.regionID = p.regionID \
          WHERE r.regionID IN (",regions_query, ") \
          GROUP BY r.regionID")
    df <- dbGetQuery(conn, query)
    return(df)
  })
  
  genes_details_df <- reactive({
    query = paste0("SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
    g.gene_symbol, m.method, t.cell_or_tissue \
    FROM GRCh38_Regions AS r \
    LEFT JOIN genes AS g ON r.regionID = g.regionID \
    LEFT JOIN tissues AS t ON g.interactionID = t.regionID \
    LEFT JOIN methods AS m ON g.interactionID = m.regionID \
    WHERE r.regionID IN (",regions_query, ")")
    df <- dbGetQuery(conn, query)
    return(df)
  })

  pheno_details_df <- reactive({
    query = paste0("SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
    p.phenotype, p.method, p.DB_source \
    FROM GRCh38_Regions AS r \
    LEFT JOIN phenos AS p ON r.regionID = p.regionID \
    WHERE r.regionID IN (",regions_query, ")")
    df <- dbGetQuery(conn, query)
    return(df)
  })
  
  TFBS_details_df <- reactive({
    query = paste0("SELECT r.regionID AS regionID, \
        t.chromosome, t.start, t.stop, t.name, \
        tt.cell_or_tissue AS cell_or_tissue \
        FROM GRCh38_Regions AS r \
        INNER JOIN GRCh38_regionID_to_TFBS AS link ON r.regionID = link.regionID \
        INNER JOIN GRCh38_TFBS AS t ON link.link_ID = t.regionID \
        INNER JOIN GRCh38_TFBS_tissues AS tt ON t.regionID = tt.regionID \
        WHERE r.regionID IN (",regions_query, ")")
    df <- dbGetQuery(conn, query)
    return(df)
  })
  
  output$regions <- DT::renderDataTable({
    datatable(regions_df(), 
              selection="none",
              options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200)))
  })
  
  output$genes_details <- DT::renderDataTable({
    datatable(genes_details_df(), 
              selection="none",
              options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200)))
  })
  
  output$pheno_details <- DT::renderDataTable({
    datatable(pheno_details_df(), 
              selection="none",
              options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200)))
  })
  
  output$TFBS_details <- DT::renderDataTable({
    datatable(TFBS_details_df(), 
              selection="none",
              options = list(scrollX = TRUE, pageLength = 20, lengthMenu = c(20, 50, 100, 200)))
  })
  
  output$region_plot <- renderPlotly({
    req(input$region_select)
    regions_plot <- regions_df()[,c("chromosome","start","stop","regionID","cell_or_tissue")]
    colnames(regions_plot)[4] <- "name"
    plot_df <- rbind(TFBS_details_df() %>% select(chromosome,start,stop,name,cell_or_tissue) %>% mutate(group=paste0("TFBS - ",name)),
                     regions_plot %>% mutate(group="region"))
    if (input$region_select == "by tissue") {
      p <- ggplot(plot_df, aes(ymin=start, ymax=stop, x=group, colour=cell_or_tissue, label=name)) +
        geom_linerange(size=1, position=position_dodge(1)) +
        geom_text(data=plot_df %>% filter(group=="region"), aes(y=start+(stop-start)/2), position=position_dodge(1), color="black") +
        geom_hline(yintercept = var_position, linetype="dashed") +
        labs(y=paste0(unique(plot_df$chromosome), " genomic position"), x="") +
        coord_flip() +
        theme(axis.text.y = element_text(size=10), axis.text.x=element_text(angle=45, hjust=1), legend.position = "none")
      tooltip <- c("label","colour")
    } else if (input$region_select == "global") {
      p <- ggplot(plot_df, aes(ymin=start, ymax=stop, x=group, label=name)) +
        geom_linerange(size=2, position=position_dodge(0.6), color="deepskyblue3") +
        geom_text(data=plot_df %>% filter(group=="region"), aes(y=start+(stop-start)/2), position=position_dodge(0.6), color="black") +
        geom_hline(yintercept = var_position, linetype="dashed") +
        labs(y=paste0(unique(plot_df$chromosome), " genomic position"), x="") +
        coord_flip() +
        theme(axis.text.y = element_text(size=10), axis.text.x=element_text(angle=45, hjust=1), legend.position = "none")
      tooltip <- c("label","ymin","ymax")
    }
     ggplotly(p, tooltip=tooltip)
  
      
  })
  
  #RSQLite::dbDisconnect(conn)
}
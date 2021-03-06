##########################
### SEGREGATION MODULE ###
##########################

operations <- c("equal" = "==", "greater than" = ">=", "less than" = "<=")

RV <- reactiveValues()

#This function set the expression on segregation column
#This is set to >= for aff and <= for unaff columns when a number is provided
#When NOT_ALLOWED is provided as value the filter is set to == 0 for the corresponding column
selectFilter <- function(var,operator,value) {
  switch(operator,
         "==" = { myexpr <- expr(!!as.name(var) == !!as.numeric(value)) },
         ">" = { myexpr <- expr(!!as.name(var) > !!as.numeric(value)) },
         "<" = { myexpr <- expr(!!as.name(var) < !!as.numeric(value)) },
         ">=" = { myexpr <- expr(!!as.name(var) >= !!as.numeric(value)) },
         "<=" = { myexpr <- expr(!!as.name(var) <= !!as.numeric(value)) }
  )

  return(myexpr)
}

segregationUI <- function(id, choices_affected, choices_unaffected) {
  ns <- NS(id)
  
  #A list of settings for seg models as a named list
  #Each model contain a 5 values in the following order
  #het_affected, het_unaffected, hom_affected, hom_unaffected, comphet_affected, sup_dnm
  if (max(choices_unaffected) >= 2) {value_unaffected <- 2} else {value_unaffected <- max(choices_unaffected)}
  RV$segregation_presets_values <- list(
    perfect_recessive = c(0,value_unaffected,max(choices_affected),0,max(choices_affected),0),
    perfect_dominant = c(max(choices_affected),0,0,0,0,0),
    dnm_permissive = c(max(choices_affected),0,0,0,0,0),
    dnm_highconf = c(max(choices_affected),0,0,0,0,max(choices_affected)),
    free = c(0,0,0,0,0,0)
  )
  
  segregation_presets <- tagList(
    selectInput(ns("segregation_preset"), "Segregation model", choices = names(RV$segregation_presets_values), multiple=F, selected="free")
  )
  
  homozygous_controls <- tagList(
    h4("Homozygous"),
    selectInput(ns("hom_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("hom_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "=="),
    selectInput(ns("hom_unaffected"), "Unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0"),
    selectInput(ns("hom_unaffected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "==")
  )
  
  heterozygous_controls <- tagList(
    h4("Heterozygous"),
    selectInput(ns("het_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("het_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "=="),
    selectInput(ns("het_unaffected"), "Unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0"),
    selectInput(ns("het_unaffected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "==")
  )
  
  comphet_controls <- tagList(
    h4("Compound hets"),
    selectInput(ns("comphet_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("comphet_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "==")
  )
  
  dnm_controls <- tagList(
    h4("High-confident denovo"),
    selectInput(ns("dnm_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "0"),
    selectInput(ns("dnm_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = ">=")
  
  )
  
  callModule(updateSegregationPresets, "segregation")
  
  output <- tagList(
    fluidRow(column(6, segregation_presets)),
    hr(),
    fluidRow(
              column(4,homozygous_controls),
              column(4,heterozygous_controls),
              column(4,comphet_controls) ),
              fluidRow(column(4), column(4, align="center", dnm_controls), column(4))
            )
  #output <- tagList(fluidRow(output))
}

updateSegregationPresets <- function(input, output, session) {
  observeEvent(input$segregation_preset,{
    #het_affected, het_unaffected, hom_affected, hom_unaffected, comphet_affected, sup_dnm
    values <- RV$segregation_presets_values[[input$segregation_preset]]
    updateSelectInput(session, inputId = "het_affected",selected = values[1])
    updateSelectInput(session, inputId = "het_unaffected",selected = values[2])
    updateSelectInput(session, inputId = "hom_affected",selected = values[3])
    updateSelectInput(session, inputId = "hom_unaffected",selected = values[4])
    updateSelectInput(session, inputId = "comphet_affected",selected = values[5])
    updateSelectInput(session, inputId = "dnm_affected",selected = values[6])
    
    updateSelectInput(session, inputId = "het_affected_op",selected = "==")
    updateSelectInput(session, inputId = "het_unaffected_op",selected = "==")
    updateSelectInput(session, inputId = "hom_affected_op",selected = "==")
    updateSelectInput(session, inputId = "hom_unaffected_op",selected = "==")
    updateSelectInput(session, inputId = "comphet_affected_op",selected = "==")
    updateSelectInput(session, inputId = "dnm_affected_op",selected = ">=")
  })
}

#cols_names is a named vector of 5 elements describing columns names for segregation counts
#expected names are: het_affected, het_unaffected, hom_affected, hom_unaffected, comphet_affected
#filter structure: comphet OR (hom filters AND het filters)
segregationModule <- function(input, output, session, segregation_df, cols_names) {
  seg_vars_list <- NULL
  
  #First evalute comphet (comphet must have all het and hom cols equal to zero)
  het_expr <- expr(!!selectFilter(cols_names["het_affected"], "==", "0") &
                    !!selectFilter(cols_names["het_unaffected"], "==", "0"))
  hom_expr <- expr(!!selectFilter(cols_names["hom_affected"], "==", "0") &
                    !!selectFilter(cols_names["hom_unaffected"], "==", "0"))
  comphet_expr <- selectFilter(cols_names["comphet_affected"], input$comphet_affected_op, input$comphet_affected)
  segregation_expr <- expr(!!het_expr & !!hom_expr & !!comphet_expr)
  seg_vars <- as.data.frame(segregation_df %>% filter(!!segregation_expr))
  seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  
  #Then evaluate het and hom calls with AND logic (comphet col must be zero)
  het_expr <- expr(!!selectFilter(cols_names["het_affected"], input$het_affected_op, input$het_affected) &
                    !!selectFilter(cols_names["het_unaffected"], input$het_unaffected_op, input$het_unaffected))
  hom_expr <- expr(!!selectFilter(cols_names["hom_affected"], input$hom_affected_op, input$hom_affected) &
                    !!selectFilter(cols_names["hom_unaffected"], input$hom_unaffected_op, input$hom_unaffected))
  dnm_expr <- selectFilter(cols_names["dnm_affected"], input$dnm_affected_op, input$dnm_affected)
  comphet_expr <- selectFilter(cols_names["comphet_affected"], "==", "0")
  segregation_expr <- expr(!!het_expr & !!hom_expr & !!comphet_expr & !!dnm_expr)
  seg_vars <- as.data.frame(segregation_df %>% filter(!!segregation_expr))
  seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  
  return(seg_vars_list)
}

getDF_segregation <- function(input, output, session) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)

    for (n in ctrl_names) {
      newline <- c("segregation", n, input[[n]])
      filters_df[nrow(filters_df)+1,] <- newline
    }
  filters_df <- filters_df[order(filters_df$filter),]
  return(filters_df)
}

getJSON_segregation <- function(input, output, session) {
  filters_json <- list()
  ctrl_names <- names(input)
  
  for (n in ctrl_names) {
    filters_json[[n]] <- input[[n]]
  }
  
  return(filters_json) 
}

loadSettings_segregation <- function(input, output, session, filters_values) {
  for (v in names(filters_values)) {
    if (v %in% names(input)) {
      updateSelectInput(session, inputId = v, selected = filters_values[[v]])
    } else {
      message("LOAD WARNING: ",v, " control not found in the segregation module")
    }
  }
}
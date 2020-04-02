##########################
### SEGREGATION MODULE ###
##########################

#This function set the expression on segregation column
#This is set to >= for aff and <= for unaff columns when a number is provided
#When NOT_ALLOWED is provided as value the filter is set to == 0 for the corresponding column
selectFilter <- function(var,operator,value) {
  if (value == 0) {
    myexpr <- quo((!!(as.name(var))) == !!(as.numeric(value)))
  } else {
    switch(operator,
           ">" = { myexpr <- quo((!!(as.name(var))) > !!(as.numeric(value))) },
           "<" = { myexpr <- quo((!!(as.name(var))) < !!(as.numeric(value))) },
           ">=" = { myexpr <- quo((!!(as.name(var))) >= !!(as.numeric(value))) },
           "<=" = { myexpr <- quo((!!(as.name(var))) <= !!(as.numeric(value))) }
    )
  }
  return(myexpr)
}


segregationUI <- function(id, choices_affected, choices_unaffected) {
  ns <- NS(id)
  
  homozygous_controls <- tagList(
    h4("Homozygous"),
    selectInput(ns("hom_affected"), "Min affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("hom_unaffected"), "Max unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0")
  )
  
  heterozygous_controls <- tagList(
    h4("Heterozygous"),
    selectInput(ns("het_affected"), "Min affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("het_unaffected"), "Max unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0")
  )
  
  comphet_controls <- tagList(
    h4("Compound hets"),
    selectInput(ns("comphet_affected"), "Min affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
  )
  
  output <- tagList(
    column(4,homozygous_controls),
    column(4,heterozygous_controls),
    column(4,comphet_controls)
  )
  output <- tagList(fluidRow(output))

}

#cols_names is a named vector of 5 elements describing columns names for segregation counts
#expected names are: het_affected, het_unaffected, hom_affected, hom_unaffected, comphet_affected
#filter structure: comphet OR (hom filters AND het filters)
segregationModule <- function(input, output, session, segregation_df, cols_names) {
  seg_vars_list <- NULL
  
  #First evalute comphet, is they are requested (input$comphet_affected > 0)
  if (input$comphet_affected > 0) {
    comphet_expr <- selectFilter(cols_names["comphet_affected"], ">=", input$comphet_affected)
    seg_vars <- as.data.frame(segregation_df %>% filter(!!comphet_expr))
    seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  }
  
  #Then evaluate het and hom calls with AND logic
  het_expr <- quo(!!selectFilter(cols_names["het_affected"], ">=", input$het_affected) &
                    !!selectFilter(cols_names["het_unaffected"], "<=", input$het_unaffected))
  hom_expr <- quo(!!selectFilter(cols_names["hom_affected"], ">=", input$hom_affected) &
                    !!selectFilter(cols_names["hom_unaffected"], "<=", input$hom_unaffected))
  segregation_expr <- quo(!!het_expr & !!hom_expr)
  
  seg_vars <- as.data.frame(segregation_df %>% filter(!!segregation_expr))
  seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  
  return(seg_vars_list)
}
##########################
### SEGREGATION MODULE ###
##########################

operations <- c("equal" = "==", "greater than" = ">=", "less than" = "<=")

#This function set the expression on segregation column
#This is set to >= for aff and <= for unaff columns when a number is provided
#When NOT_ALLOWED is provided as value the filter is set to == 0 for the corresponding column
selectFilter <- function(var,operator,value) {
  switch(operator,
         "==" = { myexpr <- quo((!!(as.name(var))) == !!(as.numeric(value))) },
         ">" = { myexpr <- quo((!!(as.name(var))) > !!(as.numeric(value))) },
         "<" = { myexpr <- quo((!!(as.name(var))) < !!(as.numeric(value))) },
         ">=" = { myexpr <- quo((!!(as.name(var))) >= !!(as.numeric(value))) },
         "<=" = { myexpr <- quo((!!(as.name(var))) <= !!(as.numeric(value))) }
  )

  return(myexpr)
}

segregationUI <- function(id, choices_affected, choices_unaffected) {
  ns <- NS(id)
  
  homozygous_controls <- tagList(
    h4("Homozygous"),
    selectInput(ns("hom_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("hom_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "equal"),
    selectInput(ns("hom_unaffected"), "Unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0"),
    selectInput(ns("hom_unaffected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "equal")
  )
  
  heterozygous_controls <- tagList(
    h4("Heterozygous"),
    selectInput(ns("het_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("het_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "equal"),
    selectInput(ns("het_unaffected"), "Unaffected carriers:", choices = choices_unaffected, multiple=FALSE, selected = "0"),
    selectInput(ns("het_unaffected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "equal")
  )
  
  comphet_controls <- tagList(
    h4("Compound hets"),
    selectInput(ns("comphet_affected"), "Affected carriers:", choices = choices_affected, multiple=FALSE, selected = "1"),
    selectInput(ns("comphet_affected_op"), "Operator:", choices = operations, multiple=FALSE, selected = "equal")
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
  
  #First evalute comphet
  comphet_expr <- selectFilter(cols_names["comphet_affected"], input$comphet_affected_op, input$comphet_affected)
  seg_vars <- as.data.frame(segregation_df %>% filter(!!comphet_expr))
  seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  
  #Then evaluate het and hom calls with AND logic
  het_expr <- quo(!!selectFilter(cols_names["het_affected"], input$het_affected_op, input$het_affected) &
                    !!selectFilter(cols_names["het_unaffected"], input$het_unaffected_op, input$het_unaffected))
  hom_expr <- quo(!!selectFilter(cols_names["hom_affected"], input$hom_affected_op, input$hom_affected) &
                    !!selectFilter(cols_names["hom_unaffected"], input$hom_unaffected_op, input$hom_unaffected))
  comphet_expr <- selectFilter(cols_names["comphet_affected"], "==", "0")
  segregation_expr <- quo(!!het_expr & !!hom_expr & !!comphet_expr)
  
  seg_vars <- as.data.frame(segregation_df %>% filter(!!segregation_expr))
  seg_vars_list <- c(seg_vars_list,seg_vars$rec_id)
  
  return(seg_vars_list)
}
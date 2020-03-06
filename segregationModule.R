##########################
### SEGREGATION MODULE ###
##########################

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
  
  seg_vars <- as.data.frame(segregation_df %>% filter(
    !!as.symbol(cols_names["comphet_affected"]) >= input$comphet_affected |
    (!!as.symbol(cols_names["het_affected"]) >= input$het_affected & 
      !!as.symbol(cols_names["het_unaffected"]) <= input$het_unaffected &
      !!as.symbol(cols_names["hom_affected"]) >= input$hom_affected &
      !!as.symbol(cols_names["hom_unaffected"]) <= input$hom_unaffected)
    )
  )
  
  return(seg_vars$ID)
}
# MODULE FOR FILTERS MANAGMENT
# Allow to configure groups of filters
# return a quotation object containing AND concatened filters 

# variables: list with names corresponding to varID (name of variable/column to apply filter).
# Each element contain a sublist with names
#   label = string to be used as label for the control
#   selected = the value selected when control is initialized
#   values = a vector of possible choices or a vector with min,max,step for numeric scales
#   filter = string containing one of <, >, >=, <=, %in% (filter to be applied)
#   type = slider or multichoice

# UI can be configured passing ncols and variables arguments
# ncols set the number of columns to arrange the controls in fluidRow

selectFilter <- function(var,operator,value) {
  switch(operator,
    ">" = { myexpr <- quo((!!(as.name(var))) > !!(as.numeric(value))) },
    "<" = { myexpr <- quo((!!(as.name(var))) < !!(as.numeric(value))) },
    ">=" = { myexpr <- quo((!!(as.name(var))) >= !!(as.numeric(value))) },
    "<=" = { myexpr <- quo((!!(as.name(var))) <= !!(as.numeric(value))) },
    "%in%" = { myexpr <- quo((!!(as.name(var))) %in% !!value) }
  )
  return(myexpr)
}

filtersGroupUI <- function(id,variables,ncols) {
  ns <- NS(id)
  
  col_size <- 12/ncols
  UI_elements <- list()
  
  for (v in names(variables)) {
    switch(variables[[v]][["type"]],
      multichoice = {input_element <- tagList(selectInput(ns(v), variables[[v]][["label"]], choices = variables[[v]][["values"]], selected = variables[[v]][["selected"]], multiple = TRUE))},
      slider = {input_element <-tagList(sliderInput(ns(v), variables[[v]][["label"]], 
                                   min = variables[[v]][["values"]][1],
                                   max = variables[[v]][["values"]][2],
                                   value = variables[[v]][["selected"]],
                                   step = variables[[v]][["values"]][3]))}  
    )
    UI_elements[[v]] <- input_element
  }
  
  row_elements <- split(names(variables), ceiling(seq_along(variables)/ncols))
  out_rows <- NULL
  for (e in row_elements) {
    myrow <- NULL
    for (c in length(r)) {
      myrow <- tagList(myrow, UI_elements[[ row_elements[[e]][c] ]])
    }
    out_rows <- tagList(out_rows, myrow)
  }
  return(out_rows)
}

# Filters are constructed based on variables list
filtersGroupModule <- function(input, output, session, variables) {
  final_expr <- selectFilter(names(variables)[1],variables[[1]][["filter"]],input[[names(variables)[1]]])
  
  for (n in 2:length(variables)) {
    myexpr <- selectFilter(names(variables)[n],variables[[n]][["filter"]],input[[names(variables)[n]]])
    final_expr <- quo(!!final_expr & !!myexpr)
  }
  
  return(final_expr)
}
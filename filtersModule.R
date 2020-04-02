# MODULE FOR VARIANT FILTERS MANAGEMENT
# Allow to configure define groups of vars and apply specific filters to each group
# Suppose a variants table containing various annotations columns
# returns a quotation object containing concatened filters for selecting FILTERED VARS
# Practical use to obtain PASS variants is to negate this in dplyr filter(!(!!!returned_object))
# Returned filters are configured as follows:
# global_filters | (var in group1 & (group1 filters)) | (var in group1 & (group1 filters)) ...

# Filter configuration is defined in 2 json files
# definitions:  contains 3 groups for numerical, factors and binary columns
#               each groups define column names and the associated label to show.
#               for numerical columns the operator of filter is also provided like >,<
#               for factors columnd %in% is used
# settings:     allow to define filters groups and presets

# UI aspect can be configured passing ncols
# ncols set the number of columns to arrange the controls in fluidRow
# One separated column is created for each filter group defined in settings
# Two additional multichoice ar added per group 
#   - AND / OR logic
#   - Preset if present in settings json

selectFilter <- function(var,operator,value) {
  switch(operator,
    ">" = { myexpr <- quo((!!(as.name(var))) > !!(as.numeric(value))) },
    "<" = { myexpr <- quo((!!(as.name(var))) < !!(as.numeric(value))) },
    ">=" = { myexpr <- quo((!!(as.name(var))) >= !!(as.numeric(value))) },
    "<=" = { myexpr <- quo((!!(as.name(var))) <= !!(as.numeric(value))) }
#    "%in%" = { myexpr <- quo((!!(as.name(var))) %in% !!value) }
  )
  return(myexpr)
}

filtersGroupUI <- function(id,filters_def,filters_set,ncols) {
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
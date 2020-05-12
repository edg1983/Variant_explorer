########################
### GQ FILTER MODULE ###
########################

#Control UI composed by a drop-down box and a slider
#maxGQ = slider max value
#defaultGQ = slider value on load
GQfilterUI <- function(id, maxGQ, defaultGQ) {
  ns <- NS(id)
  
  tagList(
    h4("GQ filter"),
    fluidRow(
      column(4,selectInput(ns("GQ_samples"), "Apply to:", choices = c("ANY", "ALL", "AFFECTED"), multiple=FALSE, selected = "AFFECTED")),
      column(8,sliderInput(ns("GQ_value"), label = "Min GQ value", min = 0, max = maxGQ, step = 1, value = defaultGQ))
    )
  )
  
}

#Server module evaluate GQ columns and then apply filter
#variants_df = dataframe of variants
#affected_cols = indexes of affected GQ cols
#GQ_cols = indexes of GQ cols for all samples
#exclude_var_type = variant types (from var_type column) to be ignored when applying the filter
#This latter is useful right now since GQ filtering can not be done effectively on SVs, so we can ignore them
GQfilterModule <- function(input, output, session, variants_df, GQ_cols, affected_cols, exclude_var_type=NULL) {
  
  variants_df$aff_GQ <- rowSums(variants_df[,..affected_cols] >= input$GQ_value,na.rm = T)
  variants_df$all_GQ <- rowSums(variants_df[,..GQ_cols] >= input$GQ_value,na.rm = T)
  n_samples <- length(GQ_cols)
  n_affected <- length(affected_cols)
  
  switch(input$GQ_samples,
         "ANY" = { pass_vars <- variants_df %>% filter(all_GQ >= 1 | var_type %in% exclude_var_type) },
         "ALL" = { pass_vars <- variants_df %>% filter(all_GQ == n_samples | var_type %in% exclude_var_type) },
         "AFFECTED" = { pass_vars <- variants_df %>% filter(aff_GQ == n_affected | var_type %in% exclude_var_type) }
  )
  
  pass_vars_list <- pass_vars$rec_id
  
  return(pass_vars_list)
}

getDF_GQ <- function(input, output, session) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)
  
  for (n in ctrl_names) {
    newline <- c("GQ", n, input[[n]])
    filters_df[nrow(filters_df)+1,] <- newline
  }
  return(filters_df)
}

getJSON_GQ <- function(input, output, session) {
  filters_json <- list()
  ctrl_names <- names(input)

    for (n in ctrl_names) {
      filters_json[[n]] <- input[[n]]
    }

  return(filters_json) 
}

loadSettings_GQ <- function(input, output, session, filters_values) {
  if (!is.null(filters_values$GQ_samples)) {
    updateSelectInput(session, inputId = "GQ_samples", selected = filters_values$GQ_samples)
  } 
  if (!is.null(filters_values$GQ_value)) {
    updateSliderInput(session, inputId = "GQ_value", value = as.numeric(filters_values$GQ_value))
  }
}
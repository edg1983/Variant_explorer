###########################
### INTERSECTBED MODULE ###
###########################

#Arguments
#label: define checkbox label
#slider: generate a slider controlling the min size of regions for filtering
#samples: generate a selectInput of samples names from ranges
#ranges: a GRanges object or named list containing 1 Grange object per sample
#NOTE: 
#- WHEN SAMPLES is TRUE RANGES IS EXPECTED TO BE A LIST
#- WHEN SLIDER is TRUE, RANGES IS REQUIRED
bedcontrolUI <- function(id, label, ranges=NULL, slider=FALSE, samples=FALSE, selected_sample="AFFECTED_SHARED") {
  ns <- NS(id)
  
  checkbox <- tagList(checkboxInput(ns("bed_filter"), label = paste0("Select only vars in ", label), value = FALSE))
  filter_control <- NULL
  sample_control <- NULL
  
  if (samples) {
    sample_control <- tagList(
      selectInput(ns("sample_select"), label = "Select sample:", choices = names(ranges), multiple = FALSE, selected = selected_sample)
    )
  }
  
  if (slider) {
    filter_control <-  uiOutput(ns("filter_control"))
  }
  
  output <- tagList(
    fluidRow(
      column(4, checkbox),
      column(4, sample_control),
      column(4, filter_control)
    ) )
    
}

#Configure the slider dynamically based on selected sample
#scale_value: value to scale the slider
configureSlider <- function(input, output, session, ranges=NULL, scale_value=1) {
  output$filter_control <- renderUI ({
    if (!is.null(input$sample_select)) {
      max_value <- max(ranges[[input$sample_select]]$value)
    } else {
      max_value <- max(ranges$value)
    }
    sliderInput( session$ns("bed_limit"),
      label = "ROH min dimension (kb)",
      min = 0,
      max = max_value / scale_value,
      step = 10,
      value = 250 )
  })
}

#The function will perform intersection between variants ranges and bed ranges
#Return a list of vars IDs located within bed regions
#variants_ranges and bed_ranges are GenomicRanges objects
#variants_ranges must have an ID column, containing the var IDs
#multiplier can be passed to multiply the input value from input control 
#Useful if values in bed are large and you want to scale them inte UI control
bedfilterModule <- function(input, output, session, variants_ranges, bed_ranges, multiplier=1) {
  #if (!is.null(input$bed_filter)) {
    if (input$bed_filter == TRUE) {
      if (!is.null(input$sample_select)) { bed_ranges <- bed_ranges[[input$sample_select]]}
      if (!is.null(input$bed_limit)) { 
        threshold <- input$bed_limit * multiplier
        bed_ranges <- bed_ranges[bed_ranges$value >= threshold,]
      }
      
      mcols(bed_ranges)$value <- NULL
      xy <- c(variants_ranges, bed_ranges)
      r <- reduce(xy, with.revmap=TRUE)
      revmap <- mcols(r)$revmap
      r_values <- extractList(mcols(xy)$ID, revmap)
      r_values <- r_values[sapply(r_values, function(x) sum(grepl("ROH", x)) > 0)]
      
      varID_list <- unique(unlist(r_values))
      ROH_ids <- grep("ROH", varID_list)
      varID_list <- varID_list[-ROH_ids]
    } else {
      varID_list <- variants_ranges$ID
    } 
  #} else {
  #  varID_list <- variants_ranges$ID
  #}
  return(varID_list)
}

getDF_regions <- function(input, output, session, regions_label) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)
  
  for (n in ctrl_names) {
    newline <- c(regions_label, n, input[[n]])
    filters_df[nrow(filters_df)+1,] <- newline
  }
  filters_df <- filters_df[order(filters_df$filter),]
  return(filters_df)
}

getJSON_regions <- function(input, output, session) {
  filters_json <- list()
  ctrl_names <- names(input)
  for (n in ctrl_names) {
    filters_json[[n]] <- input[[n]]
  }
  return(filters_json) 
}

loadSettings_regions <- function(input, output, session, filters_values) {
    if (!is.null(filters_values$bed_filter)) {
      updateCheckboxInput(session, inputId = "bed_filter", value = as.logical(filters_values$bed_filter))
    } 
    if (!is.null(filters_values$bed_limit)) {
      updateSliderInput(session, inputId = "bed_limit", value = as.numeric(filters_values$bed_limit))
    }
    if (!is.null(filters_values$sample_select)) {
      updateSelectInput(session, inputId = "sample_select", value = filters_values$sample_select)
    }
}
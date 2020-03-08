###########################
### INTERSECTBED MODULE ###
###########################

#Use label to define checkbox label
#You can set filter config to a named vector to generate a slider control
#Named vector must contain
#label = label for the slider
#max = max value
#min = min value
#step = step of slider
#value = starting value
bedcontrolUI <- function(id, label, filter_config=FALSE) {
  ns <- NS(id)
  
  checkbox <- tagList(checkboxInput(ns("bed_filter"), label = label, value = FALSE))
  filter_control <- NULL
  
  if (filter_config != FALSE){
    filter_control <- tagList(sliderInput(
      ns("bed_limit"),
      label = filter_config[["label"]],
      min = as.numeric(filter_config[["min"]]),
      max = as.numeric(filter_config[["max"]]),
      step = as.numeric(filter_config[["step"]]),
      value = as.numeric(filter_config[["value"]]),
    ))
  }
  
  output <- tagList(
    fluidRow(
      column(6, checkbox),
      column(6, filter_control)
    ) )
    
}

#The function will perform intersection between variants ranges and bed ranges
#Return a list of vars IDs located within bed regions
#variants_ranges and bed_ranges are GenomicRanges objects
#variants_ranges must have an ID column, containing the var IDs
#multiplier can be passed to multiply the input value from input control 
#Useful if values in bed are large and you want to scale them inte UI control
bedfilterModule <- function(input, output, session, variants_ranges, bed_ranges, multiplier=1) {
  
  if (!is.null(input$bed_limit)) { 
    threshold <- input$bed_limit * multiplier
    bed_ranges <- bed_ranges[bed_ranges$value >= threshold,]
  }
  
  mcols(bed_ranges)$value <- NULL
  xy <- c(variants_ranges, bed_ranges)
  r <- reduce(xy, with.revmap=TRUE)
  revmap <- mcols(r)$revmap
  r_values <- extractList(mcols(xy)$ID, revmap)
  
  varID_list <- unique(unlist(r_values))
  ROH_ids <- grep("ROH", ID_list)
  varID_list <- ID_list[-ROH_ids]

  return(varID_list)
}
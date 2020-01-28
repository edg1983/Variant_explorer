# MODULE FOR DATA PLOTTING IN SHINY APPS
# Provide functionality for both standard plots and interactive plotly 
# Plotted variables can be manually configured or user-selected by multi-choice lists
 
# FUNCTIONS
# plotFixedUI + plotFixedModule: plot using manually defined variables
# plotSelectedUI + plotSelectedModule: plot using user-defined variables
# plotly argument (TRUE/FALSE) to the UI function set the type of output 

# FUNCTIONS INPUTS
# plot_data = data.frame
# plotType = string for the plot type (see below)
# plotOptions = named list of options to be passed to configure the plots (accepted size, position_dodge, scale, width)
# additionalOptions = list containing any ggplot configuration element (like scale, theme or facet)
# variables = named list of variables to map in aes list("x"="colnamex", "y"="colnamey")

# TO MADE USER CONFIGURABLE PLOTS
# provide variables argument to plotSelectedUI describing the variable than can be configured and the available options
# example: variables = list("x" = list("var1_label"="var1_name"), "y"=list("var2_label"="var2_name"))
# when calling plotSelectedModule in the server you can provide and additional set of variables as fixed mapping beside the user selected ones
# example: variables = list("color"="color_variable", "fill"="fill_variable")


selectPlot <- function(plotType, plotOptions) {
  switch(plotType,
         scatter = {p <- geom_point(size = ifelse(is.null(plotOptions[["size"]]),1,plotOptions[["size"]]))},
         barplot = {p <- geom_bar(aes(y=(..count..)), position = position_dodge(ifelse(is.null(plotOptions[["position_dodge"]]),0,plotOptions[["position_dodge"]])))},
         density = {p <- geom_density()},
         violin = {p <- geom_violin(scale = ifelse(is.null(plotOptions[["scale"]]),"width",plotOptions[["scale"]]))},
         boxplot = {p <- geom_boxplot(width = ifelse(is.null(plotOptions[["width"]]),1,plotOptions[["width"]]))},
         jitter = {p <- geom_jitter(width = ifelse(is.null(plotOptions[["width"]]),1,plotOptions[["width"]]))}
  )
  return(p) 
}

setValue <- function(id, input, variables) {
  if (id %in% names(input)) {
    return(input[[id]])
  } else if (id %in% names(variables)) {
    return(variables[[id]])
  } else {
    return(NULL)
  }  
}

##########################################
### Plot with fixed variables for axes ###
##########################################

plotFixedUI <- function(id,plotly=FALSE) {
  ns <- NS(id)
  RV$plotly = plotly
  
  if (plotly == TRUE) {
    tagList(plotlyOutput(ns("plotly")))
  } else {
    tagList(plotOutput(ns("plot")))
  }
} 

plotFixedModule <- function(input, output, session, plot_data, plotType, variables, plotOptions=NULL, additionalOptions=NULL) {
  p <- selectPlot(plotType, plotOptions)
  
  output$plotly <- renderPlotly({
    myplot <- ggplot(plot_data, aes_string(x = variables[["x"]], y = variables[["y"]], color= variables[["color"]], size=variables[["size"]], fill=variables[["fill"]], shape=variables[["shape"]])) + p     
    for (opt in additionalOptions) {
      myplot <- myplot + opt
    }
    ggplotly( myplot )
  })    
  
  output$plot <- renderPlot({
        myplot <- ggplot(plot_data, aes_string(x = variables[["x"]], y = variables[["y"]], color= variables[["color"]], size=variables[["size"]], fill=variables[["fill"]], shape=variables[["shape"]])) + p
        for (opt in additionalOptions) {
          myplot <- myplot + opt
        }
        myplot
  })
}

####################################################
### Plot with user-selectable variables for axes ###
####################################################

plotSelectedUI <- function(id,variables,plotly=FALSE) {
  ns <- NS(id)

  col_size <- 12/length(variables)
  commands_output <- NULL
  
  if (plotly == TRUE) {
    plot_output <- tagList(plotlyOutput(ns("plotly")))
  } else {
    plot_output <- tagList(plotOutput(ns("plot")))
  }
  
  for (v in names(variables)) {
    commands_output <- tagList(commands_output, 
                               column(col_size,selectInput(ns(v), paste0(v, " variable"), choices = variables[[v]], multiple = FALSE)) )
  }
  commands_output <- tagList(fluidRow(commands_output))
  output <- tagList(commands_output, plot_output)
}

plotSelectedModule <- function(input, output, session, plot_data, plotType, variables=NULL, plotOptions=NULL, additionalOptions=NULL) {
  p <- selectPlot(plotType, plotOptions)

  output$plotly <- renderPlotly({
    myplot <- ggplot(plot_data, 
                     aes_string(
                       x = setValue("x", input, variables),
                       y = setValue("y", input, variables),
                       color = setValue("color", input, variables),
                       size = setValue("size", input, variables),
                       fill = setValue("fill", input, variables),
                       shape = setValue("shape", input, variables)
                     )) + p
    for (opt in additionalOptions) {
      myplot <- myplot + opt
    }
    ggplotly ( myplot )
  })

  output$plot <- renderPlot({
    myplot <- ggplot(plot_data, 
                     aes_string(
                       x = setValue("x", input, variables),
                       y = setValue("y", input, variables),
                       color = setValue("color", input, variables),
                       size = setValue("size", input, variables),
                       fill = setValue("fill", input, variables),
                       shape = setValue("shape", input, variables)
                     )) + p
    for (opt in additionalOptions) {
      myplot <- myplot + opt
    }
    myplot
  }) 
}
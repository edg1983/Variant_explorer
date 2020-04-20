# MODULE FOR DATA PLOTTING IN SHINY APPS
# Provide functionality for both standard plots and interactive plotly 
# Plotted variables can be manually configured or user-selected by multi-choice lists
 
# FUNCTIONS
# plotFixedUI: plot using manually defined variables
# plotSelectedUI: plot using user-defined variables
# plotly argument (TRUE/FALSE) to the UI function set the type of output to static ggplot or plotly
# plotModule, function to be called on server

# FUNCTIONS INPUTS
# plot_data = data.frame
# plotType = string for the plot type (see below)
# plotOptions = named list of options to be passed to configure the plots (accepted size, position_dodge, scale, width)
# additionalOptions = list containing any ggplot configuration element (like scale, theme or facet)
# variables = named list of variables to map in aes list("x"="colnamex", "y"="colnamey")
# missingValues = vector of values representing missing data, all points with these values will be removed

# TO MADE USER CONFIGURABLE PLOTS
# provide variables argument to plotSelectedUI describing the variable than can be configured and the available options
# example: variables = list("x" = list("var1_label"="var1_name"), "y"=list("var2_label"="var2_name"))
# when calling plotModule in the server you can provide an additional set of variables as fixed mapping beside the user selected ones
# example: variables = list("color"="color_variable", "fill"="fill_variable")

# BIGDATA PLOTS
# A special plotType = "bigdata" is implemented to quickly plot scatterplots with a lot of points.
# This require scattermore library and is not compatible with plotly
# https://github.com/exaexa/scattermore

colors_palette <- c("red", "darkgreen", "blue", "purple", "magenta", "green", "orange", "brown")

`%nin%` = Negate(`%in%`)

selectPlot <- function(plotType, plotOptions) {
  switch(plotType,
         scatter = {p <- geom_point(size = ifelse(is.null(plotOptions[["size"]]),1,plotOptions[["size"]]))},
         barplot = {if (!is.null(plotOptions[["position_dodge"]])) {
                      p <- geom_bar(aes(y=(..count..)), position = position_dodge(plotOptions[["position_dodge"]]))  
                    } else {
                      p <- geom_bar(aes(y=(..count..)) )
                    } },
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

########################################
### UI with fixed variables for axes ###
########################################

plotFixedUI <- function(id,plotly=FALSE) {
  ns <- NS(id)
  RV$plotly = plotly
  
  if (plotly == TRUE) {
    tagList(plotlyOutput(ns("plotly")))
  } else {
    tagList(plotOutput(ns("plot")))
  }
} 

#################################################
### UI with user-selectable variables for axes ###
#################################################

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

################################
### plotting server function ###
################################

plotModule <- function(input, output, session, plot_data, plotType, missingValues=NULL, variables=NULL, plotOptions=NULL, additionalOptions=NULL) {
  
  myplot <- reactive ({
    x = setValue("x", input, variables)
    y = setValue("y", input, variables)
    color = setValue("color", input, variables)
    size = setValue("size", input, variables)
    fill = setValue("fill", input, variables)
    shape = setValue("shape", input, variables)
    label = setValue("label", input, variables)
    
    #Remove data points were either x or y has missing values
    if (!is.null(missingValues)) {
      filtered_df <- plot_data[plot_data[[x]] %nin% missingValues & plot_data[[y]] %nin% missingValues,]
    } else {
      filtered_df <- plot_data
    }
    
    #This is a special scatterplot approach for large datasets
    if (plotType == "bigdata") {
      p <- ggplot()
      filtered_df <- filtered_df[,c(x,y,color)]
      
      #check if color variable has been required and how many levels are in color columns
      if (!is.null(color)) {
        colors_levels <- unique(filtered_df[[color]])
        colors_values <- colors_palette[1:length(colors_levels)]
        for (i in 1:length(colors_levels)) {
          dfplot <- filtered_df[filtered_df[[color]] == colors_levels[i],c(x,y)]
          p <- p + geom_scattermost(dfplot,pointsize = size, color=colors_values[i])
        }
      }
      out_plot <- p + labs(x=x, y=y, subtitle = paste(colors_levels, colors_values, sep=": ", collapse="\t"))
    
    } else {
      p <- selectPlot(plotType, plotOptions)
      out_plot <- ggplot(plot_data, 
                       aes_string(
                         x = x,
                         y = y,
                         color = color,
                         size = size,
                         fill = fill,
                         shape = shape,
                         label = label
                       )) + p
    }
  })
  
  output$plotly <- renderPlotly({
    out_plot <- myplot()
    for (opt in additionalOptions) {
      out_plot <- out_plot + opt
    }
    ggplotly ( out_plot )
  })

  output$plot <- renderPlot({
    out_plot <- myplot()
    for (opt in additionalOptions) {
      out_plot <- out_plot + opt
    }
    out_plot
  }) 
}
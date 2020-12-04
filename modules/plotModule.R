# MODULE FOR DATA PLOTTING IN SHINY APPS
# Provide functionality for both standard plots and interactive plotly 
# Plotted variables can be manually configured or user-selected by multi-choice lists
 
# UI FUNCTIONS
# plotFixedUI: MADE PLOTS USING FIXED VARIABLES
# plotly = TRUE/FALSE; switch between ggplot or plotly output

# plotSelectedUI: MADE PLOTS WITH USER-CONFIGURABLE VARIABLES
# variables: describe the variables than can be configured and the available options
# variables = list("x" = list("var1_label"="var1_name"), "y"=list("var2_label"="var2_name"))
# when calling plotModule in the server you can provide an additional set of variables as fixed mapping beside the user selected ones
# example: variables = list("color"="color_variable", "fill"="fill_variable")
# plotly = TRUE/FALSE; switch between ggplot or plotly output
# set_limits = c("x","y"): vector describing if x and/or y axes may have configurable limits (this creates input text boxes)

# SERVER FUNCTIONS
# plotModule: GENERATE THE PLOTS
# plot_data = data.frame with data to be plotted
# plotType = string for the plot type (scatter, barplot, density, violin, boxplot, jitter, bigdata)
# plotOptions = named list of options to be passed to configure the plots aspect (accepted size, position_dodge, scale, width)
# additionalOptions = list containing any ggplot configuration element (like scale, theme or facet)
# variables = named list of variables to map in aes list("x"="colnamex", "y"="colnamey"). This acts in addition to the ones from input controls
# missingValues = vector of values representing missing data, all points with these values will be removed

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

resetRV <- function(values) {
  for (v in values) {
    RV[[v]] <- FALSE
  }
}

RV <- reactiveValues(
  x_min = FALSE,
  x_max = FALSE,
  y_min = FALSE,
  y_max = FALSE,
  controls = NULL
  )

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

##################################################
### UI with user-selectable variables for axes ###
##################################################

plotSelectedUI <- function(id,variables,set_limits=NULL,plotly=FALSE) {
  ns <- NS(id)

  col_size <- 12/length(variables)
  commands_output <- NULL
  sliders_output <- NULL
  
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
  
  if (!is.null(set_limits)) {
    for (axes in set_limits) {
      sliders_output <- tagList(sliders_output,
        column(4, 
          textInput(ns(paste0(axes,"_min")),label = paste0(axes," min:"), placeholder = paste0(axes," min")),
          textInput(ns(paste0(axes,"_max")),label = paste0(axes," max"), placeholder = paste0(axes," max")) 
        )
      )
    }
    sliders_output <- tagList(sliders_output, 
                              column(4, actionButton(ns("Apply_scale"),label = "Scale")),
                              column(4, actionButton(ns("Reset_scale"),label = "Reset")))
    sliders_output <- tagList(fluidRow(sliders_output))
  }
  
  output <- tagList(commands_output, sliders_output, plot_output)
}

################################
### plotting server function ###
################################

plotModule <- function(input, output, session, plot_data, plotType, missingValues=NULL, variables=NULL, plotOptions=NULL, additionalOptions=NULL) {
  
  toListen <- reactive({
    output <- list()
    i <- 0
    for (n in names(input)) {
      if (n %nin% c("x_min","x_max","y_min","y_max")) { 
        i <- i + 1
        output[[i]] <- input[[n]] 
      }
    }
    return(output)
  })
  
  #observeEvent(toListen(), {
  #  resetRV(c("x_min","x_max","y_min","y_max"))
  #})
  
  observeEvent(input$Reset_scale, {
    resetRV(c("x_min","x_max","y_min","y_max"))
  })
  
  observeEvent(input$Apply_scale, {
    RV$x_min <- setValue("x_min", input, 0)
    RV$x_max <- setValue("x_max", input, 0)
    RV$y_min <- setValue("y_min", input, 0)
    RV$y_max <- setValue("y_max", input, 0)
  })
  
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
      
      tryCatch({
        if (RV$x_min != "" & RV$x_max != "" & RV$x_min != FALSE & RV$x_max != FALSE) {
          filtered_df <- filtered_df[filtered_df[[x]] >= RV$x_min & filtered_df[[x]] <= RV$x_max,] 
        }
        
        if (RV$y_min != "" & RV$y_max != "" & RV$y_min != FALSE & RV$y_max != FALSE) {
          filtered_df <- filtered_df[filtered_df[[y]] >= RV$y_min & filtered_df[[y]] <= RV$y_max,] 
        }
      }, error=function(cond) {}
      )
      
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
      out_plot <- ggplot(filtered_df, 
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
#################
### FUNCTIONS ###
#################
`%nin%` = Negate(`%in%`)

limitValue <- function(values, operation, na_value) {
  if (!is.null(na_value)) {
    switch(operation,
           min = { n <- round(min(values[values != na_value]),3) },
           max = { n <- round(max(values[values != na_value]),3) }
    )
  } else {
    switch(operation,
           min = { n <- round(min(values, na.rm = T), 3)},
           max = { n <- round(max(values, na.rm = T), 3)}
    )
  }
  return(n)
}

makeControl <- function(ctrl_type, ctrl_id, ctrl_label, ctrl_value, factorsmap, var_name, df, na_values, step=0.001, ns) {
  switch(ctrl_type,
         factors_fields = {
           if (var_name %in% names(factorsmap)) {
             choices <- factorsmap[[var_name]]
           } else {
             choices <- as.list(sort(unique(df[[var_name]])))
             names(choices) <- sort(unique(df[[var_name]])) 
           }
           input_element <- tagList(selectInput(ns(ctrl_id), ctrl_label, choices = choices, selected = ctrl_value, multiple = TRUE))
         },
         numeric_fields = {
           min <- limitValue(df[[var_name]], "min", na_values[[var_name]])
           max <- limitValue(df[[var_name]], "max", na_values[[var_name]])
           if (ctrl_value %in% c("min", "max")) {
             default <- limitValue(df[[var_name]], ctrl_value, na_values[[var_name]])
           } else {
             default <- as.numeric(ctrl_value)
           }
           input_element <-tagList(
             fluidRow(column(9, sliderInput(ns(ctrl_id), ctrl_label, min = min, max = max, value = default, step = step)),
                      column(3, textInput(ns(paste0("TXTSET_",ctrl_id)),"",placeholder = "value")) )
           )
         },
         binary_fields = {
           default <- as.logical(ctrl_value)
           input_element <-tagList(checkboxInput(ns(ctrl_id), ctrl_label, value = default))
         }
  )
  return(input_element)
}

groupDescription <- function(group_name, definition) {
  if (group_name == "global") {
    out_text <- "These filters are applied to all records"
  } else if (group_name == "comphet") {
    out_text <- "At least one of the comphet variants must pass these filters"
  } else {
    out_text <- paste0("Applied when ", definition[1], " included in: ", paste(definition[2][[1]], collapse=","))
  }
  return(out_text)
}

ctrlLayout <- function(group_controls) {
  ncols <- sum(unlist(lapply(group_controls, function(x) !is.null(x))))
  width <- 12/ncols
  controls_layout <- NULL
  for (controls in group_controls) {
    if (!is.null(controls)) {controls_layout <- tagList(controls_layout, column(width, controls))}
  }
  controls_layout <- tagList(fluidRow(controls_layout))
  return(controls_layout)
}

getObjects <- function(input_list, pattern) {
  controls <- names(input_list)[grep(pattern, names(input_list))]
  observed_input <- NULL
  for (n in controls) {
    observed_input[[n]] <- input_list[[n]]
  } 
  return(observed_input)  
}

updateValues <- function(session, preset_value, preset_config, group_vars, vars_type, vars_default) {
  #Update the logic control
  logicvalue <- preset_config$configuration[[preset_value]]$logic
  updateSelectInput(session, inputId = "LOGIC", selected = logicvalue)
  
  #Check which values are required by the presets and update them.
  #Other group values are set to default
  newvalues <- preset_config$configuration[[preset_value]]$values
  for (v in names(group_vars)) {
    if (v %in% names(newvalues)) { 
      value = newvalues[[v]] 
    } else {
      value = vars_default[[v]]
    }
    
    switch(vars_type[[v]],
           numeric_fields = {updateSliderInput(session, inputId = v, value = as.numeric(value))},
           factors_fields = {updateSelectInput(session, inputId = v, selected = value)},
           binary_fields = {updateCheckboxInput(session, inputId = v, value = as.logical(value))}
    )
  }
}

makeGroupExpression <- function(logic, filters_operations, filters_values, vars_type) {
  filters_expr <- NULL
  
  #For each variant type the expression is built with the desired logic
  #The resulting variant type expressions are then merged with AND
  for (var_type in names(vars_type)) {
    vars_in_type <- intersect(names(filters_values), names(vars_type[[var_type]]))
    type_expr <- NULL
    for (var_name in vars_in_type) {
      #Make the expression
      myexpr <- makeExpression(var_name, filters_operations[[var_name]], filters_values[[var_name]])
      
      #Check is the expression returned NULL
      #This manage unchecked boxes for which we don't want to add any expression
      if (!is.null(myexpr)) {
        
        #If the var name contains the special tag _SNP / _INDEL then expr is valid only on SNP/INDEL
        if (grepl("_SNP", var_name)) {
          myexpr <- expr((!!as.name("var_type") == "SNV" & !!myexpr))
        } else if(grepl("_INDEL", var_name)) {
          myexpr <- expr((!!as.name("var_type") == "INDEL" & !!myexpr))  
        }

        #Merge expressions using the provided logic
        if (logic == "AND") {
          if (is.null(type_expr)) {
            type_expr <- myexpr
          } else {
            type_expr <- expr(!!type_expr & !!myexpr)
          }      
        } else {
          if (is.null(type_expr)) {
            type_expr <- myexpr
          } else {
            type_expr <- expr(!!type_expr | !!myexpr)
          }  
        }
      }
    }
    #Evaluate if the expression for the variant type is not NULL
    #This avoid generate NULL expression when all checkboxes are false
    if(!is.null(type_expr)){
      if (is.null(filters_expr)) {
        filters_expr <- type_expr
      } else {
        filters_expr <- expr(!!filters_expr & !!type_expr)
      }
    }
  }
  return(filters_expr)
}

makeExpression <- function(var,operator,value) {
  #Set the quo expression for the group. 
  switch(operator,
         "grep" = { 
           value <- paste(value, collapse="|")
           myexpr <- expr(grepl(!!value,!!(as.name(var))))
         },
         "include" = { myexpr <- evaluateCheckBox(var,operator,value) },
         "exclude" = { myexpr <- evaluateCheckBox(var,operator,value) },
         "%nin%" = { myexpr <- expr(!!as.name(var) %nin% !!value) },
         "%in%" = { myexpr <- expr(!!as.name(var) %in% !!value) },
         "==" = { myexpr <- expr(!!as.name(var) == !!as.numeric(value)) },
         ">" = { myexpr <- expr(!!as.name(var) > !!as.numeric(value)) },
         "<" = { myexpr <- expr(!!as.name(var) < !!as.numeric(value)) },
         ">=" = { myexpr <- expr(!!as.name(var) >= !!as.numeric(value)) },
         "<=" = { myexpr <- expr(!!as.name(var) <= !!as.numeric(value)) }
  )
  
  return(myexpr)
}

evaluateCheckBox <- function(var, operator="include", value) {
  #Checkboxes are used to include or exclude specific regions
  #Based on desired operation true value can be 0 or 1 if box is checked (TRUE)
  #IF FALSE both zero and one values are included
  false_value <- c(0,1)
  if (operator=="include") { 
    true_value <- 1
  }
  if (operator=="exclude") { 
    true_value <- 0
  }
  
  if (value == TRUE) {
    myexpr <- expr(!!(as.name(var)) == !!(as.numeric(true_value)))
  } else {
    myexpr  = NULL
  }
  return(myexpr)
}

getDefaultValue <- function(var_name, ctrl_type, ctrl_value, na_values, df) {
  switch(ctrl_type,
         numeric_fields = {
           if (ctrl_value %in% c("min", "max")) {
             default <- limitValue(df[[var_name]], ctrl_value, na_values[[var_name]])
           } else {
             default <- as.numeric(ctrl_value)
           }
         },
         factors_fields = {
           if (ctrl_value == "ALL") {
             default <- unique(df[[var_name]])
           } else if (base::grepl("SPLIT", ctrl_value)) {
             sep = gsub("SPLIT ", "",ctrl_value)
             all_values = paste(unique(df[[var_name]]), collapse=sep)
             default = unique(unlist(strsplit(all_values,split = sep)))
           }else {
             default <- ctrl_value
           }
         },
         binary_fields = {
           default <- ctrl_value
         }
  )
  return(default)
}

#############################
### USER INTERFACE GROUPS ###
#############################

# filters_settings = VARIANTS section of filters_settings json
# variants_df = data frame of variants 
# na_values = named list with na values (usually app_settings$fill_na$fill_na_vars)
# This will generate one shinydashboard box for each defined variant group
# Within each box, 3 columns are arranged for factor, numeric and binary filters respectively
# Two additional multichoice ar added at the top for each group 
#   - AND / OR logic
#   - Preset if present in settings json

filterObjectUI <- function(id, group_name, group_definition, group_presets, vars_type, vars_default, vars_label, variants_df, na_values, factorsmap, tooltips, filters_settings) {
  ns <- NS(id)
  
  vars_type <- vars_type
  vars_default <- vars_default
  na_values <- na_values
  group_presets <- group_presets
  variants_df <- variants_df
  
  #Create all controls and store them
  UI_elements <- list()
  for (var_name in names(vars_type)) {
        UI_elements[[var_name]] <- makeControl(
          ctrl_type = vars_type[[var_name]], 
          ctrl_id = var_name, 
          ctrl_label = vars_label[[var_name]],
          ctrl_value = vars_default[[var_name]],
          factorsmap = factorsmap,
          var_name = var_name,
          df = variants_df,
          na_values = na_values,
          ns=ns)
  }

  #Build shiny dashboard box
  group_title <- paste0(group_name, " filters")
  group_description <- groupDescription(group_name, group_definition)
  group_presets <- tagList(selectInput(ns("PRESETS"), "Presets:", choices = names(group_presets$configuration), selected = group_presets$default, multiple = FALSE))
  group_logic <- selectInput(inputId = ns("LOGIC"), label = "Filters logic:", choices = c("AND", "OR"), selected = "OR", multiple = FALSE)
  
  group_controls <- list(NULL, NULL, NULL)
  for (var_name in names(vars_type)) {
    switch (vars_type[[var_name]],
            numeric_fields = {group_controls[[1]] <- tagList(group_controls[[1]], UI_elements[[var_name]], bsTooltip(ns(var_name), tooltips[[var_name]]))},
            factors_fields = {group_controls[[2]] <- tagList(group_controls[[2]], UI_elements[[var_name]], bsTooltip(ns(var_name), tooltips[[var_name]]))},
            binary_fields = {group_controls[[3]] <- tagList(group_controls[[3]], UI_elements[[var_name]], bsTooltip(ns(var_name), tooltips[[var_name]]))}
    )
  }
  controls_layout <- ctrlLayout(group_controls)
  group_box <- tagList(
    box(title = group_title, width = 12, status = "primary", solidHeader = T, collapsible = T,
        fluidRow(column(6,group_description), column(3, group_logic), column(3, group_presets)),
        hr(),
        controls_layout
    )
  )
  callModule(observeInputs, id, filters_settings=filters_settings, na_values=na_values, variants_df=variants_df)
  return(group_box)
}

##############
### SERVER ###
##############

#Manage the update of filters values when preset or text input change
observeInputs <- function(input, output, session, filters_settings, na_values, variants_df) {
  
  #Create reactive object containing the numeric text inputs
  txtControls <- reactive({
    getObjects(input, "TXTSET")
  })
  
  #Observe text inputs and update sliders values accordingly
  observeEvent(txtControls(), {
    for (n in names(txtControls())) {
      if (txtControls()[[n]] != "") {
        sliderid <- gsub("TXTSET_","",n)
        updateSliderInput(session, inputId = sliderid, value = txtControls()[[n]])
      }
    }  
  })
  
  observeEvent(input$PRESETS, {
    x <- session$ns('tmp')  # make an ID string
    group_name <- gsub("variants_filters-|-tmp", "", x)
    group_name <- gsub("genes_filters-|-tmp", "", group_name)

    vars_type <- list()
    vars_default <- list()
    for (t in names(filters_settings$DEFINITIONS)) {
      for (v in names(filters_settings$DEFINITIONS[[t]])) {
        vars_type[[v]] <- t
        vars_default[[v]] <- getDefaultValue(var_name = v,
                                             ctrl_type = t,
                                             ctrl_value = filters_settings$DEFINITIONS[[t]][[v]][[3]],
                                             na_values = na_values,
                                             df = variants_df)
      }
    }
    
    group_presets <- filters_settings$PRESETS[[group_name]]
    
    vars_type[["LOGIC"]] <- "factors_fields"
    group_vars <- names(input)
    group_vars <- group_vars[-grep("TXTSET|PRESETS",group_vars)]
    
    #Update the logic control
    logicvalue <- group_presets$configuration[[input$PRESETS]]$logic
    updateSelectInput(session, inputId = "LOGIC", selected = logicvalue)
    
    #Check which values are required by the presets and update them.
    #Other group values are set to default
    newvalues <- group_presets$configuration[[input$PRESETS]]$values
    for (v in group_vars) {
      if (v %in% names(newvalues)) { 
        value = newvalues[[v]] 
      } else {
        value = vars_default[[v]]
      }
      switch(vars_type[[v]],
             numeric_fields = {updateSliderInput(session, inputId = v, value = as.numeric(value))},
             factors_fields = {updateSelectInput(session, inputId = v, selected = value)},
             binary_fields = {updateCheckboxInput(session, inputId = v, value = as.logical(value))}
      )
    }
  })
}

#Make quo expression combining all filters for the group
getFilterExpression <- function(input, output, session, group_name, group_definition, group_vars, vars_operation, vars_definition) {
  ns <- session$ns
  
  #Build a group_def variable with the group definition
  group_def <- list()
  group_def$field <- group_definition[[1]]
  group_def$values <- group_definition[[2]]
  
  #make a named list associated to the group: NAME=var_id, VALUE=value from control
  filters_values <- list() 
  for (v in group_vars) {
    filters_values[[v]] <- input[[v]]
  }
  
  if (group_name == "global") {
    #global is a special group indicating filters that will be applied to all vars
    group_expr<- makeGroupExpression(logic = input$LOGIC,
                                     filters_operations = vars_operation,
                                     filters_values = filters_values,
                                     vars_type = vars_definition)
    pass_expr <- group_expr
    filter_expr <- expr(!(!!group_expr))
  } else if (group_name == "comphet") {
    #comphet is a special group that will generate a separate expression returned for comphet
    group_expr<- makeGroupExpression(logic = input$LOGIC,
                                     filters_operations = vars_operation,
                                     filters_values = filters_values,
                                     vars_type = vars_definition)
    pass_expr <- group_expr
    filter_expr <- expr(!(!!group_expr))
  } else {
    #for group filters a group_expr and a filter_expr are combined to evaluate vars only in the selected group
    groupid_expr <- expr(!!as.name(group_def$field) %in% !!group_def$values)
    groupfilters_expr <- makeGroupExpression(logic = input$LOGIC,
                                             filters_operations = vars_operation,
                                             filters_values = filters_values,
                                             vars_type = vars_definition)
    
    pass_expr <- expr(!!groupid_expr & (!!groupfilters_expr))
    filter_expr <- expr(!!groupid_expr & !(!!groupfilters_expr))
  }
  
  return(list(pass=pass_expr, filter=filter_expr))
}

#Return a data frame with the filters values
getDF <- function(input, output, session, group_name) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)
  remove <- grep("TXTSET|PRESETS", ctrl_names)
  ctrl_names <- ctrl_names[-remove]
    for (n in ctrl_names) {
      newline <- c(group_name, n, paste(input[[n]], collapse=","))
      filters_df[nrow(filters_df)+1,] <- newline
    }
  filters_df <- filters_df[order(filters_df$filter),]
  return(filters_df)
}

#Return a list that can be used to make a JSON of the filters settings
getJSON <- function(input, output, session) {
  filters_json <- list()
  ctrl_names <- names(input)
  remove <- grep("TXTSET", ctrl_names)
  ctrl_names <- ctrl_names[-remove]
    for (n in ctrl_names) {
      filters_json[[n]] <- input[[n]]
    }
  return(filters_json) 
}

#Read filters settings from a list of saved values
loadSettings <- function(input, output, session, filters_settings, filters_json) {
  x <- session$ns('tmp')  # make an ID string
  group_name <- gsub("variants_filters-|-tmp", "", x)
  group_name <- gsub("genes_filters-|-tmp", "", group_name)
  
  vars_type <- list()
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      vars_type[[v]] <- t
    }
  }
  vars_type[["LOGIC"]] <- "factors_fields"
  vars_type[["PRESET"]] <- "factors_fields"
  
  filters_values <- filters_json[[group_name]] 
  
    for (v in names(filters_values)) {
      value <- filters_values[[v]]
      if (v %in% names(vars_type)) {
        switch(vars_type[[v]],
               numeric_fields = {updateSliderInput(session, inputId = v, value = as.numeric(value))},
               factors_fields = {updateSelectInput(session, inputId = v, selected = value)},
               binary_fields = {updateCheckboxInput(session, inputId = v, value = as.logical(value))}
        )
      } else {
        message("LOAD WARNING: ",v, " variable not found in the DEFINITIONS from your filters settings")
      }
    }

}
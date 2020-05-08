# MODULE FOR VARIANT FILTERS MANAGEMENT
# Allow to configure groups of vars and apply specific filters to each group
# Suppose a variants table and a gene scores table containing various annotations columns and a json config file
# Filters are configured as follows:
# global_filters | (var in group1 & (group1 filters)) | (var in group1 & (group1 filters)) ...
# Returns:
# - list of vars ids passing the filters
# - data frame of the applied filters in the form of control_id, value

# Filter configuration is defined in a json file
# DEFINITIONS:  contains 3 groups: numeric_fields, factors_fields and binary_fields
#               each groups define column names and the associated label to show.
#               for numerical columns the operator of filter is also provided like >,<
#               for factors column %in% is used
#               an additional group reg_sources can be provided. This must contain field indicating the column
#               then contains subgroups defining the various reg_db categories
# GENES / VARIANTS:     allow to define filters groups and presets

#################
### FUNCTIONS ###
#################
`%nin%` = Negate(`%in%`)

#Set the quo expression
selectFilter <- function(var,operator,value) {
  switch(operator,
         "%nin%" = { myexpr <- quo((!!(as.name(var))) %in% !!(as.numeric(value))) },
         "%in%" = { myexpr <- quo((!!(as.name(var))) %nin% !!(as.numeric(value))) },
         "==" = { myexpr <- quo((!!(as.name(var))) == !!(as.numeric(value))) },
         ">" = { myexpr <- quo((!!(as.name(var))) > !!(as.numeric(value))) },
         "<" = { myexpr <- quo((!!(as.name(var))) < !!(as.numeric(value))) },
         ">=" = { myexpr <- quo((!!(as.name(var))) >= !!(as.numeric(value))) },
         "<=" = { myexpr <- quo((!!(as.name(var))) <= !!(as.numeric(value))) }
  )
  
  return(myexpr)
}

makeBoolean <- function(myvalues,myvector, logic="include") {
  if (logic=="include") { 
    true_value <- 1
  }
  if (logic=="exclude") { 
    true_value <- 0
  }
  
  outlist <- list()
  for (v in myvalues) {
    if (v %in% myvector) {
      outlist[[v]] = true_value
    } else {
      outlist[[v]] = c(0,1)
    }
  }
  return(outlist)
}

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

makeControl <- function(ctrl_type, ctrl_id, ctrl_value, var_name, df, na_values, step=0.001, ns) {
  switch(ctrl_type,
         factors_fields = {
           choices <- sort(unique(df[[var_name]]))
           names(choices) <- choices
           choices <- c("ALL" = "ALL", choices)
           ctrl_label <- ctrl_value[1]
           default <- ctrl_value[3]
           input_element <- tagList(selectInput(ns(ctrl_id), ctrl_label, choices = choices, selected = default, multiple = TRUE))
         },
         numeric_fields = {
           ctrl_label <- ctrl_value[[1]]
           min <- limitValue(df[[var_name]], "min", na_values[[var_name]])
           max <- limitValue(df[[var_name]], "max", na_values[[var_name]])
           if (ctrl_value[3] %in% c("min", "max")) {
            default <- limitValue(df[[var_name]], ctrl_value[[3]], na_values[[var_name]])
           } else {
            default <- as.numeric(ctrl_value[[3]])
           }
           input_element <-tagList(
             fluidRow(column(9, sliderInput(ns(ctrl_id), ctrl_label, min = min, max = max, value = default, step = step)),
                      column(3, textInput(ns(paste0("TXTSET_",ctrl_id)),"",placeholder = "value")) )
           )
         },
         binary_fields = {
           ctrl_label <- paste0(ctrl_value[2], " ", ctrl_value[1])
           default <- as.logical(ctrl_value[3])
           input_element <-tagList(checkboxInput(ns(ctrl_id), ctrl_label, value = default))
         }
  )
  return(input_element)
}

groupDescription <- function(group_name, definition) {
  if (group_name == "global") {
    out_text <- "These filters are applied to all variants"
  } else {
    out_text <- paste0("Applied to variants when ", definition[1], " included in: ", paste(definition[2][[1]], collapse=","))
  }
  return(out_text)
}

presetsControl <- function(id, presets_config, ns) {
  presets <- tagList(selectInput(ns(id), "Presets:", choices = names(presets_config$configuration), selected = presets_config$default, multiple = FALSE))
  return(presets)
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

updateValues <- function(session, preset_id, preset_value, preset_config, var_type, var_default, groups, pattern="PRESETS_") {
  groupid <- gsub(pattern = pattern, "", preset_id)
  group_vars <- groups[[groupid]]$associated_values
  message("### GROUP ID: ",groupid, " ###")
  logicvalue <- preset_config[[groupid]]$configuration[[preset_value]]$logic
  controlid <- paste0("LOGIC_", groupid)
  updateSelectInput(session, inputId = controlid, selected = logicvalue)
  
  newvalues <- preset_config[[groupid]]$configuration[[preset_value]]$values
  message("values to update: ", paste(names(newvalues), collapse="; "))
  for (v in group_vars) {
    controlid <- paste0(v,"_",groupid)
    if (v %in% names(newvalues)) { 
      value = newvalues[[v]] 
    } else {
      value = var_default[[v]]
    }
    message("control id: ", controlid, "; vartype: ", var_type[[v]])
    message("\n\t newvalue: ", value)
    switch(var_type[[v]],
           numeric_fields = {updateSliderInput(session, inputId = controlid, value = as.numeric(value))},
           factors_fields = {updateSelectInput(session, inputId = controlid, selected = value)},
           binary_fields = {updateCheckboxInput(session, inputId = controlid, value = as.logical(value))}
          )
  }
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

filtersVariantsUI <- function(id, filters_settings, variants_df, na_values) {
  ns <- NS(id)
  
  #Build a group_def variable identifying belonging group for each control
  group_def <- list()
  for (n in names(filters_settings$GROUPS)) {
    for (f in filters_settings$GROUPS[[n]][["associated_values"]]) {group_def[[f]] <- c(group_def[[f]], n)}
  }
  
  #Build a var_type variable identifying each variable type
  var_type <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_type[[f]] <- n}
  }
  

  
  #Create all controls and store them
  UI_elements <- list()
  for (ctrl_type in names(filters_settings$DEFINITIONS)) {
    for (var_name in names(filters_settings$DEFINITIONS[[ctrl_type]])) {
      for (ctrl_group in group_def[[var_name]]) {
        ctrl_id <- paste0(var_name, "_", ctrl_group)
        UI_elements[[ctrl_id]] <- makeControl(
        ctrl_type = ctrl_type, 
        ctrl_id = ctrl_id, 
        ctrl_value = filters_settings$DEFINITIONS[[ctrl_type]][[var_name]],
        var_name = var_name,
        df = variants_df,
        na_values = na_values,
        ns=ns)
      }
    }
  }
  
  #Build shiny dashboard boxes for each group
  variants_boxes <- NULL
  for (filter_group in names(filters_settings$GROUPS)) {
    group_title <- paste0(filter_group, " filters")
    group_description <- groupDescription(filter_group, filters_settings$GROUPS[[filter_group]][["definition"]])
    group_presets <- presetsControl(id = paste0("PRESETS_",filter_group), presets_config = filters_settings$PRESETS[[filter_group]], ns=ns)
    group_logic <- selectInput(inputId = ns(paste0("LOGIC_",filter_group)), label = "Filters logic:", choices = c("AND", "OR"), selected = "OR", multiple = FALSE)
    
    group_controls <- list(NULL, NULL, NULL)
    for (var_name in filters_settings$GROUPS[[filter_group]][["associated_values"]]) {
      ctrl_id <- paste0(var_name, "_", filter_group)
      switch (var_type[[var_name]],
        numeric_fields = {group_controls[[1]] <- tagList(group_controls[[1]], UI_elements[[ctrl_id]])},
        factors_fields = {group_controls[[2]] <- tagList(group_controls[[2]], UI_elements[[ctrl_id]])},
        binary_fields = {group_controls[[3]] <- tagList(group_controls[[3]], UI_elements[[ctrl_id]])}
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
    variants_boxes <- tagList(variants_boxes, group_box)
  }
  variants_boxes
}

##############
### SERVER ###
##############

filtersVariants <- function(input, output, session, filters_settings, variants_df, na_values) {
  test_vars <- c("var1", "var2")
  
  #Make quo expression combining all filters
  filters_expr <- reactive({
    final_expr <- NULL
    for (g in names(filters_settings$GROUPS)) {
      filters_values <- list() #make a names list containing values from controls associated to the group
      
      group_expr <- makeGroupExpression(groupid = g, #if g == global need special logic
                                  group_def = filters_settings$GROUPS[[g]]$definition,
                                  filters_values = filters_values)
      final_expr <- quo(!!final_expr & !!group_expr)
    }
    return(final_expr)
  })
  
  #Build a group_def variable identifying belonging group for each control
  group_def <- list()
  for (n in names(filters_settings$GROUPS)) {
    for (f in filters_settings$GROUPS[[n]][["associated_values"]]) {group_def[[f]] <- c(group_def[[f]], n)}
  }
  
  #Build a var_type variable identifying each variable type
  var_type <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_type[[f]] <- n}
  }
  
  #Build a var_default variable identifying each variable default value
  var_default <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_default[[f]] <- filters_settings$DEFINITIONS[[n]][[f]][3]}
  }
  
  #Create reactive object containing the numeric text inputs
  txtControls <- reactive({
    getObjects(input, "TXTSET")
  })
  
  #Create reactive object containing the preset select inputs
  presetControls <- reactive({
    getObjects(input, "PRESETS")
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
  
  observeEvent(presetControls(), {
    for (n in names(presetControls())) {
      updateValues(session = session, 
                   preset_id = n, 
                   preset_value = presetControls()[[n]], 
                   preset_config = filters_settings$PRESETS, 
                   var_type = var_type, 
                   var_default = var_default, 
                   groups = filters_settings$GROUPS)
    } 
  })
  
  return(test_vars)
}  

getFiltersDF <- function(input, output, session, filters_settings) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)
  remove <- grep("TXTSET|PRESETS", ctrl_names)
  ctrl_names <- ctrl_names[-remove]
  for (g in names(filters_settings$GROUPS)) {
    group_string <- paste0("_",g)
    ctrls_group <- grep(group_string, ctrl_names)
    for (n in ctrl_names[ctrls_group]) {
      newline <- c(g, gsub(group_string,"",n), paste(input[[n]], collapse=","))
      filters_df[nrow(filters_df)+1,] <- newline
    }
  }
  filters_df <- filters_df[order(filters_df$group, filters_df$filter),]
  return(filters_df)
}


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
  logicvalue <- preset_config[[groupid]]$configuration[[preset_value]]$logic
  controlid <- paste0("LOGIC_", groupid)
  updateSelectInput(session, inputId = controlid, selected = logicvalue)
  
  newvalues <- preset_config[[groupid]]$configuration[[preset_value]]$values
  for (v in group_vars) {
    controlid <- paste0(v,"_",groupid)
    if (v %in% names(newvalues)) { 
      value = newvalues[[v]] 
    } else {
      value = var_default[[v]]
    }
  
    switch(var_type[[v]],
           numeric_fields = {updateSliderInput(session, inputId = controlid, value = as.numeric(value))},
           factors_fields = {updateSelectInput(session, inputId = controlid, selected = value)},
           binary_fields = {updateCheckboxInput(session, inputId = controlid, value = as.logical(value))}
          )
  }
}

makeGroupExpression <- function(logic, filters_operations, filters_values) {
  filters_expr <- NULL
  for (var_name in names(filters_values)) {
    expr <- makeExpression(var_name, filters_operations[[var_name]], filters_values[[var_name]])
    #Logic is reversed on purpose since the final expression will be negated
    if (logic == "AND") {
      if (is.null(filters_expr)) {
        filters_expr <- expr
      } else {
        filters_expr <- quo(!!filters_expr | !!expr)
      }      
    } else {
      if (is.null(filters_expr)) {
        filters_expr <- expr
      } else {
        filters_expr <- quo(!!filters_expr & !!expr)
      }  
    }
  }
  return(filters_expr)
}

makeExpression <- function(var,operator,value) {
  #Set the quo expression for the group. 
  #Operators are reverted since final filter is a negation
  switch(operator,
         "grep" = { 
           value <- paste(value, collapse="|")
           myexpr <- quo(!grepl(!!value,!!(as.name(var))))
           },
         "include" = { myexpr <- evaluateCheckBox(var,operator,value) },
         "exclude" = { myexpr <- evaluateCheckBox(var,operator,value) },
         "%nin%" = { myexpr <- quo(!!(as.name(var)) %in% !!value) },
         "%in%" = { myexpr <- quo(!!(as.name(var)) %nin% !!value) },
         "==" = { myexpr <- quo(!!(as.name(var)) != !!(as.numeric(value))) },
         ">" = { myexpr <- quo(!!(as.name(var)) < !!(as.numeric(value))) },
         "<" = { myexpr <- quo(!!(as.name(var)) > !!(as.numeric(value))) },
         ">=" = { myexpr <- quo(!!(as.name(var)) <= !!(as.numeric(value))) },
         "<=" = { myexpr <- quo(!!(as.name(var)) >= !!(as.numeric(value))) }
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
    expr <- quo(!!(as.name(var)) != !!(as.numeric(true_value)))
  } else {
    expr  = quo(!!(as.name(var)) %nin% !!false_value)
  }
  return(expr)
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
             default = unique(unlist(strsplit(all_values,split = ",")))
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

filtersVariantsUI <- function(id, filters_settings, variants_df, na_values) {
  ns <- NS(id)
  
  #Build a group_vars variable identifying belonging group for each control
  group_vars <- list()
  for (n in names(filters_settings$GROUPS)) {
    for (f in filters_settings$GROUPS[[n]][["associated_values"]]) {group_vars[[f]] <- c(group_vars[[f]], n)}
  }
  
  #Build a var_type variable identifying each variable type
  var_type <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_type[[f]] <- n}
  }
  
  #Build a var_default variable identifying each variable default value
  var_default <- list()
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      var_default[[v]] <- getDefaultValue(var_name = v,
                                          ctrl_type = t,
                                          ctrl_value = filters_settings$DEFINITIONS[[t]][[v]][[3]],
                                          na_values = na_values,
                                          df = variants_df)
    }
  }
  
  #Create all controls and store them
  UI_elements <- list()
  for (ctrl_type in names(filters_settings$DEFINITIONS)) {
    for (var_name in names(filters_settings$DEFINITIONS[[ctrl_type]])) {
      for (ctrl_group in group_vars[[var_name]]) {
        ctrl_id <- paste0(var_name, "_", ctrl_group)
        UI_elements[[ctrl_id]] <- makeControl(
        ctrl_type = ctrl_type, 
        ctrl_id = ctrl_id, 
        ctrl_label = filters_settings$DEFINITIONS[[ctrl_type]][[var_name]][1],
        ctrl_value = var_default[[var_name]],
        factorsmap = filters_settings$FACTORSMAP,
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

#Manage the update of filters values when preset or text input change
observeFilters <- function(input, output, session, filters_settings, variants_df, na_values) {
  #Build a group_vars variable identifying belonging group for each control
  group_vars <- list()
  for (n in names(filters_settings$GROUPS)) {
    for (f in filters_settings$GROUPS[[n]][["associated_values"]]) {group_vars[[f]] <- c(group_vars[[f]], n)}
  }
  
  #Build a var_type variable identifying each variable type
  var_type <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_type[[f]] <- n}
  }
  
  #Build a var_default variable identifying each variable default value
  var_default <- reactive({
  default_values <- list() 
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      #var_default[[v]] <- filters_settings$DEFINITIONS[[t]][[v]][[3]]
      default_values[[v]] <- getDefaultValue(var_name = v,
                                          ctrl_type = t,
                                          ctrl_value = filters_settings$DEFINITIONS[[t]][[v]][[3]],
                                          na_values = na_values,
                                          df = variants_df)
    }
  }
  return(default_values)
  })
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
                   var_default = var_default(), 
                   groups = filters_settings$GROUPS)
    } 
  })
}  

#Return a data frame with the filters values
getDF_filters <- function(input, output, session, filters_settings) {
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

#Make quo expression combining all filters for variants and return PASS variants ids
getPASSVars_filters <- function(input, output, session, filters_settings, variants_df, comphet_df) {
  #Build a var_operation variable identifying each variable default value
  var_operation <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_operation[[f]] <- filters_settings$DEFINITIONS[[n]][[f]][[2]]}
  }
  
  #Build a group_def variable with the group definition
  group_def <- list()
  for (n in names(filters_settings$GROUPS)) {
    group_def[[n]]$field <- filters_settings$GROUPS[[n]]$definition[[1]]
    group_def[[n]]$values <- unlist(filters_settings$GROUPS[[n]]$definition[[2]])
  }
  
  final_expr <- list()
  for (g in names(filters_settings$GROUPS)) {
    #make a named list associated to the group: NAME=var_id, VALUE=value from control
    filters_values <- list() 
    for (v in filters_settings$GROUPS[[g]]$associated_values) {
      filterid <- paste0(v,"_",g)
      filters_values[[v]] <- input[[filterid]]
    }
    
    if (g == "global") {
      #global is a special group indicating filters that will be applied to all vars
      group_expr<- makeGroupExpression(logic = input[[paste0("LOGIC_",g)]],
                                       filters_operations = var_operation,
                                       filters_values = filters_values)
    } else if (g == "comphet") {
      #comphet is a special group that will generate a separate expression returned for comphet
      group_expr<- makeGroupExpression(logic = input[[paste0("LOGIC_",g)]],
                                       filters_operations = var_operation,
                                       filters_values = filters_values)
      final_expr$comphet <- group_expr
      next
    } else {
      #for group filters a group_expr and a filter_expr are combined to evaluate vars only in the selected group
      groupid_expr <- quo((!!(as.name(group_def[[g]]$field))) %in% !!group_def[[g]]$values)
      groupfilters_expr <- makeGroupExpression(logic = input[[paste0("LOGIC_",g)]],
                                               filters_operations = var_operation,
                                               filters_values = filters_values)
      
      group_expr <- quo(!!groupid_expr & (!!groupfilters_expr))
    }
    if (is.null(final_expr$vars)) {
      final_expr$vars <- group_expr
    } else {
      final_expr$vars <- quo(!!final_expr$vars | (!!group_expr))
    }
    #message("Quo vars: ", final_expr$vars)
    #message("Quo comphet: ", final_expr$comphet)
  }
  #SINGLE VARIANTS
  #Variants passing all filters
  pass_vars <- as.data.frame(variants_df %>% filter(! (!!final_expr$vars)))
  
  #COMPHET VARS
  #1. Select comphet where both vars pass the variants filters
  pass_comphet <- comphet_df %>% filter(v1 %in% pass_vars$rec_id & v2 %in% pass_vars$rec_id)
  #2. If a comphet filter is configured, get the list of required vars
  if (!is.null(final_expr$comphet)) {
    comphet_required_vars <- as.data.frame(pass_vars %>%
                                             filter(! (!!final_expr$comphet)) )
  } else {
    comphet_required_vars <- pass_vars
  }
  #3. Get the final list of accepted comphet ids
  #from comphet select combo when: 
  #both vars in the combo passed the general variants filters (pass_comphet)
  #at least one var in the combo pass the comphet consequence filter (comphet_required_vars)
  pass_comphet <- pass_comphet %>% filter( v1 %in% comphet_required_vars$rec_id | v2 %in% comphet_required_vars$rec_id)
  
  return(list(vars=pass_vars$rec_id, comphet=pass_comphet$rec_id))
}

#Make quo expression combining all filters for genes and return PASS gene symbols
getPASSGenes_filters <- function(input, output, session, filters_settings, genes_scores) {
  #Build a var_operation variable identifying each variable default value
  var_operation <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_operation[[f]] <- filters_settings$DEFINITIONS[[n]][[f]][[2]]}
  }
  
  #Build a group_def variable with the group definition
  group_def <- list()
  for (n in names(filters_settings$GROUPS)) {
    group_def[[n]]$field <- filters_settings$GROUPS[[n]]$definition[[1]]
    group_def[[n]]$values <- unlist(filters_settings$GROUPS[[n]]$definition[[2]])
  }
  
  final_expr <- NULL
  for (g in names(filters_settings$GROUPS)) {
    #make a named list associated to the group: NAME=var_id, VALUE=value from control
    filters_values <- list() 
    for (v in filters_settings$GROUPS[[g]]$associated_values) {
      filterid <- paste0(v,"_",g)
      filters_values[[v]] <- input[[filterid]]
    }
    
    if (g == "global") {
      #global is a special group indicating filters that will be applied to all vars
      group_expr<- makeGroupExpression(logic = input[[paste0("LOGIC_",g)]],
                                       filters_operations = var_operation,
                                       filters_values = filters_values)
    } else {
      #for group filters a group_expr and a filter_expr are combined to evaluate vars only in the selected group
      groupid_expr <- quo((!!(as.name(group_def[[g]]$field))) %in% !!group_def[[g]]$values)
      groupfilters_expr <- makeGroupExpression(logic = input[[paste0("LOGIC_",g)]],
                                               filters_operations = var_operation,
                                               filters_values = filters_values)
      
      group_expr <- quo(!!groupid_expr & (!!groupfilters_expr))
    }
    if (is.null(final_expr)) {
      final_expr <- group_expr
    } else {
      final_expr <- quo(!!final_expr | (!!group_expr))
    }
  }
  
  #Genes passing all filters
  pass_genes <- as.data.frame(genes_scores %>% filter(! (!!final_expr)))
  
  return(pass_genes$gene)
}

#Return a list that can be used to make a JSON of the filters settings
getJSON_filters <- function(input, output, session, filters_settings) {
  filters_json <- list()
  ctrl_names <- names(input)
  remove <- grep("TXTSET|PRESETS", ctrl_names)
  ctrl_names <- ctrl_names[-remove]
  for (g in names(filters_settings$GROUPS)) {
    filters_json[[g]] <- list()
    group_string <- paste0("_",g)
    ctrls_group <- grep(group_string, ctrl_names)
    for (n in ctrl_names[ctrls_group]) {
      var_id <- gsub(group_string,"",n)
      filters_json[[g]][[var_id]] <- input[[n]]
    }
  }
  return(filters_json) 
}

#Read filters settings from a list of saved values
loadSettings_filters <- function(input, output, session, filters_values, filters_settings) {
  #Build a var_type variable identifying each variable type
  var_type <- list()
  for (n in names(filters_settings$DEFINITIONS)) {
    for (f in names(filters_settings$DEFINITIONS[[n]])) {var_type[[f]] <- n}
  }
  var_type[["LOGIC"]] <- "factors_fields"
  
  for (g in names(filters_values)) {
    for (v in names(filters_values[[g]])) {
      controlid <- paste0(v,"_",g)
      value <- filters_values[[g]][[v]]
      if (v %in% names(var_type)) {
        switch(var_type[[v]],
               numeric_fields = {updateSliderInput(session, inputId = controlid, value = as.numeric(value))},
               factors_fields = {updateSelectInput(session, inputId = controlid, selected = value)},
               binary_fields = {updateCheckboxInput(session, inputId = controlid, value = as.logical(value))}
        )
      } else {
        message("LOAD WARNING: ",v, " variable not found in the DEFINITIONS from your filters settings")
      }
    }
  }
}
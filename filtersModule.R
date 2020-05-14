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
source("filterObject.R")

#################
### FUNCTIONS ###
#################
`%nin%` = Negate(`%in%`)

intersectLists <- function(vector1, vector2) {
  if (length(vector1) == 0) {
    return(vector2)
  } else if (length(vector2) == 0) {
    return(vector1)
  } else {
    return(intersect(vector1, vector2))
  }
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

filtersVariantsUI <- function(id, filters_settings, variants_df, na_values, tooltips) {
  ns <- NS(id)
  
  #Build list containing type, label and default value for each variable as defined in settings
  vars_type <- list()
  vars_label <- list()
  vars_default <- list()
  vars_operation <- list()
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      vars_type[[v]] <- t
      vars_label[[v]] <- filters_settings$DEFINITIONS[[t]][[v]][[1]]
      vars_operation[[v]] <- filters_settings$DEFINITIONS[[t]][[v]][[2]]
      vars_default[[v]] <- getDefaultValue(var_name = v,
                                           ctrl_type = t,
                                           ctrl_value = filters_settings$DEFINITIONS[[t]][[v]][[3]],
                                           na_values = na_values,
                                           df = variants_df)
    }
  }

  #Create one object per group
  variants_boxes <- NULL
  for (g in names(filters_settings$GROUPS)) {
    grp_vars <- unlist(filters_settings$GROUPS[[g]][["associated_values"]])
    grp_vars_type <- vars_type[grp_vars]
    grp_vars_default <- vars_default[grp_vars]
    grp_vars_label <- vars_label[grp_vars]
    group_box <- filterObjectUI(id = ns(g),
                                group_name = g,
                                group_definition = filters_settings$GROUPS[[g]][["definition"]], 
                                group_presets = filters_settings$PRESETS[[g]],
                                vars_type = grp_vars_type, 
                                vars_default = grp_vars_default, 
                                vars_label = grp_vars_label, 
                                variants_df = variants_df, 
                                na_values = na_values, 
                                factorsmap = filters_settings$FACTORSMAP,
                                tooltips = tooltips)
    variants_boxes <- tagList(variants_boxes, group_box)  
  }
  
  return(variants_boxes)
}

##############
### SERVER ###
##############

#Manage the update of filters values when preset or text input change
observeFilters <- function(input, output, session, filters_settings, variants_df, na_values, tooltips) {
  for (g in names(filters_settings$GROUPS)) {
      callModule(observeInputs, g, 
               variants_df = variants_df,
               na_values = na_values,
               filters_settings = filters_settings)
  }
}

#Make quo expression combining all filters for variants and return PASS variants ids
getPASSVars_filters <- function(input, output, session, filters_settings, variants_df, comphet_df) {
  #Make list of vars operations from config
  vars_operation <- list()
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      vars_operation[[v]] <- filters_settings$DEFINITIONS[[t]][[v]][[2]]
    }
  }
  
  #Apply filters for each group
  PASS_group <- list()
  comphet_expr <- NULL
  for (group_name in names(filters_settings$GROUPS)) {
    group_definition <- filters_settings$GROUPS[[group_name]]$definition
    group_vars <- unlist(filters_settings$GROUPS[[group_name]]$associated_values)
    grp_vars_op <- vars_operation[group_vars]
    filter_expr <- callModule(getFilterExpression, group_name,
                              group_name = group_name,
                              group_definition = group_definition,
                              group_vars = group_vars,
                              vars_operation = grp_vars_op)
    if (group_name == "comphet") {
      comphet_expr <- filter_expr
    } else {
      PASS_group[[group_name]] <- as.data.frame(variants_df %>% filter(!!filter_expr))$rec_id
      #message("FILTER EXPR ", group_name)
      #cat(as.character(filter_expr))
      #message("VARS PASS ", group_name, " ", length(PASS_group[[group_name]]))
      #message(paste(PASS_group[[group_name]], collapse=";"))
    }
  }
  
  #SINGLE VARIANTS
  #Variants passing all filters
  if ("global" %in% names(PASS_group)) {
    global_pass <- PASS_group[["global"]]
    #message("PASS VARS IN GLOBAL VARIABLE ", length(global_pass))
    PASS_group["global"] <- NULL
    pass_vars_ids <- unique(unlist(PASS_group)) #Merge vars for each group
    #message("PASS VARS AFTER UNLIST ", length(pass_vars_ids))
    pass_vars_ids <- intersectLists(pass_vars_ids, global_pass) #Intersect with vars passing global filter
    #message("PASS VARS AFTER INTERSECT ", length(pass_vars_ids))
  } else {
    pass_vars_ids <- unique(unlist(PASS_group)) #Merge vars for each group
  }
  #message("FINAL LIST PASS VARS FROM MODULE ", length(pass_vars_ids))
  pass_vars <- as.data.frame(variants_df %>% filter(rec_id %in% pass_vars_ids))
  
  #COMPHET VARS
  #1. Select comphet where both vars pass the variants filters
  pass_comphet <- comphet_df %>% filter(v1 %in% pass_vars$rec_id & v2 %in% pass_vars$rec_id)
  #2. If a comphet filter is configured, get the list of required vars
  if (!is.null(comphet_expr)) {
    comphet_required_vars <- as.data.frame(pass_vars %>%
                                          filter(!!comphet_expr) )
  } else {
    comphet_required_vars <- pass_vars
  }
  #3. Get the final list of accepted comphet ids
  #from comphet select combo when: 
  #both vars in the combo passed the general variants filters (pass_comphet)
  #at least one var in the combo pass the comphet consequence filter (comphet_required_vars)
  pass_comphet <- pass_comphet %>% filter( v1 %in% comphet_required_vars$rec_id | v2 %in% comphet_required_vars$rec_id)
  #message("FINAL LIST PASS COMPHET FROM MODULE ", length(pass_comphet$rec_id))
  return(list(vars=pass_vars$rec_id, comphet=pass_comphet$rec_id))
}

#Make quo expression combining all filters for genes and return PASS gene symbols
getPASSGenes_filters <- function(input, output, session, filters_settings, genes_scores) {
  #Make list of vars operations from config
  vars_operation <- list()
  for (t in names(filters_settings$DEFINITIONS)) {
    for (v in names(filters_settings$DEFINITIONS[[t]])) {
      vars_operation[[v]] <- filters_settings$DEFINITIONS[[t]][[v]][[2]]
    }
  }
  
  #Apply filters for each group
  PASS_group <- list()
  for (group_name in names(filters_settings$GROUPS)) {
    group_definition <- filters_settings$GROUPS[[group_name]]$definition
    group_vars <- unlist(filters_settings$GROUPS[[group_name]]$associated_values)
    grp_vars_op <- vars_operation[group_vars]
    filter_expr <- callModule(getFilterExpression, group_name,
                              group_name = group_name,
                              group_definition = group_definition,
                              group_vars = group_vars,
                              vars_operation = grp_vars_op)
      PASS_group[[group_name]] <- as.data.frame(genes_scores %>% filter(!!filter_expr))$gene
  }
  
  #Genes passing all filters
  if ("global" %in% names(PASS_group)) {
    global_pass <- PASS_group[["global"]]
    PASS_group["global"] <- NULL
    pass_genes_ids <- unique(unlist(PASS_group)) #Merge vars for each group
    pass_genes_ids <- intersectLists(pass_genes_ids, global_pass) #Intersect with vars passing global filter
  } else {
    pass_genes_ids <- unique(unlist(PASS_group)) #Merge vars for each group
  }
  
  #Genes passing all filters
  pass_genes <- as.data.frame(genes_scores %>% filter(gene %in% pass_genes_ids))
  
  return(pass_genes$gene)
}

#Return a data frame with the filters values
getDF_filters <- function(input, output, session, filters_settings) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  for (group_name in names(filters_settings$GROUPS)) {
    group_df <- callModule(getDF, group_name,
                           group_name = group_name)
    filters_df <- rbind(filters_df, group_df)
  }
  filters_df <- filters_df[order(filters_df$group, filters_df$filter),]
  return(filters_df)
}

#Return a list that can be used to make a JSON of the filters settings
getJSON_filters <- function(input, output, session, filters_settings) {
  filters_json <- list()
  for (group_name in names(filters_settings$GROUPS)) {
    group_json <- callModule(getJSON, group_name)
    filters_json[[group_name]] <- group_json
  }
  return(filters_json) 
}

#Read filters settings from a list of saved values
loadSettings_filters <- function(input, output, session, filters_settings, filters_json) {
  for (group_name in names(filters_settings$GROUPS)) {
    callModule(loadSettings, group_name, 
               filters_settings = filters_settings,
               filters_json = filters_json)
  }
}
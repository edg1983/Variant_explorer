########################
### GENTOYPES MODULE ###
########################

RV <- reactiveValues()

#GT_cols_affected, GT_cols_unaffected are names of columns reporting genotypes for affected / unaffected
genotypeUI <- function(id, GT_cols_affected, GT_cols_unaffected) {
  ns <- NS(id)
  
  #Generate 2 groups of selectInputs for affected and unaffected samples
  affected_controls <- NULL
  for (id in GT_cols_affected) {
    affected_controls <- tagList(affected_controls,
                                 fluidRow(
                                   column(8,selectInput(ns(id),label = paste0(id, " genotype"), choices = c(0,1,2), selected = 0)),
                                   column(4,checkboxInput(ns(paste0(id,"_allowNA")), label= "allow missing",value = F))
                                   ))
  }
  
  unaffected_controls <- NULL
  for (id in GT_cols_unaffected) {
    unaffected_controls <- tagList(unaffected_controls,
                                   fluidRow(
                                     column(8,selectInput(ns(id),label = paste0(id, " genotype"), choices = c(0,1,2), selected = 0)),
                                     column(4,checkboxInput(ns(paste0(id,"_allowNA")), label= "allow missing",value = F))
                                   ))
  }
  
  output <- tagList(
    fluidRow(
      column(6, h4("Affected genotypes"), affected_controls),
      column(6, h4("Unaffected genotypes"), unaffected_controls)),
    )
}

genotypeModule <- function(input, output, session, vars_df) {
  ctrl_names <- names(input)
  NA_controls <- grep("_allowNA",names(input))
  ctrl_names <- ctrl_names[-NA_controls]
  
  final_expr <- NULL
  for (n in ctrl_names) {
    if (input[[paste0(n,"_allowNA")]] == T) {
      myexpr <- expr((!!as.name(n) == !!input[[n]] | is.na(!!as.name(n))))
    } else {
      myexpr <- expr(!!as.name(n) == !!input[[n]])
    }
    if (is.null(final_expr)) {
      final_expr <- myexpr
    } else {
      final_expr <- expr(!!final_expr & !!myexpr)
    }
    message(final_expr)
    
    vars_list <- vars_df %>% filter(!!final_expr) %>% pull(rec_id)
  }
  
  return(vars_list)
}

getDF_genotype <- function(input, output, session) {
  filters_df <- data.frame(group=character(), filter=character(), value=character(), stringsAsFactors = F)
  ctrl_names <- names(input)

    for (n in ctrl_names) {
      newline <- c("genotype", n, input[[n]])
      filters_df[nrow(filters_df)+1,] <- newline
    }
  filters_df <- filters_df[order(filters_df$filter),]
  return(filters_df)
}

getJSON_genotype <- function(input, output, session) {
  filters_json <- list()
  ctrl_names <- names(input)
  
  for (n in ctrl_names) {
    filters_json[[n]] <- input[[n]]
  }
  
  return(filters_json) 
}

loadSettings_genotype <- function(input, output, session, filters_values) {
  for (v in names(filters_values)) {
    if (v %in% names(input)) {
      if (grepl("_allowNA", v)) {
        updateCheckboxInput(session, inputId = v, value = filters_values[[v]])
      } else {
        updateSelectInput(session, inputId = v, selected = filters_values[[v]])
      }
    } else {
      message("LOAD WARNING: ",v, " control not found in the genotype module")
    }
  }
}
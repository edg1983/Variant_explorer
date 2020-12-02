########################
### DOWNLOAD MODULE  ###
########################

# downloadObj function usage
# Output_prefix <- string to prefix for all files or output file for single output
# Output_data <- named list of suffixes and associated data.frames or single data.frame 
# zip_archive <- NULL to save individual files, string with a file name to compress all files and save a zip archive
# data frames are saved only if nrow > 0

# Examples #
# 1. single file output 
# callModule(downloadObj, id="download", output_prefix="myprefix.tsv", output_data=data.frame)
# ==> myprefix.tsv
# 2. single file, zip archive output 
# callModule(downloadObj, id="download", output_prefix="myprefix.tsv", output_data=data.frame, zip_archive="myarchive.zip")
# ==> myarchive.zip [myprefix.tsv]
# 4. multiple files compressed in zip archive
# callModule(downloadObj, id="download", output_prefix="myprefix", output_data=list("file1.tsv" = data.frame1, "file2.tsv" = data.frame2), zip_archive="myarchive.zip")
# ==> myarchive.zip [myprefix.file1.tsv, myprefix.file2.tsv]

saveData <- function(output_data, filename, sep=NULL, row_names=F, quote=F, col_names = T) {
  switch(str_sub(filename, start= -4),
         ".tsv" = { 
           if(is.null(sep)) {sep="\t"}
           if (nrow(output_data) > 0) {
             write.table(output_data, file=filename, sep=sep, row.names=row_names, quote=quote, col.names = col_names)
           } },
         ".csv" = { 
           if(is.null(sep)) {sep=","}
           if (nrow(output_data) > 0) {
             write.table(output_data, file=filename, sep=sep, row.names=row_names, quote=quote, col.names = col_names)
           } },
         ".xml" = { 
           sep=" "
           write.table(output_data, file=filename, sep=sep, row.names=row_names, quote=quote, col.names = col_names)
           },
         "json" = { write_json(output_data, path=filename, pretty=T, auto_unbox=T) })
}

downloadObjUI <- function(id, label = "Download") {
  ns <- NS(id)
  downloadButton(ns("download_results"), label = label)
}

downloadObj <- function(input, output, session, output_prefix, output_data, sep=NULL, zip_archive=NULL, col_names=TRUE) {
  output$download_results <- downloadHandler( 
    filename = function() { 
      if (is.null(zip_archive)) {
        output_prefix
      } else {
        zip_archive
      }
    },
    content = function(file) {
      #shiny::validate(need(
      #  (output_prefix != "" & !is.null(output_prefix)) & 
      #    (inherits(output_data, "list") | inherits(output_data,"data.frame")), FALSE))
      if (is.null(zip_archive)) {
        saveData(output_data, file=file, sep=sep, col_names = col_names)
      } else {
        tmp_prefix <- paste0(tempdir(), "/", output_prefix)
        if (inherits(output_data, "list") == TRUE) {
          files <- NULL
          for (suffix in names(output_data)) {

            tmp_file <- paste0(tmp_prefix,".", suffix)
            saveData(output_data[[suffix]], file=tmp_file, sep=sep, col_names = col_names)
            #write.table(output_data[[suffix]], file=tmp_file, sep=sep, row.names=F, quote=F)
            files <- c(files, tmp_file)

          }
          zip(file,files,flags = "-j")
        } else {
          saveData(output_data, file=tmp_prefix, sep=sep, col_names = col_names)
          #write.table(output_data, file=tmp_prefix, sep=sep, row.names=F, quote=F)
          zip(file,tmp_prefix,flags = "-j")
        }
      }
    }
  )
}

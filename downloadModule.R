########################
### DOWNLOAD MODULE  ###
########################

# downloadObj function usage
# Output_prefix <- string to prefix for all files or output file for single output
# Output_data <- named list of suffixes and associated data.frames or single data.frame 
# zip_archive <- NULL to save individual files, string with a file name to compress all files and save a zip archive
# data frames are saved only if nrow > 0

# Examples #
# 1. single file tsv output 
# callModule(downloadObj, id="download", output_prefix="myprefix.tsv", output_data=data.frame)
# ==> myprefix.tsv
# 2. single file tsv, zip archive output 
# callModule(downloadObj, id="download", output_prefix="myprefix.tsv", output_data=data.frame, zip_archive="myarchive.zip")
# ==> myarchive.zip [myprefix.tsv]
# 3. multiple tsv files
# callModule(downloadObj, id="download", output_prefix="myprefix", output_data=list("file1.tsv" = data.frame1, "file2.tsv" = data.frame2)) 
# ==> myprefix.file1.tsv, myprefix.file2.tsv
# 4. multiple tsv files compressed in zip archive
# callModule(downloadObj, id="download", output_prefix="myprefix", output_data=list("file1.tsv" = data.frame1, "file2.tsv" = data.frame2), zip_archive="myarchive.zip")
# ==> myarchive.zip [myprefix.file1.tsv, myprefix.file2.tsv]

downloadObjUI <- function(id, label = "Download") {
  ns <- NS(id)
  downloadButton(ns("download_results"), label = label)
}

downloadObj <- function(input, output, session, output_prefix, output_data, sep="\t", zip_archive=NULL, col_names=TRUE) {
  output$download_results <- downloadHandler( 
    filename = function() { 
      if (is.null(zip_archive)) {
        output_prefix
      } else {
        zip_archive
      }
    },
    content = function(file) {
      
      if (is.null(zip_archive)) {
        if (inherits(output_data, "list") == TRUE) {
          for (suffix in names(output_data)) {
           out_file <- paste0(file, ".", suffix)
            if (nrow(output_data[[suffix]]) > 0) {
              write.table(output_data[[suffix]], file=out_file, sep=sep, row.names=F, quote=F, col.names = col_names)
            }
          }
        } else {
          write.table(output_data, file=file, sep=sep, row.names=F, quote=F, col.names = col_names)
        }
      } else {
        tmp_prefix <- paste0(tempdir(), "/", output_prefix)
        if (inherits(output_data, "list") == TRUE) {
          files <- NULL
          for (suffix in names(output_data)) {
            if (nrow(output_data[[suffix]]) > 0) {
            tmp_file <- paste0(tmp_prefix,".", suffix)
            write.table(output_data[[suffix]], file=tmp_file, sep=sep, row.names=F, quote=F)
            files <- c(files, tmp_file)
            }
          }
          zip(file,files,flags = "-j")
        } else if (inherits(output_data, "data.frame")) {
          write.table(output_data, file=tmp_prefix, sep=sep, row.names=F, quote=F)
          zip(file,tmp_prefix,flags = "-j")
        }
      }
    }
  )
}

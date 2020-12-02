# Encrypt any RData object with a given password
# use readRDS, so data object are expected to contain a single element saved with saveRDS

# Symmetric sodium encryption and decryption of a datafile using cyphr
# https://github.com/ropensci/cyphr


library(cyphr)

# encrypt an R object
encrypt_datafile = function(myobj, outf, pwd) {
  k = cyphr::key_sodium(sodium::hash(charToRaw(pwd)))
  # encrypt data file
  result <- tryCatch({
    cyphr::encrypt(saveRDS(myobj, outf, compress = TRUE), k)
    return(1)
  }, error=function(cond) {
    #message("Encryption error. Unable to generate ", outf)
    if (file.exists(outf)) { file.remove(outf) }
    return(0)
  })
  
  return (result)
}

# decrypt to an r object
decrypt_datafile = function(inf, outf, pwd) {
  #pwd = .rs.askForPassword("Enter password") # asks for pwd in RStudio
  k = cyphr::key_sodium(sodium::hash(charToRaw(pwd)))
  d = tryCatch({
    cyphr::decrypt(readRDS(inf), k)
  }, error=function(cond) {
    message("Could not decrypt! pwd wrong?")
    if (file.exists(outf)) { file.remove(outf) }
    return(NULL)
  })
  return (d)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args)<3) {
  stop("Usage as follow:
       Encrypt_data.R input_dir output_dir password
       All .RData files in input_dir will be encrypted and saved to output_dir", 
       call.=FALSE)
}

input_dir <- args[1]
output_dir <- args[2]
pwd <- args[3]

if (!file_test("-d", args[1])) {
  stop("input directory not found!", call.=FALSE)
}  

if (!file_test("-d", args[2])) {
  tryCatch({
    dir.create(output_dir)
    message("Output folder do not exists and will be created")
  }, error=function(cond) {
    stop("Unable to create output folder", call.=FALSE)
  })
}

message("Output folder: ", output_dir) 

filelist <- list.files(input_dir,pattern=".RData")
nfiles <- length(filelist)
if (nfiles > 0) {
  message(nfiles, " files found in input directory")
} else {
  stop("No RData files found in the input directory", call. = FALSE)
}

message("Encrypting data...")
saved_files <- 0
failed_files <- 0
skipped_files <- 0
n <- 0
pb <- txtProgressBar(min=0, max=nfiles,style = 3)  
for (f in filelist) {
  n <- n + 1
  setTxtProgressBar(pb,value = n)
  
  #load data
  inputfile <- paste0(input_dir, "/", f)
  rdata <- readRDS(inputfile)
  
  #Save encrypted object
  out_file <- paste0(output_dir, "/", f, ".enc")
  if (!file_test("-f", out_file)) {
    save_results <- encrypt_datafile(rdata,outf=out_file,pwd = pwd)
    if (save_results == 1) {
      saved_files = saved_files + 1
    } else {
      failed_files = failed_files + 1
    } 
  } else {
    skipped_files = skipped_files + 1
  }
}

message("\nProcess completed\n",
        saved_files, " saved\n", 
        failed_files, " failed\n", 
        skipped_files, " output files already present")

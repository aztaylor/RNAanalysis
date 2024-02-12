library("tximport")
library("readr")

cwd = getwd()

dirs = list.files("./quants/", include.dirs=TRUE)
data <- list()
for (dir in dirs) {
  if (endsWith(dir, "quant")) {
    if (file_test("-f", paste(cwd, "quants", dir, "quant.sf", sep = "/"))){
      txi <- tximport(paste(cwd, "quants", dir, "quant.sf", sep = "/"), type = "salmon", txOut = TRUE)
      
      # Establish the file name convention, isolating the experiemtnal condititions.
      pattern <- "(?<=ATEY090523-).+(?=_L)"
      target <- dir
      match <- regexec(pattern, target, perl = TRUE)
      
      # Extract conditions from the file name
      name <- regmatches(target, match)[[1]]
      cond <- substr(name, 1, nchar(name)-2)
      
      # Set the conditions within the txi list to include the experimental conditions.
      nelems = length(txi$counts)
      cond_vec = rep(c(cond),each=nelems)
      txi[["condition"]] <- cond_vec
      data = append(data, txi)
    }
  }
}
print(data)

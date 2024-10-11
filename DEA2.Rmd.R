```{r}
library("tximport")
library("readr")
library("DESeq2")
```
```{r}
# First lets list all of the quant.sf files
cwd <- getwd()
quant_fp <- paste(cwd, 'quants/', sep = '/')

quant_files <- list.files(path = quant_fp, pattern="quant.sf$", recursive = TRUE)
```
```{r}
# And we can then extranct the sample names.
quant_dirs <- list.files(quant_fp, pattern = "_quant$", full.names = TRUE)
sample_names <- gsub("_quant$", "", basename(quant_dirs))
# correct for the FASTQ generation file included in the directory.
sample_names <- sample_names[1: 48]
```
```{r}
# And now we can create the txi import, note that it can take an array of file names.
txi <- tximport(paste(cwd, "quants", quant_files, sep='/'), type = "salmon", txOut = TRUE)
```
```{r}
# Here we will create a temporary Dataframe for the TXI data so that we can add the sample names as columns.
# Before feeding to DESeq2 we need to convert back to a numeric matrix.
data = list()
txi_length = length(txi)
for (i in 1:(txi_length-1)) {
  temp <- as.data.frame(txi[[i]])
  colnames(temp) <- sample_names
  data[[i]] <- as.matrix(temp)
}
```
```{r}
# Exract the conditions from the file names supplied by Jen
conditions = sub(pattern ="-notebooks
                 *_L", sample_names[1])
conditions = c("NCt1", "NCt1","NCt1","NCt1","NCt1ind", "NCt1ind","NCt1ind","NCt1ind", "NCt2","NCt2","NCt2","NCt2","NCt2ind","NCt2ind","NCt2ind","NCt2ind",
               "s4t1ind","s4t1ind","s4t1ind","s4t1ind","s4t1","s4t1","s4t1","s4t1", "s4t2","s4t2","s4t2", "s4t2", "WT")
dds <- DESeqDataFromTximport(data, coldata=smaple_names, design=conditions )
```
```
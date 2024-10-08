# Data for module 1 scRNAseq

The data has been charged in a [folder](https://www.google.com/url?q=https%3A%2F%2Fdrive.google.com%2Fdrive%2Ffolders%2F1GuEwAM6xPMdyry-4s77N_6vHiDks6jGC%3Fusp%3Dsharing) in drive. If you get lost at any point in the code, look in this folder for the processed data you need and load it to be able to continue with the course. If you want to load them directly in your colab, you can use the next code, that is already in **modules\Module1_scRNAseq\scRNAseq.ipynb** code in the section **GEO DATA/Download processed data**

```
# Function to execute shell commands in Google Colab when running R kernel
shell_call <- function(command, ...) {
  result <- system(command, intern = TRUE, ...)
  cat(paste0(result, collapse = "\n"))
}

# Download data from drive
shell_call("gdown --fuzzy https://drive.google.com/drive/folders/1GuEwAM6xPMdyry-4s77N_6vHiDks6jGC?usp=sharing --folder -O Module1_QC_Clust_DEG_CellAnot")
```

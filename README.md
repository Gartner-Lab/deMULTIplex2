# deMULTIplex2

## Installation

```
install.packages("devtools") 
devtools::install_github("Gartner-Lab/deMULTIplex2") # install
library(deMULTIplex2) 
```

## Starting with tag count matrix
Must have a count matrix where rows are individual cells and columns represent unique sample tags. You can load an example tag matrix from Stoeckius et al. by calling `data(stoeckius_pbmc);tag_mtx <- stoeckius_pbmc`.

```
res <- demultiplexTags(tag_mtx,
                       plot.path = "~/",
                       plot.name = "test",
                       plot.diagnostics = FALSE)
table(res$final_assign)
                      
```

## Starting with raw sequencing data (Illumina FASTQs)
Make sure you have your barcode library sequenced and the reads saved in FASTQ format (can be gzipped).
```
~/Experiment2
| Exp2MULTI_S3_L002_R1_001.fastq.gz
| Exp2MULTI_S3_L002_R2_001.fastq.gz
```
Provide readTags() with the location of the files, the prefix of the FASTQ file names for the library you want to process, and the type of barcode and assay you used. You may also provide a vector of cell barcodes (i.e. from the barcodes.tsv file output by cellranger) to pre-filter your barcode reads.
```
read_table <- readTags(dir = "~/Experiment2",
                       name = "Exp2MULTI",
                       barcode.type = "MULTIseq",
                       assay = "RNA",
                       filter.cells = exp2_cells)
```

Next, alignTags() will take this read table and count the number of UMIs detected per tag, per cell. Sample tag reads are error-corrected by aligning them to a provided vector of tag sequences used in the experiment. You can manually supply these sequences, or they can be subset from the full vector of MULTIseq barcodes provided with this package. 

```
tag.ref <- tags[1:24]
# tag.ref <- tags[c(paste("A",1:12,sep=""),paste("B",1:12,sep=""))]

tag_mtx <- alignTags(read_table,
                     tag.ref)
```
Then, demultiplexTags() can be run to classify each cell to its sample of origin.

```
res <- demultiplexTags(tag_mtx,
                       plot.path = "~/",
                       plot.name = "test",
                       plot.diagnostics = FALSE)
                       
table(res$final_assign)
                      
```

## Visualization Tools

```
tagHist(tag_mtx,
        minUMI = 10)
```
```
tagCallHeatmap(tag_mtx,
           res$final_assign)
```

## Troubleshooting

* Installation failed on macOS - You may need to install Xquartz (https://www.xquartz.org/) first.

## Cite deMULTIplex2

To be added.

## License

This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

## References

1. Stoeckius, M., Zheng, S., Houck-Loomis, B., Hao, S., Yeung, B. Z., Mauck, W. M., Smibert, P., & Satija, R. (2018). Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome biology, 19(1), 1-12. 




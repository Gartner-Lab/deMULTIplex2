# deMULTIplex2

deMULTIplex2 is a mechanism-guided classification algorithm for multiplexed scRNA-seq data that successfully recovers many more cells across a spectrum of challenging datasets compared to existing methods. deMULTIplex2 is built on a statistical model of tag read counts derived from the physical mechanism of tag cross-contamination. Using generalized linear models and expectation-maximization, deMULTIplex2 probabilistically infers the sample identity of each cell and classifies singlets with high accuracy. 

## Installation

```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

if (!require("devtools", quietly = TRUE))
    install.packages("devtools")  # devtools must be correctly installed to install deMULTIplex2 from github source.

devtools::install_github("Gartner-Lab/deMULTIplex2") # install deMULTIplex2
library(deMULTIplex2) # load deMULTIplex2 into session
```

It is **recommended** (not required) to install ggrastr using `install.packages("ggrastr")`, because summary and diagnostic plots may have a lot of points (cells). deMULTIplex2 will use ggrastr to plot if it is correctly installed. 

## Starting with tag count matrix
**`demultiplexTags()`** is the core function of deMULTIplex2. User must provide a tag count matrix where rows are individual cells and columns represent unique sample tags. You can load an example tag matrix from Stoeckius et al. (2018) by calling `data(stoeckius_pbmc);tag_mtx <- stoeckius_pbmc`.

```
res <- demultiplexTags(tag_mtx, # Required, the tag count matrix from your experiment, can be either dense or sparse
                       plot.path = "~/", # Where to output a summary plot
                       plot.name = "test", # text append to the name of the summary plot file
                       plot.diagnostics = FALSE) # Whether to output diagnostics plots for each tag
table(res$final_assign)
```

For 10x Genomics 3’ CellPlex Multiplexing datasets, the latest version of Seurat can grab the CMO matrix from the raw matrix folder:
```
library(Seurat);library(Matrix)
seu1 <- Read10X("lib1/outs/multi/count/raw_feature_bc_matrix/") # Change this to your cellranger output path
tag_mtx <- t(seu1$`Multiplexing Capture`[, cells_pass_filter]) # Suggest define and filter cells using transcriptome prior to demultiplexing
```
Alternatively, just grep features with "CMO" in its name to get the tag count matrix.


## Starting with raw sequencing data (Illumina FASTQs)
Make sure you have your barcode library sequenced and the reads saved in FASTQ format (can be gzipped).
```
~/Experiment2
| Exp2MULTI_S3_L002_R1_001.fastq.gz
| Exp2MULTI_S3_L002_R2_001.fastq.gz
```
Provide **`readTags()`** with the location of the files, the prefix of the FASTQ file names for the library you want to process, and the type of barcode and assay you used. You may also provide a vector of cell barcodes (i.e. from the barcodes.tsv file output by cellranger) to pre-filter your barcode reads.
```
read_table <- readTags(dir = "~/Experiment2",
                       name = "Exp2MULTI",
                       barcode.type = "MULTIseq",
                       assay = "RNA",
                       filter.cells = exp2_cells)
```

Next, **`alignTags()`** will take this read table and count the number of UMIs detected per tag, per cell. Sample tag reads are error-corrected by aligning them to a provided vector of tag sequences used in the experiment. You can manually supply these sequences, or they can be subset from the full vector of MULTI-seq barcodes provided with this package. 

```
data(multiseq_oligos) # Current MULTI-seq oligo sequence provided with the package
tag.ref <- multiseq_oligos[1:24] # Assuming the first 24 tags are used
tag_mtx <- alignTags(read_table, tag.ref)
```

The produced tag matrix can then be used as input for the `demultiplexTags` function. It is recommended to pre-filter the matrix to remove majority of the empty droplets for robust classification, i.e.,
```
tag_mtx = tag_mtx[cell_barcodes, ] # cell_barcodes can be determined using the transcriptome data, or by setting a minimum count threshold
```

## Visualization tools

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

Qin, Z., Danny, N. C., & Zev Jordan, G. (2023). Demultiplex2: robust sample demultiplexing for scRNA-seq. bioRxiv, 2023.04.11.536275. https://doi.org/10.1101/2023.04.11.536275 

## License

This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.





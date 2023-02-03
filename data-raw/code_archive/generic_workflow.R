### MULTI-seq or MULTI-ATAC Demultiplexing Workflow ###

source("~/Desktop/deMULTIplex2.R")

barcodes <- readxl::read_excel("~/Desktop/MULTI-seq_oligos_Nov2019.xlsx")$`Barcode Sequence`
barcodes <- barcodes[1:10] # select barcodes used in experiment

# Parse raw reads to build read table
readTable <- readMULTI(dir = "~/Desktop/Experiment/FASTQs", # directory where your barcode library FASTQ files are kept
                       name = "MULTI", # prefix of your barcode library FASTQ files 
                       assay = "RNA", # "RNA", "ATAC", or "Multiome"
                       barcode = "MULTIseq", # "MULTIseq" or "MULTI-ATAC" 
                       filterCells = NULL) # NULL or a vector of cell barcodes

saveRDS(readTable, "~/Desktop/Experiment/MULTI/readTable.rds")

# Tally read table to build cell x barcode count matrix
barTable <- alignMULTI(readTable,
                       barcodes,
                       filterCells = NULL)

saveRDS(barTable, "~/Desktop/Experiment/MULTI/barTable.rds")

# Visualize distribution of each barcode
barHist(barTable, 
        minUMI = 10,
        plotnUMI = T)

# Filter barTable by nUMI if necessary to remove empty droplets or very poorly labeled cells
barTable_clean <- barTable[barTable$nUMI >= 150, ]

# Run classification algorithm
calls <- classifyLoop(barTable_clean,
                      counts = "raw", # "raw", "sct_count", or "sct_res"
                      doub.ident = F) # whether to keep unique multiple types or group them as "Doublet"



# Spot-check classification with a heatmap
barHeatmap(barTable,
           calls,
           log = T,
           colHigh = "dark red")

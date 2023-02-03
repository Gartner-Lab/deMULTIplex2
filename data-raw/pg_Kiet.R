


library(deMULTIplex2)

mb_mtx <- readRDS("data-raw/MB_matrices/Kiet_T1_mb_mtx.rds")

res <- classifyMULTI(mb_mtx,
              posThresh = 0.2,
              cosineThresh = seq(0.5,0.9,0.1),
              plotUMAP = "cUMI_Cos",
              plotDiagnostics = T,
              plotPath = "data-raw/plots/",
              UMAPNeighbors = 30L,
              seed = 1)



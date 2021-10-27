## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 100),
  tidy = TRUE
)
library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------
tab <- read.csv(system.file("extdata", "dataset_size.csv", package = "morphoBlocks"))

knitr::kable(tab, caption = "Table 1. Number of specimens <i>(i)</i> and blockcs <i>(j)</i> in each of the model datasets.", align = "c", table.attr = "style='width:50%;'") %>%
  kableExtra::kable_classic(full_width = TRUE) %>%
  kableExtra::row_spec(0, bold = TRUE)

## -----------------------------------------------------------------------------
library(morphoBlocks)

# Start time
start_time <- Sys.time()

# Generate dataset 1
set.seed(1)

# Specify the number of samples
i <- 10

# Specify the number of blocks
j <- 10 

# Create an empty list in which to store the dataset
dataset.1 <- list()

# Generate the blocks for the dataset
for(k in 1 : j){
dataset.1[[k]] <- dodecBlock(n = i)
}
# Duration
end_time <- Sys.time()
end_time - start_time


# Repeat for each of the datasets
start_time <- Sys.time()
set.seed(1)
i <- 50
j <- 50 
dataset.2 <- list()
for(k in 1 : j){
dataset.2[[k]] <- dodecBlock(n = i)
}
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
set.seed(1)
i <- 200
j <- 10 
dataset.3 <- list()
for(k in 1 : j){
dataset.3[[k]] <- dodecBlock(n = i)
}
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
set.seed(1)
i <- 50
j <- 206 
dataset.4 <- list()
for(k in 1 : j){
dataset.4[[k]] <- dodecBlock(n = i)
}
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
set.seed(1)
i <- 200
j <- 206 
dataset.5 <- list()
for(k in 1 : j){
dataset.5[[k]] <- dodecBlock(n = i)
}
end_time <- Sys.time()
end_time - start_time

## -----------------------------------------------------------------------------
# Format and scale dataset 1
start_time <- Sys.time()
block.list.1 <- combineBlocks(blocks = dataset.1)
end_time <- Sys.time()
end_time - start_time

# Repeat for each of the datasets
start_time <- Sys.time()
block.list.2 <-combineBlocks(blocks = dataset.2)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
block.list.3 <-combineBlocks(blocks = dataset.3)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
block.list.4 <-combineBlocks(blocks = dataset.4)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
block.list.5 <-combineBlocks(blocks = dataset.5)
end_time <- Sys.time()
end_time - start_time

## -----------------------------------------------------------------------------
# Perform PCA on dataset 1
start_time <- Sys.time()
result.pca.1 <- analyseBlocks(block.list.1, option = "pca")
end_time <- Sys.time()
end_time - start_time

# Perform RCPCA on dataset 1
start_time <- Sys.time()
result.rcpca.1 <- analyseBlocks(block.list.1, option = "rcpca", ncomp = 3)
end_time <- Sys.time()
end_time - start_time

# Plot the scores from each analysis
layout(matrix(c(1,2), ncol = 2))
scoresPlot(result.pca.1, comp = c(1, 2), plabels = 1:block.list.1@n[1])
scoresPlot(result.rcpca.1, comp = c(1, 2), plabels = 1:block.list.1@n[1], consensus.only = TRUE)

# Repeat for dataset 2
start_time <- Sys.time()
result.pca.2 <- analyseBlocks(block.list.2, option = "pca")
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
result.rcpca.2 <- analyseBlocks(block.list.2, option = "rcpca", ncomp = 3)
end_time <- Sys.time()
end_time - start_time

layout(matrix(c(1,2), ncol = 2))
scoresPlot(result.pca.2, comp = c(1, 2), plabels = 1:block.list.2@n[1])
scoresPlot(result.rcpca.2, comp = c(1, 2), plabels = 1:block.list.2@n[1], consensus.only = TRUE)

# Repeat for dataset 3
start_time <- Sys.time()
result.pca.3 <- analyseBlocks(block.list.3, option = "pca")
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
result.rcpca.3 <- analyseBlocks(block.list.3, option = "rcpca", ncomp = 3)
end_time <- Sys.time()
end_time - start_time

layout(matrix(c(1,2), ncol = 2))
scoresPlot(result.pca.3, comp = c(1, 2), plabels = 1:block.list.3@n[1])
scoresPlot(result.rcpca.3, comp = c(1, 2), plabels = 1:block.list.3@n[1], consensus.only = TRUE)

# Repeat for dataset 4
start_time <- Sys.time()
result.pca.4 <- analyseBlocks(block.list.4, option = "pca")
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
result.rcpca.4 <- analyseBlocks(block.list.4, option = "rcpca", ncomp = 3)
end_time <- Sys.time()
end_time - start_time

layout(matrix(c(1,2), ncol = 2))
scoresPlot(result.pca.4, comp = c(1, 2), plabels = 1:block.list.4@n[1])
scoresPlot(result.rcpca.4, comp = c(1, 2), plabels = 1:block.list.4@n[1], consensus.only = TRUE)

# Repeat for dataset 5
start_time <- Sys.time()
result.pca.5 <- analyseBlocks(block.list.5, option = "pca")
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
result.rcpca.5 <- analyseBlocks(block.list.5, option = "rcpca", ncomp = 3)
end_time <- Sys.time()
end_time - start_time

layout(matrix(c(1,2), ncol = 2))
scoresPlot(result.pca.5, comp = c(1, 2), plabels = 1:block.list.5@n[1])
scoresPlot(result.rcpca.5, comp = c(1, 2), plabels = 1:block.list.5@n[1], consensus.only = TRUE)


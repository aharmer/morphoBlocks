## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 100),
  tidy = TRUE
)
library(kableExtra)
library(here)

## ----warnings = FALSE, message = FALSE----------------------------------------
library(morphoBlocks)

# Generate data blocks for the n = 10 samples in this example
set.seed(1)
block1 <- dodecBlock(n = 10, vertex_shift = c(10, 6, 1, 1, 1, 1, 1, 1, 1, 1))
block2 <- dodecBlock(n = 10, vertex_shift = c(1, 1.25, 1, 1, 1, 1, 1.25, 1.25, 1.25, 1.25))
block3 <- dodecBlock(n = 10, vertex_shift = c(1.25, 1, 1.25, 1.25, 1.25, 1.25, 1, 1, 1, 1))
block4 <- dodecBlock(n = 10, vertex_shift = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
block5 <- dodecBlock(n = 10, vertex_shift = c(25, 25, 25, 25, 25, 25, 25, 25, 25, 25))
block6 <- dodecBlock(n = 10, vertex_shift = c(25, 25, 25, 25, 25, 25, 25, 25, 25, 25))
block7 <- dodecBlock(n = 10, vertex_shift = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
block8 <- dodecBlock(n = 10, vertex_shift = c(-1.25, 1, -1.25, -1.25, -1.25, -1.25, 1, 1, 1, 1))
block9 <- dodecBlock(n = 10, vertex_shift = c(1, -1.25, 1, 1, 1, 1, -1.25, -1.25, -1.25, -1.25))
block10 <- dodecBlock(n = 10, vertex_shift = c(-8, -4, 1, 1, 1, 1, 1, 1, 1, 1))

## ----echo = FALSE-------------------------------------------------------------
table_1 <- read.csv(system.file("extdata", "z_vertex_shift.csv", package = "morphoBlocks"), colClasses = "numeric")
table_1$Sample <- as.factor(table_1$Sample)

knitr::kable(table_1, caption = "Table 1. Values added to the *z* coordinate of the first vertex in the corresponding dodecahedron for each subset within each sample of the model dataset.", col.names = gsub("[.]", " ", names(table_1)), align = "c", format.args = list(nsmall = 2)) %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::row_spec(0, bold = TRUE)

## -----------------------------------------------------------------------------
# Combine data blocks into a single block list
block.list <- combineBlocks(blocks = c(
  block1, block2, block3, block4, block5, block6,
  block7, block8, block9, block10
))

## ----echo = FALSE-------------------------------------------------------------
# Generate correlation data for Table 2
table_2 <- matrix(ncol = 10, nrow = 10)
colnames(table_2) <- paste("Samp", seq(1:10))
rownames(table_2) <- paste("Samp", seq(1:10))
for (i in 1:10) {
  sample_a <- as.vector(block.list@block.list[[11]][i, ])
  for (j in 1:10) {
    sample_b <- as.vector(block.list@block.list[[11]][j, ])
    table_2[i, j] <- cor(sample_a, sample_b)
  }
}

table_2 <- read.csv(system.file("extdata", "cor_matrix.csv", package = "morphoBlocks"))
names(table_2)[1] <- ""

knitr::kable(table_2, caption = "Table 2. Correlation between multi-part samples after the data blocks in each sample were Procrustes-transformed and scaled with the normalised centroid size method described by Collyer et al. (2020).", col.names = gsub("[.]", " ", colnames(table_2)), align = "c", digits = 3) %>%
  kableExtra::kable_classic(full_width = TRUE) %>%
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::column_spec(1, bold = TRUE) %>%
  kableExtra::column_spec(2, color = "black", background = ifelse(table_2[, 2] == 1, "#A0DFAF", ifelse(table_2[, 2] == 0.998, "#99DBBB", ifelse(table_2[, 2] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(3, color = "black", background = ifelse(table_2[, 3] == 1, "#A0DFAF", ifelse(table_2[, 3] == 0.998, "#99DBBB", ifelse(table_2[, 3] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(4, color = "black", background = ifelse(table_2[, 4] == 1, "#A0DFAF", ifelse(table_2[, 4] == 0.998, "#99DBBB", ifelse(table_2[, 4] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(5, color = "black", background = ifelse(table_2[, 5] == 1, "#A0DFAF", ifelse(table_2[, 5] == 0.998, "#99DBBB", ifelse(table_2[, 5] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(6, color = "black", background = ifelse(table_2[, 6] == 1, "#A0DFAF", ifelse(table_2[, 6] == 0.998, "#99DBBB", ifelse(table_2[, 6] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(7, color = "black", background = ifelse(table_2[, 7] == 1, "#A0DFAF", ifelse(table_2[, 7] == 0.998, "#99DBBB", ifelse(table_2[, 7] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(8, color = "black", background = ifelse(table_2[, 8] == 1, "#A0DFAF", ifelse(table_2[, 8] == 0.998, "#99DBBB", ifelse(table_2[, 8] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(9, color = "black", background = ifelse(table_2[, 9] == 1, "#A0DFAF", ifelse(table_2[, 9] == 0.998, "#99DBBB", ifelse(table_2[, 9] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(10, color = "black", background = ifelse(table_2[, 10] == 1, "#A0DFAF", ifelse(table_2[, 10] == 0.998, "#99DBBB", ifelse(table_2[, 10] >= 0.985, "#8FC7C6", "#A17FA9")))) %>%
  kableExtra::column_spec(11, color = "black", background = ifelse(table_2[, 11] == 1, "#A0DFAF", ifelse(table_2[, 11] == 0.998, "#99DBBB", ifelse(table_2[, 11] >= 0.985, "#8FC7C6", "#A17FA9"))))

## -----------------------------------------------------------------------------
# Analyse block list with principal component analysis
result.pca <- analyseBlocks(block.list, option = "pca")

## -----------------------------------------------------------------------------
# Generate explained variance values
(round(result.pca$result$sdev^2 / (sum(result.pca$result$sdev^2)), 3))

## ----eval = FALSE-------------------------------------------------------------
#  # Figure 1
#  loadingsPlot(result.pca, comp = 1, cex.3d = 10)
#  
#  # Figure 2
#  loadingsPlot(result.pca, comp = 2, cex.3d = 10)

## ----echo = FALSE, loadingsPlot01, fig.align = "center", out.width = "90%", fig.cap = "Figure 1. Component one loadings from a principal component analysis of simulated shape data. Beginning top left and reading right and down, each group of shows mean position of landmarks in each data block. Orange represents points that have larger loadings and blue represents points that have smaller loadings."----
knitr::include_graphics(here::here("vignettes", "loadingsPlotExample02.png"))

## ----echo = FALSE, loadingsPlot02, fig.align = "center", out.width = "90%", fig.cap = "Figure 2. Component two loadings from a principal component analysis of simulated shape data. See Figure 1 for key."----
knitr::include_graphics(here::here("vignettes", "loadingsPlotExample03.png"))

## ----fig.align = "center", out.width = "90%", fig.height = 4, fig.cap = "Figure 3. Principal component one (PC1) and PC2 score values (left), and PC2 and PC3 score values (right), for the dataset comprised of ten blocks that were scaled, transformed and combined. Samples numbers shown with labels."----
# Figure 3
par(mfrow = c(1, 2))
scoresPlot(result.pca, comp = c(1, 2), plabels = 1:10)
scoresPlot(result.pca, comp = c(2, 3), plabels = 1:10)

## -----------------------------------------------------------------------------
# Analyse block list
result.rcpca <- analyseBlocks(block.list, option = "rcpca", ncomp = 20)

## -----------------------------------------------------------------------------
# Average variance explained for the consensus space of the RCPCA
result.rcpca$result$AVE[[1]][[11]]

## ----eval = FALSE-------------------------------------------------------------
#  # Figure 4
#  loadingsPlot(result.rcpca, comp = 1, cex.3d = 10)
#  
#  # Figure 5
#  loadingsPlot(result.rcpca, comp = 2, cex.3d = 10)

## ----echo = FALSE, loadingsPlot03, fig.align = "center", out.width = "90%", fig.cap = "Figure 4. Global component one loadings from a regularised consensus principal component analysis of simulated shape data. See Figure 1 for key."----
knitr::include_graphics(here::here("vignettes", "loadingsPlotExample04.png"))

## ----echo = FALSE, loadingsPlot04, fig.align = "center", out.width = "90%", fig.cap = "Figure 5. Global component two loadings from a regularised consensus principal component analysis of simulated shape data. See Figure 1 for key."----
knitr::include_graphics(here::here("vignettes", "loadingsPlotExample05.png"))

## ----fig.align = "center", out.width = "90%", fig.height = 4, fig.cap = "Figure 6. Global component one (GC1) and two (GC2) score values (left), and GC2 and GC3 score values (right), for the dataset comprised of ten blocks that were scaled, transformed and combined. Samples numbers shown with labels."----
# Figure 6
par(mfrow = c(1, 2))
scoresPlot(result.rcpca, comp = c(1, 2), plabels = 1:10, consensus.only = TRUE)
scoresPlot(result.rcpca, comp = c(2, 3), plabels = 1:10, consensus.only = TRUE)


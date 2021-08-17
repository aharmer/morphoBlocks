## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 90),
  tidy = TRUE
)

## ----warnings = FALSE, message = FALSE----------------------------------------
library(morphoBlocks)
data(penguinWings)

## -----------------------------------------------------------------------------
# Extract and average the landmark configurations
hum_av = (hum1@raw + hum2@raw + hum3@raw)/3
rad_av = (rad1@raw + rad2@raw + rad3@raw)/3
uln_av = (uln1@raw + uln2@raw + uln3@raw)/3

## -----------------------------------------------------------------------------
# Format the averaged landmark configurations into data blocks
block1 = formatBlock(hum_av, curves = hum1@curves, k = 3, gpa = TRUE)
block2 = formatBlock(rad_av, curves = rad1@curves, k = 3, gpa = TRUE)
block3 = formatBlock(uln_av, curves = uln1@curves, k = 3, gpa = TRUE)

## -----------------------------------------------------------------------------
# Scale and combine data blocks into a single list of blocks
blocklist = combineBlocks(blocks = c(block1, block2, block3))

## -----------------------------------------------------------------------------
# Analyse the list of data blocks using RCPCA
result = analyseBlocks(blocklist, ncomp = 10)

## -----------------------------------------------------------------------------
# Setup colour vector to show different ages of fossil penguins. Paleocene (brown), stem-lineage penguins from the Oligocene (light brown), and extant penguins (white).
pcol = c("#ffffff", "#ffffff", "#ffffff", "#ffffff", "#e6b481", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#feebd3", "#feebd3", "#ffffff", "#e6b481", "#feebd3", "#ffffff")

# Plot consensus space showing global component one (GC1) and global component one (GC2) 
scoresPlot(result, pcol = pcol)


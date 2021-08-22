## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 90),
  tidy = TRUE
)
library(kableExtra)

## ----warnings = FALSE, message = FALSE----------------------------------------
library(morphoBlocks)
data(penguinWings)

## -----------------------------------------------------------------------------
# Extract and average the landmark configurations

hum_av <- (hum1@raw + hum2@raw + hum3@raw) / 3
rad_av <- (rad1@raw + rad2@raw + rad3@raw) / 3
uln_av <- (uln1@raw + uln2@raw + uln3@raw) / 3

## -----------------------------------------------------------------------------
# Format the averaged landmark configurations into data blocks

block1 <- formatBlock(hum_av, curves = hum1@curves, k = 3, gpa = TRUE)
block2 <- formatBlock(rad_av, curves = rad1@curves, k = 3, gpa = TRUE)
block3 <- formatBlock(uln_av, curves = uln1@curves, k = 3, gpa = TRUE)

## -----------------------------------------------------------------------------
# Scale and combine data blocks into a single list of blocks

blocklist <- combineBlocks(blocks = c(block1, block2, block3))

## -----------------------------------------------------------------------------
# Analyse the list of data blocks using RCPCA

result <- analyseBlocks(blocklist, ncomp = 10)

## ----fig.align = "center"-----------------------------------------------------
# Setup colour vector to show different ages of fossil penguins. Paleocene (brown), stem-lineage penguins from the Oligocene (light brown), and extant penguins (white).

pcol <- c("#ffffff", "#ffffff", "#ffffff", "#ffffff", "#e6b481", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#feebd3", "#feebd3", "#ffffff", "#e6b481", "#feebd3", "#ffffff")


# Plot consensus space showing global component one (GC1) and global component one (GC2)

scoresPlot(result, pcol = pcol)

## ----eval = FALSE-------------------------------------------------------------
#  # Plot loadings for global component one (GC1)
#  
#  loadingsPlot(result, cex.3d = 15)

## ----echo = FALSE-------------------------------------------------------------
dat <- read.csv(system.file("extdata", "penguin_bone_metadata.csv", package = "morphoBlocks"))

knitr::kable(dat, caption = "Table 1. 3D digital replicas produced from bones of modern and fossil penguins.", align = "c") %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::column_spec(2:3, italic = TRUE) %>%
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::footnote("Institution abbreviations: CM, Canterbury Museum, Christchurch, New Zealand; NMNZ, Museum of New Zealand Te Papa Tongarewa, Wellington, New Zealand; OM, Otago Museum, Dunedin, New Zealand, New Zealand; UC, University of Canterbury, Christchurch, New Zealand (specimen held at OU, Geology Museum, University of Otago, Dunedin, New Zealand). Specimen ages from or compiled by Slack et al. (2006), Ksepka and Ando (2011), and Ksepka et al. (2012) and literature reviewed therein.")


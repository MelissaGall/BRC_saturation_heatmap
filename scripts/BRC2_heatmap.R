# Heatmap of mean log2(ratio) per amino acid substitution and treatment
# Usage: Rscript scripts/BRC2_heatmap.R [counts_BRC2.xlsx] [BRC2_heatmap.pdf]

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

args    <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "counts_BRC2.xlsx"
outfile <- if (length(args) >= 2) args[2] else "BRC2_heatmap.pdf"

treatments <- c("DMSO", "CIS", "OLA")
ref_aa     <- c("F", "Y", "S", "A")
aa_types   <- list("Positive"      = c("H", "K", "R"),
                   "Negative"      = c("D", "E"),
                   "Polar neutral" = c("T", "N", "Q", "S"),
                   "Non-polar"     = c("A", "V", "C", "L", "I", "M"),
                   "Aromatic"      = c("F", "W", "Y"),
                   "Unique"        = c("G", "P"))

aa_type <- function(aa) {
  out <- rep("ERROR", length(aa))
  for (t in names(aa_types)) out[aa %in% aa_types[[t]]] <- t
  out[aa == "*"] <- "Term"
  out
}

# --- Data ---
data <- read.xlsx(infile)

# Mean of replicates per treatment, then mean per AA x AAref
for (t in treatments) {
  data[[t]] <- rowMeans(data[, paste0(t, "/D3.RP", 1:3, ".ratio")])
}

hm <- data %>%
  mutate(Type = aa_type(AA), AAref = factor(AAref, levels = ref_aa)) %>%
  group_by(AA, AAref, Type) %>%
  summarise(across(all_of(treatments), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_wider(names_from = AAref, values_from = all_of(treatments))

# --- Matrix ---
cols <- paste(rep(treatments, each = length(ref_aa)), ref_aa, sep = "_")
mat  <- log2(as.matrix(hm[, cols]))
rownames(mat) <- hm$AA

# "Silent" row = WT codon (AA == AAref) for each column
silent <- sapply(cols, function(cl) mat[sub(".*_", "", cl), cl])
mat <- rbind(mat, Silent = silent)
mat[is.infinite(mat)] <- NA

# --- Annotations ---
row_type <- c(hm$Type, "")
row_type[rownames(mat) == "*"] <- " "
row_type <- factor(row_type, levels = c(names(aa_types), "", " "))
rownames(mat)[rownames(mat) == "*"] <- "Stop"

col_split <- factor(rep(treatments, each = length(ref_aa)), levels = treatments)
colnames(mat) <- rep(ref_aa, length(treatments))

col_fun <- colorRamp2(c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
                      c("#DC3220", "#ebe9ef", "#005AB5"))

# --- Plot ---
ht <- Heatmap(mat,
              name = "log2(ratio)",
              col = col_fun, na_col = "grey",
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = row_type, column_split = col_split,
              row_names_side = "left", column_names_rot = 0,
              width = unit(ncol(mat) * 0.9, "cm"),
              height = unit(nrow(mat) * 0.9, "cm"),
              # white dot = WT residue
              cell_fun = function(j, i, x, y, w, h, fill) {
                if (rownames(mat)[i] == colnames(mat)[j]) {
                  grid.points(x, y, pch = 21, size = unit(6, "mm"),
                              gp = gpar(fill = "white", col = "black", lwd = 0.5))
                }
              })

pdf(outfile, height = nrow(mat) * 0.4)
draw(ht)
invisible(dev.off())
message("Written: ", outfile)


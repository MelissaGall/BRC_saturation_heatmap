# Libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(tidyr)
library(circlize)

# Version
vers <- 2

# Load data
data <- read.xlsx("counts_BRC2.xlsx")

# Add AA infos
data <- data %>%
  mutate(Consequence = case_when(AA == "*" ~ "Term",
                                 AA == AAref ~ "Silent",
                                 TRUE ~ "Missense")) %>%
  mutate(Type = case_when(AA %in% c("H", "K", "R") ~ "Positive",
                          AA %in% c("D", "E") ~ "Negative",
                          AA %in% c("A", "V", "C", "L", "I", "M") ~ "Non-polar",
                          AA %in% c("T", "N", "Q", "S") ~ "Polar neutral",
                          AA %in% c("F", "W", "Y") ~ "Aromatic",
                          AA %in% c("G", "P") ~ "Unique",
                          AA == "*" ~ "Term",
                          TRUE ~ "ERROR"))

# Prepare data 
heatmap_data <- data.frame("AA" = data$AA,
                           "AAref" = data$AAref,
                           "Consequence" = data$Consequence,
                           "Type" = data$Type,
                           "DMSO" = rowMeans(data[, c("DMSO/D3.RP1.ratio",
                                                      "DMSO/D3.RP2.ratio",
                                                      "DMSO/D3.RP3.ratio")]),
                           "CIS" = rowMeans(data[, c("CIS/D3.RP1.ratio",
                                                     "CIS/D3.RP2.ratio",
                                                     "CIS/D3.RP3.ratio")]),
                           "OLA" = rowMeans(data[, c("OLA/D3.RP1.ratio",
                                                     "OLA/D3.RP2.ratio",
                                                     "OLA/D3.RP3.ratio")]))

heatmap_data$AAref <- factor(heatmap_data$AAref, levels = c("F", "Y", "S", "A"))

# Get mean per compound
heatmap_data <- heatmap_data %>%
  group_by(AA, AAref, Type) %>%
  summarize(DMSO = mean(DMSO, na.rm = TRUE),
            CIS = mean(CIS, na.rm = TRUE),
            OLA = mean(OLA, na.rm = TRUE))

heatmap_data <- pivot_wider(heatmap_data, names_from = AAref, values_from = c(DMSO, CIS, OLA))

# Convert as matrix
mat <- as.matrix(heatmap_data[, -c(1,2)])
mat <- log2(mat)
rownames(mat) <- c(heatmap_data$AA)
mat <- rbind(mat, Silent = NA)

mat["Silent",] <- c(mat["F", "DMSO_F"],
                    mat["Y", "DMSO_Y"],
                    mat["S", "DMSO_S"],
                    mat["A", "DMSO_A"],
                    mat["F", "CIS_F"],
                    mat["Y", "CIS_Y"],
                    mat["S", "CIS_S"],
                    mat["A", "CIS_A"],  
                    mat["F", "OLA_F"],
                    mat["Y", "OLA_Y"],
                    mat["S", "OLA_S"],
                    mat["A", "OLA_A"])

mat[mat == -Inf] <- NA

# Color scale
col_fun <- colorRamp2(
  c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), 
  c("#DC3220", "#ebe9ef", "#005AB5"),
)

# Row annotation
annot <- data.frame("AA" = rownames(mat))
annot <- merge(annot, heatmap_data[, c("AA", "Type")], all.x = T, sort = F)
annot$Type[annot$AA == "Silent"] <- ""
annot$Type[annot$AA == "*"] <- " "
annot$AA[annot$AA == "*"] <- "Stop"

annot$Type <- factor(annot$Type, levels = c("Positive", "Negative", "Polar neutral", "Non-polar", "Aromatic", "", " "))

colnames(mat) <- c("F", "Y", "S", "A", "F", "Y", "S", "A", "F", "Y", "S", "A")

# Column annotation
col_split <- c(rep("DMSO", 4), rep("CIS", 4), rep("OLA", 4))
col_split <- factor(col_split, levels = unique(col_split))

# Dot matrix
dot_matrix <- outer(
  rownames(mat),
  colnames(mat),
  FUN = "=="
)

# Save
pdf(paste0("BRC2_heatmap_", vers, ".pdf"), height = nrow(mat) * 0.4)
Heatmap(mat,
        col = col_fun, 
        na_col = "grey", 
        cluster_columns = F, 
        split = annot$Type, 
        column_split = col_split,
        name = "log2(ratio)",
        width = unit(ncol(mat) * 0.9, "cm"),
        height = unit(nrow(mat) * 0.9, "cm"),
        row_names_side = "left",
        column_names_rot = 0, cluster_rows = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (rownames(mat)[i] == colnames(mat)[j]) {
            grid.points(
              x, y,
              pch = 21,
              gp = gpar(fill = "white", col = "black", lwd = 0.5),
              size = unit(6, "mm")
            )
          }  
        })
dev.off()

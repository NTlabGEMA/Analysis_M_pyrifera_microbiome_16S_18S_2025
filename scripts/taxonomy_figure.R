# ===============================
# Analysis of Macrocystis pyrifera 16S microbiome
# Script cleaned for GitHub upload
# ===============================

# Load required packages
library(readxl)
library(stringr)
library(tibble)
library(phyloseq)
library(ggplot2)

# Set working directory (relative paths recommended)
# setwd("your_project_folder") # ajusta según tu repo

# -------------------------------
# Load data
# -------------------------------
otu_mat <- read_excel("asv_count_16s.xlsx") %>%
  tibble::column_to_rownames("ASVNumber")

tax_mat <- read_excel("asv_taxonomy_16s.xlsx") %>%
  tibble::column_to_rownames("ASVNumber")

samples_df <- read_excel("sample_inf_16s.xlsx") %>%
  tibble::column_to_rownames("Samples")

# Convert to matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Create phyloseq object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)
ps <- phyloseq(OTU, TAX, samples)

# Save phyloseq object
saveRDS(ps, file = "phyloseq_16s.rds")

# -------------------------------
# Transform relative abundances
# -------------------------------
ps.rel2 <- transform_sample_counts(ps, function(x) x/sum(x)*100)

# Aggregate at Class level
physeq.glom <- tax_glom(ps.rel2, taxrank = "Class", NArm = TRUE)
physeq.melt <- psmelt(physeq.glom)
physeq.melt$Class <- as.character(physeq.melt$Class)
physeq.melt$Class[physeq.melt$Abundance < 0.1] <- "<10% Abundance"

# Sort classes by abundance
physeq.melt$Class <- factor(physeq.melt$Class)
physeq.melt$Class <- reorder(physeq.melt$Class, physeq.melt$Abundance)

# -------------------------------
# Load colors
# -------------------------------
class_colors <- read_excel("colors_class.xlsx", sheet = "1%")
class_colors <- structure(class_colors$color_name, .Names = class_colors$Class)

# Order samples
physeq.melt$Sample <- factor(physeq.melt$Sample, levels = c(
  "San Antonio (IS)", "Ilque (IS)", "Los Chonos (SP)", "Pargua (SP)",
  "Topocalma (C)", "Navidad (C)", "Algarrobo (N)", "Las Docas (N)"
))

# Order classes according to colors
physeq.melt$Class <- factor(physeq.melt$Class, levels = names(class_colors))

# -------------------------------
# Taxonomy chart at class level
# -------------------------------
bar.plot.Class <- ggplot(data = physeq.melt, 
                         aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = class_colors) +
  ylab("Relative abundance (%)") +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = NULL, fill = "Class") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 30),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 15),
    axis.title.x = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 30, face = "bold", angle = 0, vjust = 0.2),
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5, vjust = 0.2)
  ) +
  coord_flip()

# Save plot
ggsave("class_16S1_2025.png", plot = bar.plot.Class, units = "in", width = 20, height = 15)

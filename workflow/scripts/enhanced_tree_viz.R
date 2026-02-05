#!/usr/bin/env Rscript
# Vibrion Sentinel: Enhanced Phylogenetic Tree Visualization
# Generates publication-quality trees similar to Haiti cholera network phylogeny
# Uses ggtree for aesthetics and network-style layout

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
  library(treeio)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: enhanced_tree_viz.R <tree.nwk> <metadata.json> <output.png> [width] [height]\n")
  quit(status = 1)
}

tree_path <- args[1]
metadata_path <- args[2]
output_path <- args[3]
fig_width <- if (length(args) >= 4) as.numeric(args[4]) else 14
fig_height <- if (length(args) >= 5) as.numeric(args[5]) else 10

cat("Enhanced Tree Visualization\n")
cat("===========================\n")
cat(sprintf("Tree: %s\n", tree_path))
cat(sprintf("Metadata: %s\n", metadata_path))
cat(sprintf("Output: %s\n", output_path))

# Load tree
tryCatch({
  tree <- read.tree(tree_path)
  cat(sprintf("✓ Loaded tree with %d tips\n", length(tree$tip.label)))
}, error = function(e) {
  cat(sprintf("✗ Error loading tree: %s\n", e$message))
  quit(status = 1)
})

# Load metadata
metadata <- NULL
if (file.exists(metadata_path)) {
  tryCatch({
    metadata <- fromJSON(metadata_path)
    cat(sprintf("✓ Loaded metadata for %d strains\n", length(metadata)))
  }, error = function(e) {
    cat(sprintf("⚠ Warning: Could not load metadata: %s\n", e$message))
  })
}

# Create metadata data frame
if (!is.null(metadata)) {
  meta_df <- do.call(rbind, lapply(names(metadata), function(name) {
    data.frame(
      label = name,
      type = ifelse(is.null(metadata[[name]]$type), "unknown", metadata[[name]]$type),
      event = ifelse(is.null(metadata[[name]]$event), "", metadata[[name]]$event),
      year = ifelse(is.null(metadata[[name]]$year), NA, metadata[[name]]$year),
      stringsAsFactors = FALSE
    )
  }))
} else {
  # Default metadata based on naming patterns
  meta_df <- data.frame(
    label = tree$tip.label,
    type = sapply(tree$tip.label, function(x) {
      if (grepl("2010|2011|2012_clin|2014", x)) return("clinical")
      if (grepl("2012_Env|2015|2020", x)) return("environmental")
      if (grepl("2013|Adaptation", x)) return("mixed")
      return("unknown")
    }),
    event = sapply(tree$tip.label, function(x) {
      if (grepl("2010", x)) return("Initiale")
      if (grepl("2012_clin", x)) return("Clin")
      if (grepl("2012_Env", x)) return("Eau")
      if (grepl("2013|2014", x)) return("Adaptation")
      if (grepl("2014.*Diversification", x)) return("Diversification")
      if (grepl("2015", x)) return("Persistante")
      if (grepl("2020", x)) return("Ancêtre 2022")
      if (grepl("2022", x)) return("Résurgence")
      return("")
    }),
    year = as.numeric(gsub(".*?(\\d{4}).*", "\\1", tree$tip.label)),
    stringsAsFactors = FALSE
  )
  cat("⚠ Using auto-generated metadata from tip labels\n")
}

# Color scheme matching Haiti phylogeny image
color_map <- c(
  "clinical" = "#e74c3c",      # Red
  "environmental" = "#2ecc71",  # Green
  "mixed" = "#f39c12",         # Yellow/Orange
  "unknown" = "#95a5a6"        # Gray
)

# Node sizes based on events
meta_df$node_size <- ifelse(meta_df$event != "", 6, 3)
meta_df$color <- color_map[meta_df$type]

# Create basic tree plot
p <- ggtree(tree, layout = "daylight", branch.length = "none") %<+% meta_df

# Add styled nodes
p <- p + 
  geom_tippoint(aes(color = type, size = node_size), alpha = 0.8) +
  scale_color_manual(
    values = color_map,
    name = "Type de Souche",
    labels = c("Clinical", "Environmental", "Mixed", "Unknown")
  ) +
  scale_size_continuous(range = c(3, 8), guide = "none")

# Add labels with events
p <- p + 
  geom_tiplab(
    aes(label = ifelse(event != "", paste0(label, "\n(", event, ")"), label)),
    size = 3.5,
    fontface = "bold",
    offset = 0.1
  )

# Add title
p <- p + 
  ggtitle("Phylogenetic Network of V. cholerae Strains") +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Save high-resolution PNG
tryCatch({
  ggsave(
    output_path,
    plot = p,
    width = fig_width,
    height = fig_height,
    dpi = 300,
    bg = "white"
  )
  cat(sprintf("✓ Saved enhanced tree to %s\n", output_path))
}, error = function(e) {
  cat(sprintf("✗ Error saving plot: %s\n", e$message))
  quit(status = 1)
})

# Also create a circular layout version
output_circular <- gsub("\\.png$", "_circular.png", output_path)
p_circular <- ggtree(tree, layout = "circular", branch.length = "none") %<+% meta_df +
  geom_tippoint(aes(color = type, size = node_size), alpha = 0.8) +
  scale_color_manual(values = color_map, name = "Type de Souche") +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  geom_tiplab(
    aes(label = ifelse(event != "", paste0(label, "\n(", event, ")"), label)),
    size = 3,
    fontface = "bold",
    offset = 0.15
  ) +
  ggtitle("Phylogenetic Network (Circular Layout)") +
  theme_tree() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(
  output_circular,
  plot = p_circular,
  width = fig_width,
  height = fig_height,
  dpi = 300,
  bg = "white"
)
cat(sprintf("✓ Saved circular layout to %s\n", output_circular))

cat("\n✅ Enhanced tree visualization complete!\n")

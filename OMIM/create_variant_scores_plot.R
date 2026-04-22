#################################### Configuration ########################################################################

# Base paths - centralized configuration
BASE_PATH <- "/home/joern/Aktuell/FP7NGS_RareVariants"
PATH_ORIGINALE <- file.path(BASE_PATH, "09Originale")
PATH_R <- file.path(BASE_PATH, "08AnalyseProgramme", "R")

# File paths
FILE_SCORES <- file.path(PATH_ORIGINALE, "RareVars_Scores_Tab.xlsx")

#################################### Libraries ########################################################################

# Suppress package startup messages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(viridisLite)
  library(ggnewscale)
  library(readxl)
})

#################################### Helper Functions ########################################################################

#' Check if file exists and stop with informative message if not
check_file_exists <- function(filepath, description = "File") {
  if (!file.exists(filepath)) {
    stop(sprintf("%s not found at: %s", description, filepath))
  }
  invisible(TRUE)
}

######################### Load Data ######################################################

message("Loading variant scores data...")
check_file_exists(FILE_SCORES, "Scores Excel file")

RareVars_Scores_Tab <- readxl::read_excel(FILE_SCORES)

######################### Data Cleaning and Preparation ######################################################

# Remove rows with missing scores
RareVars_Scores_Tab <- RareVars_Scores_Tab %>%
  dplyr::filter(!is.na(Score))

message(sprintf("  - Total scores with data: %d", nrow(RareVars_Scores_Tab)))

# Define factor levels for consistent ordering
RareVars_Scores_Tab <- RareVars_Scores_Tab %>%
  dplyr::mutate(
    Reference = stringr::str_replace_all(Reference, "[()]", "")
  ) %>%
  dplyr::filter(!is.na(Score)) %>%
  dplyr::mutate(
    Category = factor(
      Category,
      levels = unique(na.omit(RareVars_Scores_Tab$Category))
    ),
    `Model type` = factor(
      `Model type`,
      levels = unique(na.omit(RareVars_Scores_Tab$`Model type`))
    ),
    `Variant types` = factor(
      `Variant types`,
      levels = unique(na.omit(RareVars_Scores_Tab$`Variant types`))
    ),
    `Genomic scope` = factor(
      `Genomic scope`,
      levels = unique(na.omit(RareVars_Scores_Tab$`Genomic scope`))
    ),
    Score = factor(Score)
  ) %>%
  dplyr::arrange(Category, `Model type`, `Variant types`, `Genomic scope`, Score) %>%
  dplyr::mutate(Score_id = row_number())

######################### Prepare Label Data ######################################################

# Create outer label data with Score and Reference
label_data <- RareVars_Scores_Tab %>%
  dplyr::distinct(Score, Reference, Score_id) %>%
  dplyr::arrange(Score_id)

n_bar <- nrow(label_data)
message(sprintf("  - Unique scores for plotting: %d", n_bar))

# Calculate label angles for circular plot
label_data <- label_data %>%
  dplyr::mutate(
    angle_raw = 90 - 360 * (Score_id - 0.5) / n_bar,
    hjust = ifelse(angle_raw < -90, 1, 0),
    angle = ifelse(angle_raw < -90, angle_raw + 180, angle_raw)
  )

######################### Create Circular Plot ######################################################

message("Creating circular variant scores plot...")

p_variant_scores <- ggplot() +
  # Genomic scope ring (innermost)
  geom_rect(
    data = RareVars_Scores_Tab,
    aes(
      xmin = Score_id - 0.45,
      xmax = Score_id + 0.45,
      ymin = 0.8,
      ymax = 2.0,
      fill = `Genomic scope`
    ),
    color = "white", size = 0.3
  ) +
  scale_fill_viridis_d(option = "magma", name = "Genomic scope") +
  ggnewscale::new_scale_fill() +

  # Variant types ring
  geom_rect(
    data = RareVars_Scores_Tab,
    aes(
      xmin = Score_id - 0.45,
      xmax = Score_id + 0.45,
      ymin = 2.3,
      ymax = 3.5,
      fill = `Variant types`
    ),
    color = "white", size = 0.3
  ) +
  scale_fill_viridis_d(option = "inferno", name = "Variant types") +
  ggnewscale::new_scale_fill() +

  # Model type ring
  geom_rect(
    data = RareVars_Scores_Tab,
    aes(
      xmin = Score_id - 0.45,
      xmax = Score_id + 0.45,
      ymin = 3.8,
      ymax = 5.0,
      fill = `Model type`
    ),
    color = "white", size = 0.3
  ) +
  scale_fill_viridis_d(option = "viridis", name = "Model type") +
  ggnewscale::new_scale_fill() +

  # Category ring (outermost)
  geom_rect(
    data = RareVars_Scores_Tab,
    aes(
      xmin = Score_id - 0.45,
      xmax = Score_id + 0.45,
      ymin = 5.3,
      ymax = 6.5,
      fill = Category
    ),
    color = "white", size = 0.3
  ) +
  scale_fill_viridis_d(option = "plasma", name = "Category") +

  # Score labels (outer)
  geom_text(
    data = label_data,
    aes(x = Score_id, y = 7.1, label = Score, angle = angle, hjust = hjust),
    size = 3.2,
    fontface = "bold"
  ) +

  # Reference labels (outermost)
  geom_text(
    data = label_data,
    aes(x = Score_id, y = 11.2, label = Reference, angle = angle, hjust = hjust),
    size = 2.6
  ) +

  # Ring category labels with semi-transparent background
  geom_label(
    data = data.frame(
      x = 1,
      y = c(1.4, 2.9, 4.4, 5.9),
      label = c("Genomic scope", "Variant types", "Model type", "Category")
    ),
    aes(x = x, y = y, label = label),
    size = 3.8,
    fontface = "bold",
    color = "black",
    fill = "white",
    alpha = 0.7,
    label.size = 0,
    label.padding = unit(0.15, "lines")
  ) +

  coord_polar(start = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  guides(
    `Genomic scope` = guide_legend(order = 1, ncol = 1),
    Category = guide_legend(order = 2, ncol = 1),
    `Variant types` = guide_legend(order = 3, ncol = 1),
    `Model type` = guide_legend(order = 4, ncol = 2)
  )

print(p_variant_scores)

######################### Save Plot ######################################################

message("Saving plot...")

# Save as SVG
ggsave(
  filename = file.path(PATH_R, "Variant_Scores_Circular.svg"),
  plot     = p_variant_scores,
  width    = 16,
  height   = 14,
  units    = "in"
)

# Save as PNG
ggsave(
  filename = file.path(PATH_R, "Variant_Scores_Circular.png"),
  plot     = p_variant_scores,
  width    = 16,
  height   = 14,
  units    = "in",
  dpi      = 300
)

message(sprintf("  - Saved plots to: %s", PATH_R))

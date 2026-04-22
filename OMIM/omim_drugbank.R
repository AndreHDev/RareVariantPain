#################################### Configuration ########################################################################

# Base paths - centralized configuration
BASE_PATH <- "/home/joern/Aktuell/FP7NGS_RareVariants"
PATH_DRUGBANK <- file.path(BASE_PATH, "11DrugBank")
PATH_R <- file.path(BASE_PATH, "08AnalyseProgramme", "R")
PATH_OMIM <- file.path(BASE_PATH, "10OMIM")

# File paths
FILE_DRUGBANK_XML <- file.path(PATH_DRUGBANK, "full database.xml")
FILE_OMIM_TSV <- file.path(PATH_OMIM, "OMIM-Entry-Search.tsv")
FILE_WORKSPACE <- file.path(PATH_R, "DrugbankAuslesen.img")

read_DrugBank <- FALSE  # To skip reading if run accidentally

#################################### Libraries ########################################################################
# Suppress package startup messages
suppressPackageStartupMessages({
  library(dbparser)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Only load plotting libraries if needed
# library(ggplot2)
# library(ComplexHeatmap)
# library(viridis)
# library(ggthemes)
# library(ggwordcloud)
# library(ggridges)

#################################### Helper Functions ########################################################################

#' Check if file exists and stop with informative message if not
check_file_exists <- function(filepath, description = "File") {
  if (!file.exists(filepath)) {
    stop(sprintf("%s not found at: %s", description, filepath))
  }
  invisible(TRUE)
}

#' Report join statistics
report_join_stats <- function(before, after, join_type = "join") {
  message(sprintf("%s: %d rows before, %d rows after",
                  join_type, nrow(before), nrow(after)))
}

######################### Load DrugBank Data ######################################################
if (read_DrugBank) {
  # Reference: https://cran.r-project.org/web/packages/dbparser/vignettes/dbparser.html
  
  # Validate file existence before parsing
  check_file_exists(FILE_DRUGBANK_XML, "DrugBank XML")
  
  # Parse XML database (this loads data into dbparser's internal state)
  message("Parsing DrugBank XML database...")
  read_drugbank_xml_db(FILE_DRUGBANK_XML)
  
  # Extract all relevant DrugBank tables
  message("Loading DrugBank tables...")
  drugs <- drugs()
  drug_groups <- drug_groups()
  
  # Drug targets and related tables
  drug_targets <- targets()
  message(sprintf("  - Unique drug targets: %d", length(unique(drug_targets$id))))
  
  drug_targets_actions <- targets_actions()
  drug_targets_polypeptides <- targets_polypeptides()
  message(sprintf("  - Unique target polypeptides (parent_id): %d",
                  length(unique(drug_targets_polypeptides$parent_id))))
  message(sprintf("  - Unique target polypeptides (id): %d",
                  length(unique(drug_targets_polypeptides$id))))
  
  # Additional target annotations (loaded but not actively used)
  drug_targets_polypep_ex_ident <- targets_polypep_ex_ident()
  drug_targets_polypeptides_go <- targets_polypeptides_go()
  drug_targets_polypeptides_pfams <- targets_polypeptides_pfams()
  drug_targets_polypeptides_syn <- targets_polypeptides_syn()
  
  # Enzymes and transporters (loaded but not actively used)
  drug_enzymes <- enzymes()
  drug_enzymes_actions <- enzymes_actions()
  drug_transporters <- transporters()
  drug_transporters_actions <- transporters_actions()
  
  # Optional: Save/load workspace for faster subsequent runs
  # save.image(file = FILE_WORKSPACE)
}
# load(file = FILE_WORKSPACE)

######################### Load and Process OMIM Data ######################################################

# Load OMIM pain-related entries
check_file_exists(FILE_OMIM_TSV, "OMIM TSV")
message("Loading OMIM data...")

omim_raw <- read.csv(
  FILE_OMIM_TSV,
  sep = "\t",
  skip = 4,
  stringsAsFactors = FALSE
)

message(sprintf("  - Total hits: %d", nrow(omim_raw)))

# Identify disease entries (marked with #)
# OMIM prefix meanings: # = disease, * = gene, % = obsolete, + = gene with associated disease
omim_diseases <- omim_raw %>%
  dplyr::mutate(
    is_disease = grepl("^#", MIM.Number),
    OMIM_clean = str_replace_all(MIM.Number, "^[#*%+]+", "")
  ) %>%
  dplyr::filter(is_disease)

message(sprintf("  - Disease entries identified: %d", nrow(omim_diseases)))

# Map OMIM IDs to ENTREZ gene IDs
message("Mapping OMIM IDs to ENTREZ gene IDs...")
omim_keys <- unique(na.omit(omim_diseases$OMIM_clean))

# Suppress expected 1:many mapping warnings
gene2entrez <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = omim_keys,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "OMIM"
  )
)

# Join gene information to OMIM data
omim_with_genes <- omim_diseases %>%
  dplyr::left_join(gene2entrez, by = c("OMIM_clean" = "OMIM"))

n_mapped <- sum(!is.na(omim_with_genes$ENTREZID))
message(sprintf("  - OMIM entries with ENTREZ mapping: %d/%d (%.1f%%)",
                n_mapped, nrow(omim_diseases),
                100 * n_mapped / nrow(omim_diseases)))

######################### Map Drug Targets to ENTREZ IDs ######################################################

# Add ENTREZ IDs to drug target polypeptides based on gene symbols
message("Mapping drug target genes to ENTREZ IDs...")

drug_targets_with_entrez <- drug_targets_polypeptides %>%
  dplyr::mutate(
    ENTREZ_ID = suppressMessages(
      as.character(
        mapIds(
          org.Hs.eg.db,
          keys = gene_name,
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "first"  # Takes first match when multiple exist
        )
      )
    )
  )

n_drug_targets_mapped <- sum(!is.na(drug_targets_with_entrez$ENTREZ_ID))
message(sprintf("  - Drug targets with ENTREZ mapping: %d/%d (%.1f%%)",
                n_drug_targets_mapped, nrow(drug_targets_polypeptides),
                100 * n_drug_targets_mapped / nrow(drug_targets_polypeptides)))

######################### Filter to OMIM-Relevant Drugs ######################################################

# Get list of ENTREZ IDs from OMIM diseases
omim_entrez_ids <- unique(na.omit(as.character(omim_with_genes$ENTREZID)))
message(sprintf("Filtering drugs targeting OMIM disease genes (%d genes)...",
                length(omim_entrez_ids)))

# Filter drug targets to those matching OMIM genes
drug_targets_in_omim <- drug_targets_with_entrez %>%
  dplyr::filter(ENTREZ_ID %in% omim_entrez_ids)

message(sprintf("  - Drug target polypeptides matching OMIM: %d", nrow(drug_targets_in_omim)))

# Get corresponding drug target records
drug_targets_filtered <- drug_targets %>%
  dplyr::filter(id %in% drug_targets_in_omim$parent_id)

message(sprintf("  - Drug targets matching OMIM: %d", nrow(drug_targets_filtered)))

# Get drug information for these targets
drugs_filtered_initial <- drugs$general_information %>%
  dplyr::filter(primary_key %in% drug_targets_filtered$parent_key)

message(sprintf("  - Initial drugs targeting OMIM genes: %d", nrow(drugs_filtered_initial)))

# Filter out drugs with unwanted groups (withdrawn, vet_approved, nutraceutical, illicit, experimental)
excluded_groups <- c("withdrawn", "vet_approved", "nutraceutical", "illicit", "experimental")

# Get drug IDs that have any of the excluded groups
drugs_to_exclude <- drugs$groups %>%
  dplyr::filter(group %in% excluded_groups) %>%
  dplyr::pull(`drugbank-id`) %>%
  unique()

message(sprintf("  - Drugs to exclude (withdrawn/vet/nutraceutical/illicit): %d",
                length(drugs_to_exclude)))

# Filter out excluded drugs
drugs_filtered <- drugs_filtered_initial %>%
  dplyr::filter(!primary_key %in% drugs_to_exclude)

message(sprintf("  - Filtered drugs targeting OMIM genes: %d", nrow(drugs_filtered)))

######################### Create Integrated OMIM-Drug Dataset ######################################################

message("Building integrated OMIM-Drug dataset...")

# Join OMIM diseases with drug information through ENTREZ IDs
omim_drug_integrated <- omim_with_genes %>%
  # Link to drug target polypeptides via ENTREZ ID
  dplyr::left_join(
    drug_targets_in_omim %>%
      dplyr::select(ENTREZ_ID, parent_id),
    by = c("ENTREZID" = "ENTREZ_ID"),
    relationship = "many-to-many"
  ) %>%
  # Link to drug targets via parent_id
  dplyr::left_join(
    drug_targets_filtered %>%
      dplyr::select(id, parent_key),
    by = c("parent_id" = "id"),
    relationship = "many-to-many"
  ) %>%
  # Link to drug information via parent_key
  dplyr::left_join(
    drugs_filtered %>%
      dplyr::select(primary_key, drug_name = name, drug_description = description),
    by = c("parent_key" = "primary_key"),
    relationship = "many-to-many"
  ) %>%
  # Add ATC codes and classifications
  dplyr::left_join(
    drugs$atc_codes %>%
      dplyr::select(drugbank_id = `drugbank-id`, atc_code, atc_level_1 = level_1,
                    atc_level_2 = level_2, atc_level_3 = level_3),
    by = c("parent_key" = "drugbank_id"),
    relationship = "many-to-many"
  ) %>%
  # Clean up intermediate columns
  dplyr::select(-parent_id) %>%
  # Rename for clarity
  dplyr::rename(drug_id = parent_key)

# Report dataset size
n_diseases_with_drugs <- omim_drug_integrated %>%
  dplyr::filter(!is.na(drug_id)) %>%
  dplyr::pull(OMIM_clean) %>%
  dplyr::n_distinct()

message(sprintf("  - OMIM diseases with drug matches: %d/%d",
                n_diseases_with_drugs, dplyr::n_distinct(omim_with_genes$OMIM_clean)))
message(sprintf("  - Total OMIM-drug associations: %d", nrow(omim_drug_integrated)))

######################### Exploratory Analysis ######################################################

# Define search patterns for pain/analgesia-related terms
pain_pattern <- "\\b(pain|analge|anaesth|anesth)\\b"

# Search in drug descriptions
pain_in_desc_idx <- grep(pain_pattern, tolower(omim_drug_integrated$drug_description), ignore.case = TRUE)
message(sprintf("  - Drugs with pain/analgesic mentions in description: %d",
                length(pain_in_desc_idx)))

# Search in ATC classifications
atc_combined <- tolower(c(
  omim_drug_integrated$atc_level_1,
  omim_drug_integrated$atc_level_2,
  omim_drug_integrated$atc_level_3
))
pain_in_atc_idx <- grep(pain_pattern, atc_combined, ignore.case = TRUE)
message(sprintf("  - ATC classification entries with pain/analgesic terms: %d",
                length(pain_in_atc_idx)))

# Summary of ATC codes
message("\nUnique ATC classifications found:")
message(sprintf("  - ATC codes: %d", length(unique(na.omit(omim_drug_integrated$atc_code)))))
message(sprintf("  - ATC level 1: %d", length(unique(na.omit(omim_drug_integrated$atc_level_1)))))
message(sprintf("  - ATC level 2: %d", length(unique(na.omit(omim_drug_integrated$atc_level_2)))))
message(sprintf("  - ATC level 3: %d", length(unique(na.omit(omim_drug_integrated$atc_level_3)))))

# Add ATC system classification (top-level category)
ATC_SYSTEM_LOOKUP <- c(
  A = "Alimentary tract and metabolism",
  B = "Blood and blood forming organs",
  C = "Cardiovascular system",
  D = "Dermatologicals",
  G = "Genito-urinary system and sex hormones",
  H = "Systemic hormonal preparations, excluding sex hormones and insulins",
  J = "Antiinfectives for systemic use",
  L = "Antineoplastic and immunomodulating agents",
  M = "Musculo-skeletal system",
  N = "Nervous system",
  P = "Antiparasitic products, insecticides and repellents",
  R = "Respiratory system",
  S = "Sensory organs",
  V = "Various"
)

omim_drug_integrated <- omim_drug_integrated %>%
  dplyr::mutate(
    atc_system = dplyr::case_when(
      !is.na(atc_code) & nchar(atc_code) > 0 ~ ATC_SYSTEM_LOOKUP[substr(atc_code, 1, 1)],
      TRUE ~ "Unclassified/No ATC"
    )
  )

# Summary of drug therapeutic areas
atc_system_summary <- omim_drug_integrated %>%
  dplyr::filter(!is.na(atc_system)) %>%
  dplyr::count(atc_system, sort = TRUE)

# Add whether the pain or anesthesia is mentioned for the drug 
omim_drug_integrated$pain.use <- dplyr::if_else(
  !is.na(omim_drug_integrated$drug_description) & grepl("\\b(pain|analg|anaesth|anesth)\\b",
                                                        omim_drug_integrated$drug_description, ignore.case = TRUE),
  "Pain or analgesia mentioned",
  "Pain or analgesia not mentioned"
)

message("\nTop therapeutic areas (ATC systems):")
print(head(atc_system_summary, 10))

message("\nTop lines of complete dataframe:")
print(names(omim_drug_integrated))
print(head(omim_drug_integrated, 10))

######################### Create Sankey Plot ######################################################

library(dplyr)
library(ggplot2)
library(ggsankey)

# 1. Clean and prepare
sankey_clean <- omim_drug_integrated %>%
  select(MIM.Number, SYMBOL, drug_name, atc_system, pain.use) %>%
  filter(!is.na(drug_name), !is.na(SYMBOL))

sankey_stages <- sankey_clean %>%
  rename(
    OMIM = MIM.Number,
    Gene = SYMBOL,
    Drug = drug_name,
    ATC  = atc_system,
    Pain = pain.use
  )

# 2. Long format (this ONLY returns x, node, next_x, next_node)
df_long <- sankey_stages %>%
  make_long(OMIM, Gene, Drug, ATC, Pain)

str(df_long)
# x, node, next_x, next_node

# 3. Plot: fill by node (works with the columns you actually have)
p_omim_drug_sankey <- ggplot(
  df_long,
  aes(
    x         = x,
    next_x    = next_x,
    node      = node,
    next_node = next_node,
    fill      = factor(node)
  )
) +
  geom_sankey(flow.alpha = 0.5, node.color = "grey30", show.legend = FALSE) +
  geom_sankey_label(
    aes(label = node),
    size  = 2.5,
    color = "black",
    fill  = "white",
    label.padding = unit(0.1, "lines"),  # Minimal padding
    label.size = 0  # Remove border around labels
  ) +
  # scale_fill_viridis_d(option = "A", alpha = 0.95, name = "Node") +
  theme_sankey(base_size = 14) +
  labs(
    title = "OMIM Disease-Gene-Drug-Therapeutic Area Flow",
    x = NULL, y = NULL
  )

print(p_omim_drug_sankey)

# 4. Save
ggsave(
  filename = file.path(PATH_R, "OMIM_Drug_Sankey_ggplot.pdf"),
  plot     = p_omim_drug_sankey,
  width    = 14,
  height   = 10,
  units    = "in"
)

ggsave(
  filename = file.path(PATH_R, "OMIM_Drug_Sankey_ggplot.png"),
  plot     = p_omim_drug_sankey,
  width    = 14,
  height   = 10,
  units    = "in",
  dpi      = 300
)

######################### Alternative: Plotly Sankey ######################################################
# Load required libraries
library(plotly)
library(htmlwidgets)
library(viridisLite)
library(dplyr)

message("\nCreating Plotly Sankey diagram...")

# Prepare data for Plotly
sankey_clean <- omim_drug_integrated %>%
  dplyr::select(MIM.Number, SYMBOL, drug_name, atc_system, pain.use) %>%
  dplyr::filter(!is.na(drug_name), !is.na(SYMBOL))

# OMIM ORDER - establish once upfront
unique_omim <- sort(unique(sankey_clean$MIM.Number))
n_omim <- length(unique_omim)

if (n_omim <= 12) {
  omim_colors <- viridis(max(3, n_omim), option = "C")
} else {
  omim_colors <- viridis(n_omim, option = "C")
}
names(omim_colors) <- unique_omim

# # Generate colors using a color palette
# library(RColorBrewer)
# if (n_omim <= 12) {
#   omim_colors <- brewer.pal(max(3, n_omim), "Set3")
# } else if (n_omim <= 74) {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# } else {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# }
# 
# library(ggthemes)
# if (n_omim <= 8) {
#   omim_colors <- ggthemes::colorblind_pal()(max(3, n_omim))
# } else if (n_omim <= 74) {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# } else {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# }

names(omim_colors) <- unique_omim


# Create OMIM position column FIRST
sankey_clean$omim_pos <- match(sankey_clean$MIM.Number, unique_omim)

# Generate links maintaining omim_pos order - FIXED first() error
links1 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = MIM.Number, target = SYMBOL) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

links2 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = SYMBOL, target = drug_name) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

links3 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = drug_name, target = atc_system) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

# links4 REVERSED by omim_pos per source (compensates Plotly flip)
links4 <- sankey_clean %>%
  dplyr::group_by(source = atc_system, omim_pos, target = pain.use) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(source, desc(omim_pos), target) %>%
  dplyr::ungroup()

# Combine maintaining relative order
all_links <- dplyr::bind_rows(links1, links2, links3, links4)

# Stable node ordering
nodes <- data.frame(
  name = unique(c(as.character(all_links$source), as.character(all_links$target))), 
  stringsAsFactors = FALSE
)

all_links$source_id <- match(all_links$source, nodes$name) - 1
all_links$target_id <- match(all_links$target, nodes$name) - 1

message(sprintf("  - Nodes: %d", nrow(nodes)))
message(sprintf("  - Links: %d", nrow(all_links)))

# Colors
link_colors_rgba <- sapply(1:nrow(all_links), function(i) {
  omim_id <- all_links$omim_origin[i]
  if (!is.na(omim_id) && omim_id %in% names(omim_colors)) {
    rgb_vals <- col2rgb(omim_colors[omim_id])
    sprintf("rgba(%d,%d,%d,0.5)", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  } else {
    "rgba(150,150,150,0.4)"
  }
})

node_origin <- sapply(nodes$name, function(node_name) {
  rows <- which(all_links$source == node_name | all_links$target == node_name)
  if (length(rows) == 0) return(NA_character_)
  tab <- table(all_links$omim_origin[rows])
  names(tab)[which.max(tab)]
})

node_colors <- sapply(node_origin, function(omim_id) {
  if (!is.na(omim_id) && omim_id %in% names(omim_colors)) {
    omim_colors[omim_id]
  } else {
    "#D3D3D3"
  }
})

# Plot
fig_sankey <- plot_ly(
  type = "sankey", orientation = "h", arrangement = "snap",
  node = list(
    label = nodes$name, pad = 15, thickness = 20,
    line = list(color = "black", width = 0.5),
    color = node_colors
  ),
  link = list(
    source = all_links$source_id, target = all_links$target_id,
    value = all_links$value, color = link_colors_rgba,
    line = list(color = "rgba(0,0,0,0.1)", width = 0.5)
  )
) %>%
  layout(
    title = list(text = "OMIM Disease-Gene-Drug-Therapeutic Area Flow", font = list(size = 18)),
    font = list(size = 8, color = "black"),
    xaxis = list(showgrid = FALSE, zeroline = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE),
    height = 1500, width = 1500,
    paper_bgcolor = "white", plot_bgcolor = "white",
    margin = list(b = 100, t = 80, l = 50, r = 50),
    annotations = list(
      list(x = 0.05, y = -0.05, text = "OMIM", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.28, y = -0.05, text = "Gene", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.51, y = -0.05, text = "Drug", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.74, y = -0.05, text = "ATC System", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.95, y = -0.05, text = "Pain/Analgesia", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold"))
    )
  )

print(fig_sankey)


# Save as interactive HTML
html_file <- file.path(PATH_R, "OMIM_Drug_Sankey_plotly.html")
saveWidget(fig_sankey, html_file, selfcontained = TRUE)

library(plotly)
library(reticulate)

# Set environment FIRST, before library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/joern/.Datenplatte/Subversion/.virtualenvs/r-reticulate/bin/python")

# Now load reticulate - it will use the env var
library(reticulate)

# Check - MUST show r-reticulate path
py_config()
py_module_available("plotly")
py_module_available("kaleido")

virtualenv_root()  # shows base dir
virtualenv_list()  # lists all virtualenvs

Sys.setenv(RETICULATE_PYTHON = "/home/joern/.Datenplatte/Subversion/.virtualenvs/r-reticulate/bin/python")

library(plotly)
library(reticulate)

svg_file <- file.path(PATH_R, "OMIM_Drug_Sankey_plotly.svg")

reticulate::py_install("kaleido==0.2.1", pip = TRUE)  # older stable version
plotly::save_image(fig_sankey, svg_file, width = 1500, height = 1500)

message(sprintf("  - Saved interactive HTML: %s", html_file))
message("\nTo export PNG/SVG, oopen it in a browser and export from there")

######################### Alternative: Chord Diagram ######################################################

# Load required library
library(circlize)

message("\nCreating Chord diagram...")

# Prepare data for chord diagram
# Chord diagrams work best with pairwise connections, so we'll create adjacency matrix
# We'll show connections between consecutive stages

# Create adjacency matrices for each transition
chord_data <- omim_drug_integrated %>%
  dplyr::select(MIM.Number, SYMBOL, drug_name, atc_system, pain.use) %>%
  dplyr::filter(!is.na(drug_name), !is.na(SYMBOL))

chord_data$drug_name <- stringr::str_wrap(chord_data$drug_name, width = 25)

# Create a combined connection matrix
# Transition 1: OMIM -> Gene
mat1 <- chord_data %>%
  dplyr::group_by(MIM.Number, SYMBOL) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 2: Gene -> Drug
mat2 <- chord_data %>%
  dplyr::group_by(SYMBOL, drug_name) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 3: Drug -> ATC
mat3 <- chord_data %>%
  dplyr::group_by(drug_name, atc_system) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 4: ATC -> Pain
mat4 <- chord_data %>%
  dplyr::group_by(atc_system, pain.use) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Combine all transitions into one data frame
colnames(mat1) <- c("from", "to", "value")
colnames(mat2) <- c("from", "to", "value")
colnames(mat3) <- c("from", "to", "value")
colnames(mat4) <- c("from", "to", "value")

chord_links <- dplyr::bind_rows(mat1, mat2, mat3, mat4)

message(sprintf("  - Total connections: %d", nrow(chord_links)))
message(sprintf("  - Unique nodes: %d", length(unique(c(chord_links$from, chord_links$to)))))

# Create color scheme for different node types and track categories
all_nodes <- unique(c(chord_links$from, chord_links$to))
node_colors <- rep("#999999", length(all_nodes))
names(node_colors) <- all_nodes

# Create category assignment for each node
node_categories <- rep("Other", length(all_nodes))
names(node_categories) <- all_nodes

# Color by category and assign category labels
omim_nodes <- all_nodes %in% unique(chord_data$MIM.Number)
gene_nodes <- all_nodes %in% unique(chord_data$SYMBOL)
drug_nodes <- all_nodes %in% unique(chord_data$drug_name)
atc_nodes <- all_nodes %in% unique(chord_data$atc_system)
pain_nodes <- all_nodes %in% unique(chord_data$pain.use)

node_colors[omim_nodes] <- "#E41A1C"  # Red for OMIM
node_colors[gene_nodes] <- "#377EB8"   # Blue for genes
node_colors[drug_nodes] <- "#4DAF4A"   # Green for drugs
node_colors[atc_nodes] <- "#FF7F00"    # Orange for ATC
node_colors[pain_nodes] <- "#984EA3"   # Purple for pain

node_categories[omim_nodes] <- "OMIM"
node_categories[gene_nodes] <- "Gene"
node_categories[drug_nodes] <- "Drug"
node_categories[atc_nodes] <- "ATC"
node_categories[pain_nodes] <- "Pain"

# Create a color function for links based on source node
link_colors <- node_colors[chord_links$from]

# Save as PDF
pdf_file <- file.path(PATH_R, "OMIM_Drug_Chord.pdf")
pdf(pdf_file, width = 8, height = 8)

# Create chord diagram
circos.clear()
chordDiagram(
  chord_links,
  grid.col = node_colors,
  col = link_colors,
  transparency = 0,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.lwd = 2
)

# Add labels with category indicators
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    category <- node_categories[sector_name]
    
    # Add category label in parentheses for first few items of each category
    label_text <- paste0(sector_name, " [", category, "]")
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      label_text,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.3
    )
  },
  bg.border = NA
)

# Add title
title("OMIM Disease-Gene-Drug-Therapeutic Area Connections", cex.main = 1)

# Add comprehensive legend
legend(
  "topleft",
  legend = c(
    "OMIM Disease (Red)",
    "Gene Symbol (Blue)",
    "Drug Name (Green)",
    "Therapeutic Area (Orange)",
    "Pain/Analgesia (Purple)"
  ),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"),
  cex = 0.5,
  bty = "n",
  title = "Node Categories"
)

dev.off()
circos.clear()
message(sprintf("  - Saved: %s", pdf_file))

# Save as SVG
svg_file <- file.path(PATH_R, "OMIM_Drug_Chord.svg")
svg(svg_file, width = 8, height = 8)

circos.clear()
chordDiagram(
  chord_links,
  grid.col = node_colors,
  col = link_colors,
  transparency = 0.2,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.lwd = 2
)

circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    category <- node_categories[sector_name]
    
    # Add category label in parentheses for first few items of each category
    label_text <- paste0(sector_name, " [", category, "]")
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      label_text,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.3
    )
  },
  bg.border = NA
)

title("OMIM Disease-Gene-Drug-Therapeutic Area Connections", cex.main = 1)

legend(
  "topleft",
  legend = c(
    "OMIM Disease (Red)",
    "Gene Symbol (Blue)",
    "Drug Name (Green)",
    "Therapeutic Area (Orange)",
    "Pain/Analgesia (Purple)"
  ),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"),
  cex = 0.5,
  bty = "n",
  title = "Node Categories"
)

dev.off()
circos.clear()
message(sprintf("  - Saved: %s", svg_file))

message("\nChord diagram shows all connections in circular layout.")
message("Note: Sequential flow is less clear than in Sankey/parallel sets.")
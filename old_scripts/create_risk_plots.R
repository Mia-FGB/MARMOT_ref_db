setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database")

# Packages 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readr)
library(fs)  # for file path safety
library(stringr)
library(RColorBrewer)
library(viridis)

# Set paths and read in the data -----

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Expecting: Rscript create_risk_plots.R <input_dir>
if (length(args) < 1) {
  stop("Usage: Rscript create_risk_plots.R <input_dir>")
}

# input_dir is where the downloaded outputs from the pipeline have been saved
input_dir <- args[1]

message("Set paths and read in the data, input_dir ", input_dir)

# Alternatively provide the path, when not running on command line 
# input_dir <- "~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/MARMOT_Outputs/marmot_cf_23"

# Construct output directory path
output_dir <- path(input_dir, "Graphs")

# Check and create the graph output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

read_input_file <- function(filename, input_dir, delim = "\t") {
  readr::read_delim(fs::path(input_dir, filename), delim = delim, show_col_types = FALSE)
}

# Read all required files
lcaparse_perread <- read_input_file("lcaparse_perread.txt", input_dir)
genome_coverage <- read_input_file("genome_coverage_all.txt", input_dir)
read_numbers     <- read_input_file("read_numbers.tsv", input_dir)

# Created with ~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/scripts/generate_risk_table.py
# Same risk table for every analysis which lives in Pathogen_Database
risk_table <- read.csv("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/risk_table.csv")

# Set plotting theme
custom_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.line.y = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Process data ------------------------------------------------------------

# Risk data ------
# Filter out viruses as their species names aren't correct
risk_table <- risk_table %>%
  filter(Type_of_pest != ("Virus or Viroid"))

# Look at the risk data - see how many species are present more than once
risk_counts <- risk_table %>%
  count(Species, sort = TRUE) %>% 
  filter(n > 1)

# Candidatus Phytoplasma is in there 25 times, then 11 species present 2 - 6 times

# Collapsing to have one species per row, taking the max value of each
# If there are 5 strains of a Bacteria we are assuming worst case scenario
collapsed_risk_table <- risk_table %>%
  group_by(Species) %>%
  summarise(
    Regulation = paste(unique(na.omit(Regulation)), collapse = ";"),
    
    Likelihood_unmitigated = if (all(is.na(Likelihood_unmitigated))) NA_real_ else max(Likelihood_unmitigated, na.rm = TRUE),
    Impact_unmitigated = if (all(is.na(Impact_unmitigated))) NA_real_ else max(Impact_unmitigated, na.rm = TRUE),
    Risk_Rating_unmitigated = if (all(is.na(Risk_Rating_unmitigated))) NA_real_ else max(Risk_Rating_unmitigated, na.rm = TRUE),
    
    Likelihood_mitigated = if (all(is.na(Likelihood_mitigated))) NA_real_ else max(Likelihood_mitigated, na.rm = TRUE),
    Impact_mitigated = if (all(is.na(Impact_mitigated))) NA_real_ else max(Impact_mitigated, na.rm = TRUE),
    Risk_Rating_mitigated = if (all(is.na(Risk_Rating_mitigated))) NA_real_ else max(Risk_Rating_mitigated, na.rm = TRUE),
    
    Regulated = if_else(any(Regulated == "Yes", na.rm = TRUE), "Yes", "No"),
    Natural_Spread = if_else(any(Natural_Spread == "Yes", na.rm = TRUE), "Yes", "No"),
    
    UK = if (length(unique(na.omit(UK))) == 1) unique(na.omit(UK)) else "N/A",
    
    .groups = "drop"
  )


#Lcaparse data ----

# Make anything below species level species instead e.g. subspecies to species 
# This does not handle species group / complex
lcaparse_perread <- lcaparse_perread %>%
  mutate(
    Species = case_when(
      Taxon_Rank == "subspecies" ~ str_extract(Taxon_Name, "^[^ ]+\\s[^ ]+(-[^ ]+)?"),
      TRUE ~ Taxon_Name
    )
  )

# Summarise and filter the file to have one row per species per barcode
# Loose the taxonID in this step as it is complicated by the combined subspecies 
lcaparse <- lcaparse_perread %>%
  group_by(Barcode, Species) %>%
  summarise(
    Taxon_Name = first(Taxon_Name),
    Taxon_Rank = first(Taxon_Rank),
    Read_Count = n_distinct(Read_ID),
    Taxon_ID = first(Taxon_ID),
    Avg_Mean_Identity = mean(Mean_Identity, na.rm = TRUE),
    .groups = "drop"
  )

# merge read numbers & lcaparse on Barcode
read_numbers <- read_numbers %>% 
  rename(TotalReadCount = Read_Count)

lcaparse <- lcaparse %>% 
  left_join(read_numbers, by = "Barcode")  

#Create normalised read count columns 
lcaparse <- lcaparse %>% 
  mutate(
    HP100k = (Read_Count / TotalReadCount) * 100000, 
    Filtered_HP100k = (Read_Count / FilterReadCount) * 100000)


# Combining Read & DEFRA Risk Data    -----
# This will loose the unassigned reads so if I want them in the graph would have to add back in later
lca_risk <- lcaparse %>%
  left_join(collapsed_risk_table, by = "Species") %>% 
  filter(Read_Count >10) # Filter low level species that could be artefacts, may need to play with filter 

# Add a column that groups by risk factor - groupings based on DEFRA documentation
lca_risk  <- lca_risk  %>%
  mutate(
    Risk_Category = case_when(
      Risk_Rating_mitigated <= 14                     ~ "Blue",
      Risk_Rating_mitigated >= 15 & Risk_Rating_mitigated <= 29  ~ "Green",
      Risk_Rating_mitigated >= 30 & Risk_Rating_mitigated <= 44  ~ "Yellow",
      Risk_Rating_mitigated >= 45 & Risk_Rating_mitigated <= 59  ~ "Orange",
      Risk_Rating_mitigated >= 60                              ~ "Red",
      TRUE ~ "Unclassified"
    )
  )

# Improve UK column - unknown to be n/a
lca_risk$UK[is.na(lca_risk$UK) | lca_risk$UK == "Unknown"] <- "N/A"


lca_risk$Risk_Category <- factor(
  lca_risk$Risk_Category,
  levels = c("Red", "Orange", "Yellow", "Green", "Blue", "Unclassified")
)

lca_risk$Barcode <- 
  factor(lca_risk$Barcode, levels = sort(unique(lca_risk$Barcode)))

#Summarise the data to only include species with a DEFRA defined Risk category
risk_only <- lca_risk %>% 
  filter(Risk_Category != ("Unclassified"))

# Defra Risk plots ------
colours <- c(
  Blue = "#66c2a5",
  Green = "#a6d854",
  Yellow = "#ffd92f",
  Orange = "#fc8d62",
  Red = "#e31a1c",
  Unclassified = "grey"
)

# Plot stacked bar chart of defra risk categories --
message("Plotting risk plots...")

# Create new output subdirectory 
risk_dir <- fs::path(output_dir, "defra_risk")
if (!dir.exists(risk_dir)) {
  dir.create(risk_dir, recursive = TRUE)
}

plot_risk_stacked <- function(data, y_col, save_path, colours) {
  message("Plotting graph: ", save_path)
  p <- ggplot(data, aes(x = Barcode, y = .data[[y_col]], fill = Risk_Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours) +
    labs(
      title = "Read Distribution by Defra Risk Category",
      x = "Barcode",
      y = "Reads per 100,000",
      fill = "Risk Category"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    custom_theme
  
  # print(p)
  
  ggsave(save_path, p, width = 12, height = 6)
}

# Plot with the unclassified and without, also with different read normalisation 
plot_risk_stacked(risk_only, y="HP100k",
                  save_path = fs::path(risk_dir, "HP100k_risk.svg"),
                  colours = colours)
plot_risk_stacked(lca_risk, y="HP100k",
                  save_path = fs::path(risk_dir, "HP100k_all.svg"),
                  colours = colours)
plot_risk_stacked(risk_only, y="Filtered_HP100k",
                  save_path = fs::path(risk_dir, "Filtered_HP100k_risk.svg"),
                  colours = colours)
plot_risk_stacked(lca_risk, y="Filtered_HP100k", 
                  save_path = fs::path(risk_dir, "Filtered_HP100k_all.svg"),
                  colours = colours)

# Faceted by Presence
plot_risk_facet <- function(data, save_path, y_col, colours) {
  # Ensure factor levels are consistent
  data$UK <- factor(
    data$UK,
    levels = c("Absent", "Present (Limited)", "Present (Unknown Distribution)",
               "Present (Widespread)", "N/A")
  )
  message("Plotting graph: ", save_path)
  p <- ggplot(data, aes(x = Barcode, y = .data[[y_col]], fill = Risk_Category)) +
    facet_wrap(~ UK, scales = "fixed", ncol =2) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours) +
    labs(
      title = "Read Distribution by Defra Risk Category",
      x = "Barcode",
      y = "Reads per 100,000",
      fill = "Risk Category"
    )  +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    custom_theme
  
  # print(p)
  
  ggsave(save_path, p, width = 14, height = 6)
}

# Plot with the unclassified and without, also with different read normalisation 
plot_risk_facet(risk_only, y="HP100k",
                  save_path = fs::path(risk_dir, "HP100k_risk_facet.svg"),
                  colours = colours)
plot_risk_facet(lca_risk, y="HP100k",
                  save_path = fs::path(risk_dir, "HP100k_all_facet.svg"),
                  colours = colours)
plot_risk_facet(risk_only, y="Filtered_HP100k",
                  save_path = fs::path(risk_dir, "Filtered_HP100k_risk_facet.svg"),
                  colours = colours)
plot_risk_facet(lca_risk, y="Filtered_HP100k", 
                  save_path = fs::path(risk_dir, "Filtered_HP100k_all_facet.svg"),
                  colours = colours)


# 9 Pathogens --------

# Pathogens of interest
target_genera <- c(
  "Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
  "Pyrenophora", "Claviceps", "Parastagonospora")


# Add genus column
lcaparse_genus <- lcaparse %>%
  mutate(
    Genus = if_else(Taxon_Rank == "genus", Taxon_Name, word(Taxon_Name, 1))
  )

# Summarise read count per Barcode × Genus (for ALL genera)
genus_summary_all <- lcaparse_genus %>%
  group_by(Barcode, Genus) %>%
  summarise(
    Read_Count = sum(Read_Count, na.rm = TRUE),
    .groups = "drop"
  )

# Get all the unique barcodes in the df
all_barcodes <- unique(lcaparse$Barcode)

# Create every possible combination of each Barcode with each genus in target_genera
genus_summary <- expand_grid(
  Barcode = all_barcodes,
  Genus = target_genera
) %>%
  # Then join to real data
  left_join(genus_summary_all, by = c("Barcode", "Genus")) %>%
  mutate(Read_Count = replace_na(Read_Count, 0)) # replacing NA with 0

# Add TotalReadCount and FilterReadCount to calculate normalised values
genus_summary <- genus_summary %>%
  left_join(
    lcaparse %>% select(Barcode, TotalReadCount, FilterReadCount) %>% distinct(),
    by = "Barcode"
  ) %>%
  mutate(
    HP100k = (Read_Count / TotalReadCount) * 100000,
    Filtered_HP100k = (Read_Count / FilterReadCount) * 100000
  )

# Export pathogen tsv  ----
# Clean up the output
genus_export <- genus_summary %>%
  select(-TotalReadCount, -FilterReadCount)

# Define save path - input directory so it will live with the other summary 
save_path <- fs::path(input_dir, "genus_summary.tsv")

# Write to TSV
write_tsv(genus_export, save_path)


message("Genus summary saved to: ", save_path)




# Pathogen Plots ----
pathogen_colours <- setNames(brewer.pal(9, "Set3"), target_genera)

# Create subdirectory for these to be saved into
pathogen_dir <- fs::path(output_dir, "pathogen_graphs")
if (!dir.exists(pathogen_dir)) {
  dir.create(pathogen_dir, recursive = TRUE)
}

# Individual plot per pathogen 
message("Plotting pathogen plots...")

plot_pathogen_bar <- function(data,
                              pathogen, 
                              y_col, 
                              save_path, 
                              pathogen_colours,
                              y_lab = "Reads per 100,000",
                              error_col = NULL) {
  # Filter for the specified pathogen
  pathogen_data <- data %>%
    filter(Genus == pathogen)
  
  # Extract fill colour for this pathogen
  fill_colour <- pathogen_colours[[pathogen]]
  
  # Determine y_axis upper limit
  y_max <- max(pathogen_data[[y_col]], na.rm = TRUE)
  y_upper <- if (y_max == 0) 1 else NA
  
  # Plot
  p <- ggplot(pathogen_data, aes(x = Barcode, y = .data[[y_col]])) +
    geom_bar(stat = "identity", fill = fill_colour) 
  
  # Add error bars if specified
  if (!is.null(error_col) && error_col %in% colnames(pathogen_data)) {
    p <- p + geom_errorbar(
      aes(
        ymin = .data[[y_col]] - .data[[error_col]],
        ymax = .data[[y_col]] + .data[[error_col]]
      ),
      width = 0.3,
      colour = "black"
    )
  }
  
  p <- p +  labs(
      title = pathogen,
      x = "Barcode",
      y = y_lab
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_upper)) +
    custom_theme
  
  ggsave(save_path, p, width = 10, height = 6)
  message("Plot saved: ", save_path)
}

#Run for the two y cols and all genera
y_cols <- c("HP100k", "Filtered_HP100k")

for (pathogen in target_genera) {
  for (y in y_cols) {
    filename <- paste0(pathogen, "_", y, ".svg")
    save_path <- fs::path(pathogen_dir, filename)
    
    plot_pathogen_bar(
      data = genus_summary,
      pathogen = pathogen,
      y_col = y,
      save_path = save_path,
      pathogen_colours = pathogen_colours
    )
  }
}

# Facet of all those plots for both y cols 
scales_options <- c("fixed", "free_y")

for (y in y_cols) {
  for (sc in scales_options) {
    
    # Create the plot
    p <- ggplot(genus_summary, aes(x = Barcode, y = .data[[y]], fill = Genus)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ Genus, scales = sc, ncol = 3) +
      scale_fill_manual(values = pathogen_colours) +
      labs(
        title = paste("Pathogen Abundance by Genus (", y, ", ", sc, ")", sep = ""),
        x = "Barcode",
        y = "Reads per 100,000",
        fill = "Genus"
      ) +
      scale_y_continuous(expand = c(0, 0), limits = if (sc == "fixed") c(0, NA) else NULL) +
      custom_theme +  
      theme(legend.position = "none")
    
    # Construct filename and path
    filename <- paste0("facet_pathogen_", y, "_", sc, ".svg")
    save_path <- fs::path(pathogen_dir, filename)
    
    # Save
    ggsave(save_path, p, width = 16, height = 8)
    message("Facet plot saved: ", save_path)
  }
}

# Genome coverage of pathgoens ----

# Create a single regex pattern that matches any of the genera
genera_pattern <- str_c(target_genera, collapse = "|")

# Filter the lcaparse_perread to only include the target genera
perread_target <- lcaparse_perread %>%
  filter(str_detect(Taxon_Name, genera_pattern)) %>%
  filter(!is.na(Taxon_ID), Taxon_ID != 0)


# Prepare unique Barcode + Taxon_ID rows with taxon metadata
# one row per barcode & taxaID
perread_taxa_unique <- perread_target %>%
  group_by(Barcode, Taxon_ID, Taxon_Name, Taxon_Rank) %>%
  summarise(
    Mean_Identity = mean(Mean_Identity, na.rm = TRUE),
    .groups = "drop"
  )

# Join with genome_coverage using Barcode + taxaID
target_coverage <- genome_coverage %>%
  inner_join(
    perread_taxa_unique,
    by = c("Barcode" = "Barcode", "taxaID" = "Taxon_ID")
  )

# Summarise to get mean + standard error of coverage per Barcode × Genus
# For cases when there are multiple species within the genera
coverage_summary <- target_coverage %>%
  mutate(Genus = word(Taxon_Name, 1)) %>%
  group_by(Barcode, Genus) %>%
  summarise(
    Mean_Coverage = mean(coverage_percentage, na.rm = TRUE),
    SE_Coverage = sd(coverage_percentage, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Expand to include missing Barcode × Genus pairs with 0 coverage
barcode_genus_grid <- expand_grid(
  Barcode = unique(coverage_summary$Barcode),
  Genus = target_genera
)

# Full dataset 
coverage_plot_data <- barcode_genus_grid %>%
  left_join(coverage_summary, by = c("Barcode", "Genus")) %>%
  mutate(
    Mean_Coverage = replace_na(Mean_Coverage, 0),
    SE_Coverage = replace_na(SE_Coverage, 0)
  )

# Plot coverage of pathogens ----
message("Plotting pathogen genome coverage plots...")
# Calling the earlier function
for (genus in target_genera) {
  save_path <- fs::path(pathogen_dir, paste0(genus, "_genome_coverage.svg"))
  
  plot_pathogen_bar(
    data = coverage_plot_data,
    pathogen = genus,
    y_col = "Mean_Coverage",
    y_lab = "Average Genome Coverage (%)",
    error_col = "SE_Coverage",
    save_path = save_path,
    pathogen_colours = pathogen_colours
  )
}

# Facetted

for (sc in scales_options) {
  
  # Fallback for zero-only data
  y_max <- max(coverage_plot_data$Mean_Coverage, na.rm = TRUE)
  y_upper <- if (y_max == 0) 1 else NA
  
  # Create the plot
  p <- ggplot(coverage_plot_data, aes(x = Barcode, y = Mean_Coverage, fill = Genus)) +
    geom_bar(stat = "identity") +
    geom_errorbar(
      aes(ymin = Mean_Coverage - SE_Coverage,
          ymax = Mean_Coverage + SE_Coverage),
      width = 0.3,
      colour = "black"
    ) +
    facet_wrap(~ Genus, scales = sc, ncol = 3) +
    scale_fill_manual(values = pathogen_colours) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = if (sc == "fixed") c(0, y_upper) else NULL
    ) +
    labs(
      title = paste("Genome Coverage by Genus (", sc, " scale)", sep = ""),
      x = "Barcode",
      y = "Average Genome Coverage (%)",
      fill = "Genus"
    ) +
    custom_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
  
  # Construct filename and save
  filename <- paste0("facet_coverage_", sc, ".svg")
  save_path <- fs::path(pathogen_dir, filename)
  
  ggsave(save_path, p, width = 16, height = 8)
  message("Facet plot saved: ", save_path)
}



# Genus level heatmap --------
message("Plotting heatmaps...")
heat_dir <- fs::path(output_dir, "heatmap")
if (!dir.exists(heat_dir)) {
  dir.create(heat_dir, recursive = TRUE)
}

plot_genus_heatmap <- function(data,
                               read_numbers,
                               y_col = "HP100k",
                               min_reads = 50,
                               output_dir = heat_dir,
                               palette_option = "D") {
  
  # Filter and join total reads
  filtered_data <- data %>%
    filter(!is.na(Genus), Genus != "Unassigned") %>%
    group_by(Genus) %>%
    filter(sum(Read_Count, na.rm = TRUE) >= min_reads) %>%
    ungroup() %>%
    left_join(read_numbers, by = "Barcode") %>%
    mutate(
      HP100k = (Read_Count / TotalReadCount) * 100000,
      Filtered_HP100k = (Read_Count / FilterReadCount) * 100000
    )
  
  # Reorder genera by abundance
  genus_order <- filtered_data %>%
    group_by(Genus) %>%
    summarise(total = sum(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Genus)
  
  plot_data <- filtered_data %>%
    mutate(Genus = factor(Genus, levels = genus_order))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Barcode, y = Genus, fill = .data[[y_col]])) +
    geom_tile() +
    scale_fill_viridis_c(
      trans = "log1p",
      name = "Reads per 100k",
      breaks = c(1, 10, 100, 1000, 10000),
      labels = scales::comma_format()(c(1, 10, 100, 1000, 10000)),
      option = palette_option
    ) +
    labs(
      title = paste("Genus Heatmap (", y_col, ", min reads: ", min_reads, ")", sep = ""),
      x = "Barcode",
      y = "Genus"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      legend.key.size = unit(1.2, "cm")
    )
  
  # Save plot
  filename <- paste0("genus_heatmap_", y_col, "_min", min_reads, ".svg")
  save_path <- fs::path(heat_dir, filename)
  ggsave(save_path, p, width = 12, height = 10)
  message("Saved heatmap: ", save_path)
  
  return(p)
}

# Plot both y_cols with both thresholds
for (y in c("HP100k", "Filtered_HP100k")) {
  for (min_reads in c(50, 100)) {
    plot_genus_heatmap(
      data = genus_summary_all,
      read_numbers = read_numbers,
      y_col = y,
      min_reads = min_reads,
      output_dir = heat_dir
    )
  }
}


# Pull out the species and maybe ReadIDs to consider in more depth -----
extract_red_risk_reads <- function(data_risk,
                                   data_perread,
                                   output_dir,
                                   exclude_widespread = FALSE,
                                   filename = "RedRisk_ReadIDs.tsv") {
  # Base filtering
  filtered <- data_risk %>%
    filter(
      Read_Count > 10,
      Risk_Category == "Red",
      Avg_Mean_Identity > 90
    )
  
  # Optional additional filter
  if (exclude_widespread) {
    filtered <- filtered %>% filter(UK != "Present (Widespread)")
  }
  
  # Get read IDs with metadata
  barcode_taxa <- filtered %>%
    distinct(Barcode, Taxon_Name)
  
  read_ids <- data_perread %>%
    semi_join(barcode_taxa, by = c("Barcode", "Taxon_Name")) %>%
    select(Barcode, Read_ID, Taxon_Name) %>%
    left_join(
      filtered %>% select(Barcode, Taxon_Name, Risk_Category, UK),
      by = c("Barcode", "Taxon_Name")
    )
  
  # Save file
  save_path <- fs::path(output_dir, filename)
  write_tsv(read_ids, save_path)
  
  message("Saved ", nrow(read_ids), " read IDs to: ", save_path)
}

# Red risk not considering presence in the UK
extract_red_risk_reads(
  data_risk = risk_only,
  data_perread = lcaparse_perread,
  output_dir = output_dir,
  exclude_widespread = FALSE,
  filename = "RedRisk_ReadIDs_all.tsv"
)

# Red Risk, filtering out those already known to be widespread present in the UK
extract_red_risk_reads(
  data_risk = risk_only,
  data_perread = lcaparse_perread,
  output_dir = output_dir,
  exclude_widespread = TRUE,
  filename = "RedRisk_ReadIDs_noWidespread.tsv"
)


# Create a summary file of species of interest -----
# Filter risk_only for red-risk species of interest 
red_risk_summary <- risk_only %>%
  filter(
    Risk_Category == "Red",
    Read_Count > 10,
    Avg_Mean_Identity > 90
  ) %>%
  select(
    Barcode, Taxon_Name, Taxon_Rank, Taxon_ID,
    Risk_Category, UK,
    Read_Count, Avg_Mean_Identity
  ) %>%
  distinct() 

#  Join on genome coverage using Barcode + Taxon_ID
red_risk_summary <- red_risk_summary %>%
  left_join(
    genome_coverage %>%
      select(Barcode, taxaID, coverage_percentage, num_reads),
    by = c("Barcode" = "Barcode", "Taxon_ID" = "taxaID")
  ) %>%
  rename(
    Genome_Coverage = coverage_percentage,
    Coverage_Reads = num_reads
  )

# Save
write_tsv(red_risk_summary, fs::path(output_dir, "RedRisk_Species_Summary.tsv"))

message("Saved ", red_risk_summary, " table to: ", output_dir, "/RedRisk_Species_Summary.tsv")

# To finish 
message("R script complete")

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database")

# Packages 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readr)
library(fs)  # for file path safety
library(stringr)

# Set paths and read in the data -----
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Expecting: Rscript create_risk_plots.R <input_dir>
if (length(args) < 1) {
  stop("Usage: Rscript create_risk_plots.R <input_dir>")
}

# input_dir is where the downloaded outputs from the pipeline have been saved
input_dir <- args[1]

# Alternatively provide the path, when not running on command line 
# input_dir <- "~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/MARMOT_Outputs/test"

# Construct output directory path
output_dir <- path(input_dir, "Graphs")

# Check and create the graph output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

read_input_file <- function(filename, input_dir, delim = "\t") {
  readr::read_delim(fs::path(input_dir, filename), delim = delim)
}

# Read all required files
lcaparse_perread <- read_input_file("lcaparse_perread.txt", input_dir)
genome_coverage <- read_input_file("genome_coverage_all.txt", input_dir)
read_numbers     <- read_input_file("read_numbers.tsv", input_dir)

# Created with ~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/scripts/generate_risk_table.py
# Same risk table for every analysis which lives in Pathogen_Database
risk_table <- read.csv("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/risk_table.csv")


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
    Avg_Mean_Identity = mean(Mean_Identity, na.rm = TRUE)
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
  filter(Read_Count >1) # Filter low level species that could be artefacts, may need to play with filter 

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
    theme_minimal()
  
  # print(p)
  
  ggsave(save_path, p, width = 12, height = 6)
}

# Plot with the unclassified and without, also with different read normalisation 
plot_risk_stacked(risk_only, y="HP100k",
                  save_path = fs::path(output_dir, "HP100k_risk.svg"),
                  colours = colours)
plot_risk_stacked(lca_risk, y="HP100k",
                  save_path = fs::path(output_dir, "HP100k_all.svg"),
                  colours = colours)
plot_risk_stacked(risk_only, y="Filtered_HP100k",
                  save_path = fs::path(output_dir, "Filtered_HP100k_risk.svg"),
                  colours = colours)
plot_risk_stacked(lca_risk, y="Filtered_HP100k", 
                  save_path = fs::path(output_dir, "Filtered_HP100k_all.svg"),
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
    facet_wrap(~ UK, scales = "free_y") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours) +
    labs(
      title = "Read Distribution by Defra Risk Category",
      x = "Barcode",
      y = "Reads per 100,000",
      fill = "Risk Category"
    ) +
    theme_minimal()
  
  # print(p)
  
  ggsave(save_path, p, width = 16, height = 6)
}

# Plot with the unclassified and without, also with different read normalisation 
plot_risk_facet(risk_only, y="HP100k",
                  save_path = fs::path(output_dir, "HP100k_risk_facet.svg"),
                  colours = colours)
plot_risk_facet(lca_risk, y="HP100k",
                  save_path = fs::path(output_dir, "HP100k_all_facet.svg"),
                  colours = colours)
plot_risk_facet(risk_only, y="Filtered_HP100k",
                  save_path = fs::path(output_dir, "Filtered_HP100k_risk_facet.svg"),
                  colours = colours)
plot_risk_facet(lca_risk, y="Filtered_HP100k", 
                  save_path = fs::path(output_dir, "Filtered_HP100k_all_facet.svg"),
                  colours = colours)


# 9 Pathogens I am particularly interested in --------

# Pathogens of interest
target_genera <- c(
  "Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
  "Pyrenophora", "Claviceps", "Parastagonospora", "Phaeosphaeria")

# Create a single regex pattern that matches any of the genera
genera_pattern <- str_c(target_genera, collapse = "|")

# Filter lcaparse to keep only rows matching any of the target genera
lcaparse_target <- lcaparse %>%
  filter(str_detect(Taxon_Name, genera_pattern))

# Hard to know how I would want to plot this with such limited data 
# Maybe another bar chart for each pathogen by barcode or a heatmap just of these pathogens 









# # Older code -----
# Using the summary file, no longer using this but the per read
# 
# # Species data 
# #This specific file was copied 
# lca_parse <- read.csv("lcaparse_summary.txt", sep = "\t")
# 
# # Extract species from lca_parse
# species_df <- lca_parse %>%
#   filter(Taxon_Rank == "species") %>%
#   mutate(Species = str_split(Taxon_Path, ",", simplify = FALSE) %>%
#            map_chr(~ tail(.x, 1)))
# 
# # Get all non-species rows (higher taxa)
# higher_taxa_df <- lca_parse %>%
#   filter(Taxon_Rank != "species") %>%
#   group_by(Barcode) %>%
#   summarise(
#     Read_Count = sum(as.numeric(Read_Count), na.rm = TRUE),
#     Percentage_of_Reads = sum(as.numeric(Percentage_of_Reads), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Species = "Higher Taxa",
#     Taxon_Rank = "higher",
#     Taxon_Path = NA  # or "Higher-level taxon"
#   )
# 
# # Combine species & higher taxa dfs
# combined_df <- bind_rows(species_df, higher_taxa_df) 


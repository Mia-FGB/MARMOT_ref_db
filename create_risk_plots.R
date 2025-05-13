setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database")

# Packages 
library(ggplot2)
library(dplyr)
library(tidyverse)


# Risk data ------

# Created with ~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/scripts/generate_risk_table.py
# Should match the pathogen database used for the analysis
risk_table <- read.csv("risk_table.csv")

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
# Noticed an issue with this file that stems back to the pipeline so need to go and fix that!
lcaparse_perread <- 
  read_delim("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/MARMOT_Outputs/cf_23_subsamp/lcaparse_perread.txt",
             delim = "\t")



# Merge the lcaparse & risk data   -------       
merged_df <- left_join(combined_df, collapsed_risk_table, by = "Species") %>% 
  filter(Read_Count >10) %>% # Filter low level species that could be artefacts, may need to play with filter 
  filter(Barcode <=39) # Filter to jsut be barcode 1 - 39 as this is CF 2023 data (won't do this in future)

# Add a column that groups by risk factor - groupings based on DEFRA documentation

merged_df  <- merged_df  %>%
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
merged_df$UK[is.na(merged_df$UK) | merged_df$UK == "Unknown"] <- "N/A"


merged_df$Risk_Category <- factor(
  merged_df$Risk_Category,
  levels = c("Red", "Orange", "Yellow", "Green", "Blue", "Unclassified")
)

merged_df$Barcode <- 
  factor(merged_df$Barcode, levels = sort(unique(merged_df$Barcode)))

#dataframe with only data that have a risk category
risk_only <- merged_df %>% 
  filter(Risk_Category != ("Unclassified"))

# Outputs ----
colours <- c(
  Blue = "#66c2a5",
  Green = "#a6d854",
  Yellow = "#ffd92f",
  Orange = "#fc8d62",
  Red = "#e31a1c",
  Unclassified = "grey"
)

# Plot stacked bar chart
stacked <- ggplot(
  merged_df, # To plot with unclassified 
  # risk_only,   # To plot without unclassified 
       aes(x = Barcode, y = Percentage_of_Reads, fill = Risk_Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  labs(
    title = "Read Distribution by Defra Risk Category",
    x = "Barcode",
    y = "Percentage of Reads",
    fill = "Risk Category"
  ) +
  theme_minimal() 

ggsave("test_plot.svg", stacked, width = 12, height =6 )

# Faceted by Presence 
risk_only$UK<- factor(
  risk_only$UK,
  levels = c("Absent", "Present (Limited)", "Present (Unknown Distribution)",
             "Present (Widespread)", "N/A")
)

presence <- ggplot(risk_only,
  aes(x = Barcode, y = Percentage_of_Reads, fill = Risk_Category)) +
  facet_wrap(~ UK, scales = "free_y") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  labs(
    title = "Read Distribution by Defra Risk Category",
    x = "Barcode",
    y = "Percentage of Reads",
    fill = "Risk Category"
  ) +
  theme_minimal()

ggsave("test_facet_plot.svg", presence, width = 16, height =6 )


# Maybe a file with some general stats about the data e.g. which had the highest risk rating 







# # Older code - to be deleted once all working
# # Was using the summary file instead want to be per read
# 
# # Species data -------
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


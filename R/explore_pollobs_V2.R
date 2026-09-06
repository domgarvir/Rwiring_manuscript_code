library(dplyr)
library(tidyr)
library(iNEXT)

# -------------------------------------------------------------------------
# 1. Load and inspect data
# -------------------------------------------------------------------------

weekly_data <- readRDS(
  "DATA/PollObs/weekly_data_for_rewiring.RDS"
)

colnames(weekly_data)

weekly_data %>%
  count(Sampling_type, Botanical_garden)

# -------------------------------------------------------------------------
# 2. Pool focal and complementary observations
# -------------------------------------------------------------------------

weekly_interactions <- weekly_data %>%
  filter(
    !is.na(Botanical_garden),
    !is.na(Week),
    !is.na(Plant),
    !is.na(Pollinator),
    !is.na(Total_pair_interactions),
    Total_pair_interactions > 0
  ) %>%
  group_by(
    Botanical_garden,
    Week,
    Plant,
    Pollinator
  ) %>%
  summarise(
    Total_pair_interactions = sum(
      Total_pair_interactions,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Pair = paste(Plant, Pollinator, sep = "___")
  ) %>%
  arrange(Botanical_garden, Week, Plant, Pollinator)

weekly_interactions

# -------------------------------------------------------------------------
# 3. Overall dataset summary
# -------------------------------------------------------------------------

overall_summary <- weekly_interactions %>%
  summarise(
    Plant_species = n_distinct(Plant),
    Floral_visitor_species = n_distinct(Pollinator),
    Interaction_events = sum(Total_pair_interactions),
    Unique_interactions = n_distinct(Plant, Pollinator),
    Weekly_networks = n_distinct(Botanical_garden, Week)
  )

overall_summary

# -------------------------------------------------------------------------
# 4. Summary by botanical garden
# -------------------------------------------------------------------------

pollobs_summary <- weekly_interactions %>%
  group_by(Botanical_garden) %>%
  summarise(
    Plant_species = n_distinct(Plant),
    Floral_visitor_species = n_distinct(Pollinator),
    Unique_interactions = n_distinct(Plant, Pollinator),
    Interaction_events = sum(Total_pair_interactions),
    Weekly_networks = n_distinct(Week),
    .groups = "drop"
  ) %>%
  arrange(Botanical_garden)

pollobs_summary

# -------------------------------------------------------------------------
# 5. Sample coverage for each weekly network
# -------------------------------------------------------------------------

get_sample_coverage <- function(pair_abundances) {
  
  pair_abundances <- pair_abundances[
    is.finite(pair_abundances) &
      pair_abundances > 0
  ]
  
  if (length(pair_abundances) == 0) {
    return(NA_real_)
  }
  
  coverage_info <- iNEXT::DataInfo(
    pair_abundances,
    datatype = "abundance"
  )
  
  as.numeric(coverage_info$SC[[1]])
}

weekly_coverage <- weekly_interactions %>%
  group_by(Botanical_garden, Week) %>%
  summarise(
    Observed_interactions = n(),
    Interaction_events = sum(Total_pair_interactions),
    Sample_coverage = get_sample_coverage(
      Total_pair_interactions
    ),
    .groups = "drop"
  ) %>%
  arrange(Botanical_garden, Week)

weekly_coverage

# -------------------------------------------------------------------------
# 6. Overall coverage summary
# -------------------------------------------------------------------------

coverage_summary <- weekly_coverage %>%
  summarise(
    Networks = n(),
    Mean_coverage = mean(Sample_coverage, na.rm = TRUE),
    SD_coverage = sd(Sample_coverage, na.rm = TRUE),
    Median_coverage = median(Sample_coverage, na.rm = TRUE),
    Minimum_coverage = min(Sample_coverage, na.rm = TRUE),
    Maximum_coverage = max(Sample_coverage, na.rm = TRUE)
  )

coverage_summary

# -------------------------------------------------------------------------
# 7. Coverage summary by botanical garden
# -------------------------------------------------------------------------

coverage_by_garden <- weekly_coverage %>%
  group_by(Botanical_garden) %>%
  summarise(
    Networks = n(),
    Mean_coverage = mean(Sample_coverage, na.rm = TRUE),
    SD_coverage = sd(Sample_coverage, na.rm = TRUE),
    Median_coverage = median(Sample_coverage, na.rm = TRUE),
    Minimum_coverage = min(Sample_coverage, na.rm = TRUE),
    Maximum_coverage = max(Sample_coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Botanical_garden)

coverage_by_garden

# -------------------------------------------------------------------------
# 8. Coverage table for the manuscript
# -------------------------------------------------------------------------

coverage_table <- weekly_coverage %>%
  select(
    Botanical_garden,
    Week,
    Sample_coverage
  ) %>%
  mutate(
    Sample_coverage = round(Sample_coverage, 2)
  ) %>%
  pivot_wider(
    names_from = Week,
    values_from = Sample_coverage,
    names_sort = TRUE
  ) %>%
  left_join(
    coverage_by_garden %>%
      transmute(
        Botanical_garden,
        mean = round(Mean_coverage, 2),
        std = round(SD_coverage, 2)
      ),
    by = "Botanical_garden"
  ) %>%
  arrange(Botanical_garden)

coverage_table

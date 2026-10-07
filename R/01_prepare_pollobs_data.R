library(dplyr)
library(lubridate)
library(iNEXT)
library(tidyr)

# -------------------------------------------------------------------------
# 1. Load data and define temporal variables
# -------------------------------------------------------------------------

interaction_data <- readRDS(
  "DATA/PollObs/raw/interaction_data.RDS"
) %>%
  mutate(
    ISO_year = isoyear(Date),
    Week = isoweek(Date),
    Y_day = yday(Date)
  )

# Inspect values before filtering.
interaction_data %>%
  count(Sampling)

interaction_data %>%
  count(
    Plant_rank,
    Pollinator_rank
  )


# -------------------------------------------------------------------------
# 2. Define the sampling schedule
# -------------------------------------------------------------------------

# Include all sampling occasions, even those with no recorded interactions.
sampling_schedule <- interaction_data %>%
  filter(
    Sampling %in% c(
      "Focal",
      "Random_census"
    ),
    !is.na(Botanical_garden),
    !is.na(Date)
  ) %>%
  distinct(
    Botanical_garden,
    ISO_year,
    Week,
    Date,
    Y_day
  )

# Check whether a garden was sampled on more than one date in the same week.
sampling_dates_check <- sampling_schedule %>%
  count(
    Botanical_garden,
    ISO_year,
    Week,
    name = "n_dates"
  ) %>%
  filter(n_dates > 1)

sampling_dates_check

stopifnot(nrow(sampling_dates_check) == 0)

# -------------------------------------------------------------------------
# 3. Retain valid plant--floral visitor interactions
# -------------------------------------------------------------------------

pollobs_interactions <- interaction_data %>%
  filter(
    Sampling %in% c(
      "Focal",
      "Random_census"
    ),
    Plant_rank == "SPECIES",
    Pollinator_rank == "SPECIES",
    !is.na(Plant_accepted_name),
    !is.na(Pollinator_accepted_name),
    Pollinator != "None",
    !is.na(Interactions),
    Interactions > 0
  )

# -------------------------------------------------------------------------
# 4. Pool focal and complementary observations within weekly networks
# -------------------------------------------------------------------------

weekly_interactions <- pollobs_interactions %>%
  group_by(
    Botanical_garden,
    ISO_year,
    Week,
    Date,
    Y_day,
    Plant = Plant_accepted_name,
    Floral_visitor = Pollinator_accepted_name
  ) %>%
  summarise(
    n_visits = sum(
      Interactions,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(n_visits > 0)

# -------------------------------------------------------------------------
# 5. Dataset summary
# -------------------------------------------------------------------------

network_summary <- weekly_interactions %>%
  group_by(Botanical_garden) %>%
  summarise(
    Plant_species = n_distinct(Plant),
    Floral_visitor_species = n_distinct(
      Floral_visitor
    ),
    Unique_interactions = n_distinct(
      Plant,
      Floral_visitor
    ),
    Interaction_events = sum(n_visits),
    Weekly_networks = n_distinct(
      ISO_year,
      Week
    ),
    .groups = "drop"
  )

network_summary

overall_summary <- weekly_interactions %>%
  summarise(
    Plant_species = n_distinct(Plant),
    Floral_visitor_species = n_distinct(
      Floral_visitor
    ),
    Unique_interactions = n_distinct(
      Plant,
      Floral_visitor
    ),
    Interaction_events = sum(n_visits),
    Weekly_networks = n_distinct(
      Botanical_garden,
      ISO_year,
      Week
    )
  )

overall_summary

# -------------------------------------------------------------------------
# 6. Prepare the binary edge list for Rwiring
# -------------------------------------------------------------------------

interactions_rwiring <- weekly_interactions %>%
  transmute(
    Site = Botanical_garden,
    Time = Y_day,
    Origin_node = Plant,
    Destiny_node = Floral_visitor,
    Origin_guild = "plant",
    Destiny_guild = "floral_visitor"
  ) %>%
  distinct()

prepared <- Rwiring::prepare_interactions(
  interactions_rwiring,
  directed = FALSE,
  multipartite = TRUE
)

summary(prepared)

network_index <- weekly_interactions %>%
  distinct(
    Site = Botanical_garden,
    ISO_year,
    Week,
    Date,
    Time = Y_day
  ) %>%
  arrange(
    Site,
    Date
  )

# -------------------------------------------------------------------------
# 7. Define comparisons between successive non-empty networks
# -------------------------------------------------------------------------

comparisons <- network_index %>%
  group_by(Site) %>%
  transmute(
    ISO_year,
    Week_1 = Week,
    Week_2 = lead(Week),
    Date_1 = Date,
    Date_2 = lead(Date),
    Time_1 = Time,
    Time_2 = lead(Time)
  ) %>%
  ungroup() %>%
  filter(!is.na(Time_2)) %>%
  mutate(
    Elapsed_days = Time_2 - Time_1,
    Mid_day = (Time_1 + Time_2) / 2,
    Transition = paste0(
      Week_1,
      "\u2192",
      Week_2
    )
  )

comparisons

# -------------------------------------------------------------------------
# 8. Sanity checks
# -------------------------------------------------------------------------

# Number of non-empty weekly networks and comparisons.
comparison_summary <- comparisons %>%
  group_by(Site) %>%
  summarise(
    Networks = n() + 1L,
    Comparisons = n(),
    Minimum_interval = min(Elapsed_days),
    Median_interval = median(Elapsed_days),
    Maximum_interval = max(Elapsed_days),
    .groups = "drop"
  )

comparison_summary

# Unusually short or long intervals.
comparisons %>%
  filter(
    Elapsed_days < 5 |
      Elapsed_days > 10
  )

# Each edge must occur only once per network.
duplicated_edges <- interactions_rwiring %>%
  count(
    Site,
    Time,
    Origin_node,
    Destiny_node
  ) %>%
  filter(n > 1)

duplicated_edges

stopifnot(nrow(duplicated_edges) == 0)

# Sampling occasions that did not produce a non-empty network.
empty_sampling_periods <- sampling_schedule %>%
  anti_join(
    network_index,
    by = c(
      "Botanical_garden" = "Site",
      "ISO_year",
      "Week",
      "Date",
      "Y_day" = "Time"
    )
  )

empty_sampling_periods




# -------------------------------------------------------------------------
# Sample coverage of each non-empty network
# -------------------------------------------------------------------------

get_sample_coverage <- function(abundances) {
  
  abundances <- abundances[
    is.finite(abundances) &
      abundances > 0
  ]
  
  if (length(abundances) == 0) {
    return(NA_real_)
  }
  
  coverage_information <- iNEXT::DataInfo(
    abundances,
    datatype = "abundance"
  )
  
  coverage_information$SC[[1]]
}

weekly_coverage <- weekly_interactions %>%
  group_by(
    Botanical_garden,
    ISO_year,
    Week,
    Date
  ) %>%
  summarise(
    Observed_interactions = n(),
    Interaction_events = sum(n_visits),
    Sample_coverage = get_sample_coverage(
      n_visits
    ),
    .groups = "drop"
  )

weekly_coverage

coverage_by_garden <- weekly_coverage %>%
  group_by(Botanical_garden) %>%
  summarise(
    Networks = n(),
    Mean = mean(
      Sample_coverage,
      na.rm = TRUE
    ),
    SD = sd(
      Sample_coverage,
      na.rm = TRUE
    ),
    Median = median(
      Sample_coverage,
      na.rm = TRUE
    ),
    Minimum = min(
      Sample_coverage,
      na.rm = TRUE
    ),
    Maximum = max(
      Sample_coverage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

coverage_by_garden

coverage_overall <- weekly_coverage %>%
  summarise(
    Networks = n(),
    Mean = mean(
      Sample_coverage,
      na.rm = TRUE
    ),
    SD = sd(
      Sample_coverage,
      na.rm = TRUE
    ),
    Median = median(
      Sample_coverage,
      na.rm = TRUE
    ),
    Minimum = min(
      Sample_coverage,
      na.rm = TRUE
    ),
    Maximum = max(
      Sample_coverage,
      na.rm = TRUE
    )
  )

coverage_overall

coverage_table <- weekly_coverage %>%
  select(
    Botanical_garden,
    Week,
    Sample_coverage
  ) %>%
  mutate(
    Sample_coverage = sprintf(
      "%.2f",
      Sample_coverage
    )
  ) %>%
  pivot_wider(
    names_from = Week,
    values_from = Sample_coverage
  ) %>%
  left_join(
    coverage_by_garden %>%
      transmute(
        Botanical_garden,
        mean = sprintf("%.2f", Mean),
        SD = sprintf("%.2f", SD)
      ),
    by = "Botanical_garden"
  )

coverage_table


# -------------------------------------------------------------------------
# 9. Save processed data for descriptive analyses
# -------------------------------------------------------------------------

output_dir <- file.path(
  "DATA",
  "PollObs",
  "processed"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pollobs_descriptive_data <- list(
  
  # Sampling occasions, including the two occasions without interactions.
  sampling_schedule = sampling_schedule,
  
  # Valid species-level records before weekly aggregation.
  pollobs_interactions = pollobs_interactions,
  
  # Canonical weekly interaction table, including frequencies.
  weekly_interactions = weekly_interactions,
  
  # Binary edge list supplied to Rwiring.
  interactions_rwiring = interactions_rwiring,
  
  # Correspondence between network identifiers, dates and calendar weeks.
  network_index = network_index,
  
  # Comparisons between successive non-empty networks.
  comparisons = comparisons,
  
  # Sampling occasions without floral visits.
  empty_sampling_periods = empty_sampling_periods,
  
  # Network-level sample coverage in numeric format.
  weekly_coverage = weekly_coverage,
  
  # Dataset summaries.
  network_summary = network_summary,
  overall_summary = overall_summary,
  comparison_summary = comparison_summary,
  coverage_by_garden = coverage_by_garden,
  coverage_overall = coverage_overall
)

saveRDS(
  pollobs_descriptive_data,
  file.path(
    output_dir,
    "pollobs_descriptive_data.RDS"
  )
)

saveRDS(
  weekly_interactions,
  file.path(
    output_dir,
    "pollobs_weekly_interactions.RDS"
  )
)

saveRDS(
  sampling_schedule,
  "DATA/PollObs/processed/pollobs_sampling_schedule.RDS"
)

write.csv(
  network_summary,
  file.path(
    output_dir,
    "pollobs_network_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparison_summary,
  file.path(
    output_dir,
    "pollobs_comparison_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  weekly_coverage,
  file.path(
    output_dir,
    "pollobs_weekly_coverage.csv"
  ),
  row.names = FALSE
)

write.csv(
  coverage_by_garden,
  file.path(
    output_dir,
    "pollobs_coverage_by_garden.csv"
  ),
  row.names = FALSE
)

write.csv(
  empty_sampling_periods,
  file.path(
    output_dir,
    "pollobs_empty_sampling_periods.csv"
  ),
  row.names = FALSE
)

# Prepare the Clough plant--pollinator networks
#
# This script creates the canonical site-level weighted networks used in the
# spatial case study. Observations from all sampling rounds and both years are
# pooled within each of the 18 sites. Link strength is the total number of
# recorded visits for a plant-pollinator pair.

library(dplyr)
library(tidyr)
library(lubridate)
library(iNEXT)
library(geosphere)

# =============================================================================
# 1. Files and analysis settings
# =============================================================================

input_file <- file.path(
  "DATA",
  "Clough",
  "clough_data.rds"
)

# Optional. If the EuPPollNet flower-count table is available, it provides a
# more complete sampling schedule because sampling occasions with no recorded
# interactions are absent from the interaction table.
flower_counts_file <- file.path(
  "DATA",
  "Clough",
  "Flower_counts.csv"
)

output_dir <- file.path(
  "DATA",
  "Clough",
  "processed"
)

# "all_accepted" retains every accepted taxonomic identification and matches
# the current Clough processing. "species_only" retains records for which both
# partners were identified to species level.
taxonomic_scope <- "all_accepted"

stopifnot(
  taxonomic_scope %in% c(
    "all_accepted",
    "species_only"
  )
)

# =============================================================================
# 2. Load data and check the expected structure
# =============================================================================

clough_raw <- readRDS(input_file)

required_columns <- c(
  "Study_id",
  "Network_id",
  "Latitude",
  "Longitude",
  "Date",
  "Interaction",
  "Plant_accepted_name",
  "Plant_rank",
  "Pollinator_accepted_name",
  "Pollinator_rank"
)

missing_columns <- setdiff(
  required_columns,
  names(clough_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Ensure that Network_id contains the site code only. This is harmless if the
# year suffix was already removed by clough_data_processing.R.
clough_raw <- clough_raw %>%
  mutate(
    Site = sub("_.*", "", Network_id),
    Date = as.Date(Date),
    Year = year(Date)
  )

stopifnot(
  all(clough_raw$Study_id == "11_Clough"),
  all(!is.na(clough_raw$Site)),
  all(!is.na(clough_raw$Date))
)

# Inspect taxonomic resolution before filtering. Event_records counts rows,
# not distinct taxa.
taxonomic_resolution <- clough_raw %>%
  count(
    Plant_rank,
    Pollinator_rank,
    name = "Event_records"
  ) %>%
  arrange(
    Plant_rank,
    Pollinator_rank
  )

taxonomic_resolution

# =============================================================================
# 3. Define the sampling schedule
# =============================================================================

# The interaction-derived schedule contains only occasions with at least one
# retained interaction. Use the flower-count table when it is available.
interaction_schedule <- clough_raw %>%
  distinct(
    Site,
    Year,
    Date
  ) %>%
  arrange(
    Site,
    Date
  )

if (file.exists(flower_counts_file)) {
  flower_counts <- read.csv(
    flower_counts_file,
    stringsAsFactors = FALSE
  )
  
  flower_required_columns <- c(
    "Study_id",
    "Network_id",
    "Day",
    "Month",
    "Year"
  )
  
  flower_missing_columns <- setdiff(
    flower_required_columns,
    names(flower_counts)
  )
  
  if (length(flower_missing_columns) > 0) {
    stop(
      "Missing required flower-count columns: ",
      paste(flower_missing_columns, collapse = ", ")
    )
  }
  
  sampling_schedule <- flower_counts %>%
    filter(Study_id == "11_Clough") %>%
    transmute(
      Site = sub("_.*", "", Network_id),
      Year = as.integer(Year),
      Date = make_date(
        year = Year,
        month = Month,
        day = Day
      )
    ) %>%
    filter(
      !is.na(Site),
      !is.na(Date)
    ) %>%
    distinct() %>%
    arrange(
      Site,
      Date
    )
  
  schedule_source <- "flower_counts"
} else {
  sampling_schedule <- interaction_schedule
  schedule_source <- "interaction_records"
  
  warning(
    "Flower_counts.csv was not found. The sampling schedule excludes ",
    "occasions with no recorded interactions."
  )
}

sampling_summary <- sampling_schedule %>%
  count(
    Site,
    Year,
    name = "Sampling_rounds"
  ) %>%
  complete(
    Site,
    Year = c(2020L, 2021L),
    fill = list(Sampling_rounds = 0L)
  ) %>%
  mutate(
    Expected_rounds = if_else(
      Year == 2020L,
      3L,
      4L
    ),
    Complete = Sampling_rounds == Expected_rounds
  ) %>%
  arrange(
    Site,
    Year
  )

incomplete_sampling <- sampling_summary %>%
  filter(!Complete)

sampling_summary
incomplete_sampling

# =============================================================================
# 4. Retain valid plant--pollinator interactions
# =============================================================================

clough_interactions <- clough_raw %>%
  filter(
    !is.na(Site),
    !is.na(Plant_accepted_name),
    !is.na(Pollinator_accepted_name),
    !is.na(Interaction),
    is.finite(Interaction),
    Interaction > 0
  )

if (taxonomic_scope == "species_only") {
  clough_interactions <- clough_interactions %>%
    filter(
      Plant_rank == "SPECIES",
      Pollinator_rank == "SPECIES"
    )
}

filtering_summary <- tibble(
  Taxonomic_scope = taxonomic_scope,
  Raw_event_records = nrow(clough_raw),
  Retained_event_records = nrow(clough_interactions),
  Retained_fraction = nrow(clough_interactions) / nrow(clough_raw)
)

filtering_summary

# =============================================================================
# 5. Pool sampling rounds and years within each site
# =============================================================================

site_interactions <- clough_interactions %>%
  group_by(
    Site,
    Plant = Plant_accepted_name,
    Pollinator = Pollinator_accepted_name
  ) %>%
  summarise(
    n_visits = sum(
      Interaction,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(n_visits > 0) %>%
  arrange(
    Site,
    Plant,
    Pollinator
  )

# Weighted edge list for Rwiring. This script prepares but does not run the
# pairwise network comparisons.
interactions_rwiring <- site_interactions %>%
  transmute(
    Site,
    Origin_node = Plant,
    Destiny_node = Pollinator,
    Origin_guild = "plant",
    Destiny_guild = "pollinator",
    Strength = n_visits
  )

# =============================================================================
# 6. Dataset and network summaries
# =============================================================================

network_summary <- site_interactions %>%
  group_by(Site) %>%
  summarise(
    Plant_taxa = n_distinct(Plant),
    Pollinator_taxa = n_distinct(Pollinator),
    Unique_interactions = n(),
    Interaction_events = sum(n_visits),
    .groups = "drop"
  ) %>%
  left_join(
    sampling_schedule %>%
      group_by(Site) %>%
      summarise(
        Sampling_rounds = n(),
        Sampling_effort_min = 15 * Sampling_rounds,
        .groups = "drop"
      ),
    by = "Site"
  ) %>%
  arrange(Site)

overall_summary <- site_interactions %>%
  summarise(
    Sites = n_distinct(Site),
    Plant_taxa = n_distinct(Plant),
    Pollinator_taxa = n_distinct(Pollinator),
    Unique_interactions = n_distinct(
      Plant,
      Pollinator
    ),
    Interaction_events = sum(n_visits)
  ) %>%
  mutate(
    Sampling_rounds = nrow(sampling_schedule),
    Sampling_effort_min = 15 * Sampling_rounds,
    Sampling_effort_h = Sampling_effort_min / 60,
    Taxonomic_scope = taxonomic_scope,
    Schedule_source = schedule_source
  )

network_summary
overall_summary

# Distinct accepted taxa by recorded taxonomic rank.
plant_taxonomic_summary <- clough_interactions %>%
  distinct(
    Plant_rank,
    Plant_accepted_name
  ) %>%
  count(
    Plant_rank,
    name = "Plant_taxa"
  )

pollinator_taxonomic_summary <- clough_interactions %>%
  distinct(
    Pollinator_rank,
    Pollinator_accepted_name
  ) %>%
  count(
    Pollinator_rank,
    name = "Pollinator_taxa"
  )

plant_taxonomic_summary
pollinator_taxonomic_summary

# =============================================================================
# 7. Sample coverage of the pooled site networks
# =============================================================================

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

network_coverage <- site_interactions %>%
  group_by(Site) %>%
  summarise(
    Observed_interactions = n(),
    Interaction_events = sum(n_visits),
    Sample_coverage = get_sample_coverage(n_visits),
    .groups = "drop"
  ) %>%
  arrange(Site)

coverage_overall <- network_coverage %>%
  summarise(
    Networks = n(),
    Mean = mean(Sample_coverage, na.rm = TRUE),
    SD = sd(Sample_coverage, na.rm = TRUE),
    Median = median(Sample_coverage, na.rm = TRUE),
    Minimum = min(Sample_coverage, na.rm = TRUE),
    Maximum = max(Sample_coverage, na.rm = TRUE)
  )

network_coverage
coverage_overall

# =============================================================================
# 8. Site coordinates and pairwise geographic distances
# =============================================================================

coordinate_check <- clough_raw %>%
  group_by(Site) %>%
  summarise(
    Coordinate_pairs = n_distinct(
      paste(Latitude, Longitude)
    ),
    .groups = "drop"
  ) %>%
  filter(Coordinate_pairs != 1L)

coordinate_check
stopifnot(nrow(coordinate_check) == 0)

site_coordinates <- clough_raw %>%
  distinct(
    Site,
    Latitude,
    Longitude
  ) %>%
  arrange(Site)

stopifnot(
  all(!is.na(site_coordinates$Latitude)),
  all(!is.na(site_coordinates$Longitude)),
  nrow(site_coordinates) == n_distinct(site_coordinates$Site)
)

pairwise_distances <- site_coordinates %>%
  rename(
    Site_1 = Site,
    Latitude_1 = Latitude,
    Longitude_1 = Longitude
  ) %>%
  cross_join(
    site_coordinates %>%
      rename(
        Site_2 = Site,
        Latitude_2 = Latitude,
        Longitude_2 = Longitude
      )
  ) %>%
  filter(Site_1 < Site_2) %>%
  rowwise() %>%
  mutate(
    Distance_km = geosphere::distVincentyEllipsoid(
      c(Longitude_1, Latitude_1),
      c(Longitude_2, Latitude_2)
    ) / 1000
  ) %>%
  ungroup() %>%
  arrange(Distance_km)

distance_summary <- pairwise_distances %>%
  summarise(
    Comparisons = n(),
    Minimum_km = min(Distance_km),
    Median_km = median(Distance_km),
    Mean_km = mean(Distance_km),
    Maximum_km = max(Distance_km)
  )

pairwise_distances
distance_summary

# -------------------------------------------------------------------------
# Triangular matrix of pairwise geographic distances
# -------------------------------------------------------------------------

site_order <- sort(
  unique(
    c(
      pairwise_distances$Site_1,
      pairwise_distances$Site_2
    )
  )
)

distance_matrix <- matrix(
  NA_real_,
  nrow = length(site_order),
  ncol = length(site_order),
  dimnames = list(
    site_order,
    site_order
  )
)

for (i in seq_len(nrow(pairwise_distances))) {
  
  site_1 <- pairwise_distances$Site_1[i]
  site_2 <- pairwise_distances$Site_2[i]
  
  row_position <- match(
    site_2,
    site_order
  )
  
  column_position <- match(
    site_1,
    site_order
  )
  
  distance_matrix[
    row_position,
    column_position
  ] <- pairwise_distances$Distance_km[i]
}

distance_table <- as.data.frame(
  distance_matrix,
  check.names = FALSE
) %>%
  mutate(
    Site = rownames(distance_matrix),
    .before = 1
  )

# Format the lower triangle.
distance_table_formatted <- distance_table

for (column in site_order) {
  distance_table_formatted[[column]] <- ifelse(
    is.na(distance_table[[column]]),
    "",
    sprintf(
      "%.2f",
      distance_table[[column]]
    )
  )
}

# Mark the diagonal.
for (site in site_order) {
  distance_table_formatted[
    distance_table_formatted$Site == site,
    site
  ] <- "--"
}

distance_table_formatted

# =============================================================================
# 9. Sanity checks
# =============================================================================

n_sites <- n_distinct(site_coordinates$Site)
expected_comparisons <- choose(n_sites, 2)

stopifnot(
  n_sites == 18L,
  n_distinct(site_interactions$Site) == n_sites,
  nrow(pairwise_distances) == expected_comparisons,
  all(pairwise_distances$Site_1 != pairwise_distances$Site_2),
  all(pairwise_distances$Distance_km > 0),
  all(site_interactions$n_visits > 0)
)

# Each edge must occur only once per pooled site network.
duplicated_edges <- interactions_rwiring %>%
  count(
    Site,
    Origin_node,
    Destiny_node
  ) %>%
  filter(n > 1)

duplicated_edges
stopifnot(nrow(duplicated_edges) == 0)

# All sites in the interaction table must have coordinates and vice versa.
interaction_sites_without_coordinates <- site_interactions %>%
  distinct(Site) %>%
  anti_join(
    site_coordinates,
    by = "Site"
  )

coordinate_sites_without_interactions <- site_coordinates %>%
  distinct(Site) %>%
  anti_join(
    site_interactions %>% distinct(Site),
    by = "Site"
  )

stopifnot(
  nrow(interaction_sites_without_coordinates) == 0,
  nrow(coordinate_sites_without_interactions) == 0
)

# =============================================================================
# 10. Save processed data and descriptive outputs
# =============================================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

clough_descriptive_data <- list(
  settings = list(
    taxonomic_scope = taxonomic_scope,
    schedule_source = schedule_source
  ),
  sampling_schedule = sampling_schedule,
  sampling_summary = sampling_summary,
  incomplete_sampling = incomplete_sampling,
  clough_interactions = clough_interactions,
  site_interactions = site_interactions,
  interactions_rwiring = interactions_rwiring,
  site_coordinates = site_coordinates,
  pairwise_distances = pairwise_distances,
  network_summary = network_summary,
  overall_summary = overall_summary,
  filtering_summary = filtering_summary,
  taxonomic_resolution = taxonomic_resolution,
  plant_taxonomic_summary = plant_taxonomic_summary,
  pollinator_taxonomic_summary = pollinator_taxonomic_summary,
  network_coverage = network_coverage,
  coverage_overall = coverage_overall,
  distance_summary = distance_summary
)

saveRDS(
  clough_descriptive_data,
  file.path(
    output_dir,
    "clough_descriptive_data.rds"
  )
)

saveRDS(
  site_interactions,
  file.path(
    output_dir,
    "clough_site_interactions.rds"
  )
)

saveRDS(
  interactions_rwiring,
  file.path(
    output_dir,
    "clough_interactions_rwiring.rds"
  )
)

saveRDS(
  pairwise_distances,
  file.path(
    output_dir,
    "clough_pairwise_distances.rds"
  )
)

write.csv(
  network_summary,
  file.path(
    output_dir,
    "clough_network_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(
    output_dir,
    "clough_overall_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  network_coverage,
  file.path(
    output_dir,
    "clough_network_coverage.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_distances,
  file.path(
    output_dir,
    "clough_pairwise_distances.csv"
  ),
  row.names = FALSE
)

write.csv(
  sampling_summary,
  file.path(
    output_dir,
    "clough_sampling_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pollinator_taxonomic_summary,
  file.path(
    output_dir,
    "clough_pollinator_taxonomic_summary.csv"
  ),
  row.names = FALSE
)

# PollObs within-season interaction reorganization
#
# This script:
#   1. reads the processed weekly plant--floral visitor interactions;
#   2. constructs binary networks for Rwiring;
#   3. compares successive non-empty networks using their sampling dates;
#   4. derives link-, mechanism-, and node-level summaries; and
#   5. saves all analytical results for the plotting script.
#
# The input table must contain one row per garden, sampling period, plant, and
# floral visitor. Interaction frequencies are retained for reference, but they
# are not supplied to Rwiring because this case study uses binary networks.


# =============================================================================
# 1. Packages and file paths
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(Rwiring)

data_file <- file.path(
  "DATA",
  "PollObs",
  "processed",
  "pollobs_weekly_interactions.RDS"
)

results_file <- file.path(
  "DATA",
  "PollObs",
  "processed",
  "pollobs_rwiring_results.RDS"
)


# =============================================================================
# 2. Load and validate the processed interaction table
# =============================================================================

weekly_interactions <- readRDS(data_file)

required_columns <- c(
  "Botanical_garden",
  "ISO_year",
  "Week",
  "Date",
  "Y_day",
  "Plant",
  "Floral_visitor",
  "n_visits"
)

missing_columns <- setdiff(
  required_columns,
  names(weekly_interactions)
)

if (length(missing_columns) > 0) {
  stop(
    "The input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# This analysis currently contains a single year. Using day of year as the
# Rwiring time identifier would not be unique if several years were combined.
stopifnot(
  dplyr::n_distinct(weekly_interactions$ISO_year) == 1,
  all(weekly_interactions$ISO_year == 2023),
  inherits(weekly_interactions$Date, "Date")
)

# Confirm that the processed input has one row per interaction and network.
duplicated_pairs <- weekly_interactions %>%
  count(
    Botanical_garden,
    ISO_year,
    Week,
    Plant,
    Floral_visitor,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1)

stopifnot(
  nrow(duplicated_pairs) == 0,
  all(is.finite(weekly_interactions$n_visits)),
  all(weekly_interactions$n_visits > 0)
)

# Confirm that each garden-week corresponds to one sampling date.
network_date_check <- weekly_interactions %>%
  distinct(
    Botanical_garden,
    ISO_year,
    Week,
    Date,
    Y_day
  ) %>%
  count(
    Botanical_garden,
    ISO_year,
    Week,
    name = "n_dates"
  ) %>%
  filter(n_dates > 1)

stopifnot(nrow(network_date_check) == 0)


# =============================================================================
# 3. Construct and prepare the binary networks
# =============================================================================

# Strength is deliberately omitted: Rwiring therefore treats every observed
# interaction as present, irrespective of its recorded frequency.
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

duplicated_edges <- interactions_rwiring %>%
  count(
    Site,
    Time,
    Origin_node,
    Destiny_node,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1)

stopifnot(nrow(duplicated_edges) == 0)

prepared <- Rwiring::prepare_interactions(
  interactions_rwiring,
  directed = FALSE,
  multipartite = TRUE
)


# =============================================================================
# 4. Define comparisons between successive non-empty networks
# =============================================================================

# The day of year identifies each Rwiring network. Calendar week and date are
# retained as metadata for interpretation and plotting.
network_index <- weekly_interactions %>%
  distinct(
    Site = Botanical_garden,
    ISO_year,
    Week,
    Date,
    Time = Y_day
  ) %>%
  arrange(Site, Date)

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
    # A Date cannot retain half days; round to the closest calendar date for
    # date-based plotting while keeping the exact midpoint in Mid_day.
    Mid_date = Date_1 + round(Elapsed_days / 2),
    Transition = paste0(Week_1, "\u2192", Week_2)
  )

comparison_counts <- comparisons %>%
  count(Site, name = "n_comparisons")

expected_comparison_counts <- tibble(
  Site = c("Halle", "Jena", "Leipzig"),
  expected = c(20L, 17L, 20L)
)

comparison_count_check <- expected_comparison_counts %>%
  left_join(comparison_counts, by = "Site")

stopifnot(
  all(comparison_count_check$expected ==
        comparison_count_check$n_comparisons),
  all(comparisons$Elapsed_days > 0)
)


# =============================================================================
# 5. Quantify changes with Rwiring
# =============================================================================

comparisons <- comparisons %>%
  mutate(
    Changes = pmap(
      list(Site, Time_1, Time_2),
      function(site, time_1, time_2) {
        Rwiring::quantify_network_changes(
          prepared,
          site_1 = site,
          time_1 = time_1,
          site_2 = site,
          time_2 = time_2
        )
      }
    )
  )


# =============================================================================
# 6. Extract link- and node-level classifications
# =============================================================================

# Time_1 and Time_2 are already stored inside the nested Rwiring tables. The
# outer metadata therefore excludes them to avoid duplicated names on unnest.
comparison_metadata <- c(
  "Site",
  "ISO_year",
  "Week_1",
  "Week_2",
  "Date_1",
  "Date_2",
  "Elapsed_days",
  "Mid_day",
  "Mid_date",
  "Transition"
)

link_changes <- comparisons %>%
  transmute(
    across(all_of(comparison_metadata)),
    Links = map(Changes, "links")
  ) %>%
  unnest(Links)

node_changes <- comparisons %>%
  transmute(
    across(all_of(comparison_metadata)),
    Nodes = map(Changes, "nodes")
  ) %>%
  unnest(Nodes)

# Check that the metadata returned by Rwiring matches the requested sites.
stopifnot(
  all(link_changes$Site == link_changes$Site_1),
  all(link_changes$Site_1 == link_changes$Site_2),
  all(node_changes$Site == node_changes$Site_1),
  all(node_changes$Site_1 == node_changes$Site_2)
)

link_changes <- link_changes %>%
  select(-any_of(c("Site_1", "Site_2", "Comparison_type")))

node_changes <- node_changes %>%
  select(-any_of(c("Site_1", "Site_2", "Comparison_type")))


# =============================================================================
# 7. Relative contribution of interaction-change mechanisms
# =============================================================================

mechanism_levels <- c(
  "rewiring",
  "one_non_permanent_node",
  "two_non_permanent_nodes"
)

mechanism_labels <- c(
  rewiring = "Rewiring",
  one_non_permanent_node = "Turnover: one transient species",
  two_non_permanent_nodes = "Turnover: two transient species"
)

# Persistent links are excluded from the denominator. For mechanism m and
# comparison j, the reported proportion is (lost_m + gained_m) / (lost + gained).
mechanism_contribution <- link_changes %>%
  filter(Direction %in% c("lost", "gained")) %>%
  count(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition,
    Mechanism,
    name = "n_links"
  ) %>%
  complete(
    nesting(
      Site,
      Time_1,
      Time_2,
      Week_1,
      Week_2,
      Date_1,
      Date_2,
      Elapsed_days,
      Mid_day,
      Mid_date,
      Transition
    ),
    Mechanism = mechanism_levels,
    fill = list(n_links = 0)
  ) %>%
  group_by(Site, Time_1, Time_2) %>%
  mutate(
    Total_interaction_changes = sum(n_links),
    Proportion_interaction_change = if_else(
      Total_interaction_changes > 0,
      n_links / Total_interaction_changes,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    Mechanism_label = recode(
      as.character(Mechanism),
      !!!mechanism_labels
    ),
    Mechanism_label = factor(
      Mechanism_label,
      levels = unname(mechanism_labels)
    )
  ) %>%
  arrange(Site, Date_1, Mechanism_label)

# Independently verify the denominator and the sum of mechanism proportions.
mechanism_denominator_check <- link_changes %>%
  group_by(Site, Time_1, Time_2) %>%
  summarise(
    Expected_total = sum(Direction %in% c("lost", "gained")),
    .groups = "drop"
  ) %>%
  left_join(
    mechanism_contribution %>%
      distinct(Site, Time_1, Time_2, Total_interaction_changes),
    by = c("Site", "Time_1", "Time_2")
  )

stopifnot(
  all(mechanism_denominator_check$Expected_total ==
        mechanism_denominator_check$Total_interaction_changes)
)

mechanism_proportion_check <- mechanism_contribution %>%
  group_by(Site, Time_1, Time_2) %>%
  summarise(
    Sum_proportions = sum(Proportion_interaction_change),
    .groups = "drop"
  )

stopifnot(
  all(abs(mechanism_proportion_check$Sum_proportions - 1) < 1e-10)
)


# =============================================================================
# 8. Seven binary link-change categories
# =============================================================================

category_levels <- c(
  "Persistent",
  "Lost: rewiring",
  "Lost: one transient species",
  "Lost: two transient species",
  "Gained: rewiring",
  "Gained: one transient species",
  "Gained: two transient species"
)

link_categories <- link_changes %>%
  mutate(
    Category = case_when(
      Direction == "persistent" & Mechanism == "permanent" ~
        "Persistent",
      Direction == "lost" & Mechanism == "rewiring" ~
        "Lost: rewiring",
      Direction == "lost" & Mechanism == "one_non_permanent_node" ~
        "Lost: one transient species",
      Direction == "lost" & Mechanism == "two_non_permanent_nodes" ~
        "Lost: two transient species",
      Direction == "gained" & Mechanism == "rewiring" ~
        "Gained: rewiring",
      Direction == "gained" & Mechanism == "one_non_permanent_node" ~
        "Gained: one transient species",
      Direction == "gained" & Mechanism == "two_non_permanent_nodes" ~
        "Gained: two transient species",
      TRUE ~ NA_character_
    )
  )

stopifnot(!anyNA(link_categories$Category))

link_category_summary <- link_categories %>%
  count(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition,
    Category,
    name = "n_links"
  ) %>%
  complete(
    nesting(
      Site,
      Time_1,
      Time_2,
      Week_1,
      Week_2,
      Date_1,
      Date_2,
      Elapsed_days,
      Mid_day,
      Mid_date,
      Transition
    ),
    Category = category_levels,
    fill = list(n_links = 0)
  ) %>%
  mutate(Category = factor(Category, levels = category_levels)) %>%
  arrange(Site, Date_1, Category)

# Every link in each comparison must belong to exactly one category.
link_category_check <- link_category_summary %>%
  group_by(Site, Time_1, Time_2) %>%
  summarise(
    Classified_links = sum(n_links),
    .groups = "drop"
  ) %>%
  left_join(
    link_changes %>%
      count(Site, Time_1, Time_2, name = "Expected_links"),
    by = c("Site", "Time_1", "Time_2")
  )

stopifnot(
  all(link_category_check$Classified_links ==
        link_category_check$Expected_links)
)


# =============================================================================
# 9. Node-status summaries
# =============================================================================

node_status_levels <- c("permanent", "lost", "gained")

node_status_summary <- node_changes %>%
  count(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition,
    Guild,
    Node_status,
    name = "n_nodes"
  ) %>%
  complete(
    nesting(
      Site,
      Time_1,
      Time_2,
      Week_1,
      Week_2,
      Date_1,
      Date_2,
      Elapsed_days,
      Mid_day,
      Mid_date,
      Transition,
      Guild
    ),
    Node_status = node_status_levels,
    fill = list(n_nodes = 0)
  ) %>%
  mutate(
    Node_status = factor(Node_status, levels = node_status_levels),
    Node_status_label = recode(
      as.character(Node_status),
      permanent = "Permanent",
      lost = "Lost",
      gained = "Gained"
    ),
    Node_status_label = factor(
      Node_status_label,
      levels = c("Permanent", "Lost", "Gained")
    )
  ) %>%
  arrange(Guild, Site, Date_1, Node_status)

# Express each node status relative to the union of species in the paired
# networks. These proportions permit comparison between the two guilds.
node_status_proportions <- node_status_summary %>%
  group_by(Site, Time_1, Time_2, Guild) %>%
  mutate(
    Species_union = sum(n_nodes),
    Proportion_species = if_else(
      Species_union > 0,
      n_nodes / Species_union,
      NA_real_
    )
  ) %>%
  ungroup()

node_summary_by_garden <- node_status_proportions %>%
  group_by(Site, Guild, Node_status_label) %>%
  summarise(
    n_comparisons = n(),
    Mean_proportion = mean(Proportion_species, na.rm = TRUE),
    SD_proportion = sd(Proportion_species, na.rm = TRUE),
    Median_proportion = median(Proportion_species, na.rm = TRUE),
    Q1 = quantile(Proportion_species, 0.25, na.rm = TRUE),
    Q3 = quantile(Proportion_species, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


# =============================================================================
# 10. Net assembly balances
# =============================================================================

# Link balance: -1 indicates only losses, +1 only gains, and 0 equal numbers.
assembly_balance <- link_changes %>%
  filter(Direction %in% c("lost", "gained")) %>%
  count(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition,
    Direction,
    name = "n_links"
  ) %>%
  complete(
    nesting(
      Site,
      Time_1,
      Time_2,
      Week_1,
      Week_2,
      Date_1,
      Date_2,
      Elapsed_days,
      Mid_day,
      Mid_date,
      Transition
    ),
    Direction = c("lost", "gained"),
    fill = list(n_links = 0)
  ) %>%
  pivot_wider(
    names_from = Direction,
    values_from = n_links,
    values_fill = 0
  ) %>%
  mutate(
    Total_change = lost + gained,
    Assembly_balance = if_else(
      Total_change > 0,
      (gained - lost) / Total_change,
      NA_real_
    )
  ) %>%
  arrange(Site, Date_1)

# Node balance is calculated separately for plants and floral visitors.
node_assembly_balance <- node_status_summary %>%
  filter(Node_status %in% c("lost", "gained")) %>%
  select(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition,
    Guild,
    Node_status,
    n_nodes
  ) %>%
  pivot_wider(
    names_from = Node_status,
    values_from = n_nodes,
    values_fill = 0
  ) %>%
  mutate(
    Total_node_change = lost + gained,
    Node_assembly_balance = if_else(
      Total_node_change > 0,
      (gained - lost) / Total_node_change,
      NA_real_
    )
  ) %>%
  arrange(Guild, Site, Date_1)


# =============================================================================
# 11. Garden-level summaries and candidate matrix comparisons
# =============================================================================

mechanism_summary_by_garden <- mechanism_contribution %>%
  group_by(Site, Mechanism_label) %>%
  summarise(
    n_comparisons = n(),
    Mean_proportion = mean(
      Proportion_interaction_change,
      na.rm = TRUE
    ),
    SD_proportion = sd(
      Proportion_interaction_change,
      na.rm = TRUE
    ),
    Median_proportion = median(
      Proportion_interaction_change,
      na.rm = TRUE
    ),
    Q1 = quantile(
      Proportion_interaction_change,
      0.25,
      na.rm = TRUE
    ),
    Q3 = quantile(
      Proportion_interaction_change,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Rank potential matrix examples by the number of represented link categories,
# the number of changed links, and the resulting matrix size.
candidate_transitions <- link_category_summary %>%
  group_by(
    Site,
    Time_1,
    Time_2,
    Week_1,
    Week_2,
    Date_1,
    Date_2,
    Elapsed_days,
    Mid_day,
    Mid_date,
    Transition
  ) %>%
  summarise(
    Total_links = sum(n_links),
    Changed_links = sum(n_links[Category != "Persistent"]),
    Persistent_links = sum(n_links[Category == "Persistent"]),
    Categories_present = sum(n_links > 0),
    .groups = "drop"
  ) %>%
  left_join(
    assembly_balance %>%
      select(Site, Time_1, Time_2, Assembly_balance),
    by = c("Site", "Time_1", "Time_2")
  ) %>%
  left_join(
    node_changes %>%
      group_by(Site, Time_1, Time_2) %>%
      summarise(
        Plant_nodes = n_distinct(Node[Guild == "plant"]),
        Visitor_nodes = n_distinct(Node[Guild == "floral_visitor"]),
        Matrix_cells = Plant_nodes * Visitor_nodes,
        .groups = "drop"
      ),
    by = c("Site", "Time_1", "Time_2")
  ) %>%
  arrange(
    desc(Categories_present),
    desc(Changed_links),
    Matrix_cells
  )


# =============================================================================
# 12. Final checks and save analytical results
# =============================================================================

stopifnot(
  nrow(network_index) == 60,
  nrow(comparisons) == 57,
  all(is.finite(comparisons$Elapsed_days)),
  all(is.finite(mechanism_contribution$Proportion_interaction_change)),
  all(assembly_balance$Assembly_balance >= -1),
  all(assembly_balance$Assembly_balance <= 1),
  all(
    node_assembly_balance$Node_assembly_balance >= -1,
    na.rm = TRUE
  ),
  all(
    node_assembly_balance$Node_assembly_balance <= 1,
    na.rm = TRUE
  )
)

pollobs_rwiring_results <- list(
  interactions_rwiring = interactions_rwiring,
  network_index = network_index,
  comparisons = comparisons,
  link_changes = link_changes,
  node_changes = node_changes,
  mechanism_contribution = mechanism_contribution,
  link_categories = link_categories,
  link_category_summary = link_category_summary,
  node_status_summary = node_status_summary,
  node_status_proportions = node_status_proportions,
  node_summary_by_garden = node_summary_by_garden,
  assembly_balance = assembly_balance,
  node_assembly_balance = node_assembly_balance,
  mechanism_summary_by_garden = mechanism_summary_by_garden,
  candidate_transitions = candidate_transitions,
  validation = list(
    duplicated_pairs = duplicated_pairs,
    duplicated_edges = duplicated_edges,
    network_date_check = network_date_check,
    comparison_counts = comparison_counts,
    mechanism_denominator_check = mechanism_denominator_check,
    mechanism_proportion_check = mechanism_proportion_check,
    link_category_check = link_category_check
  )
)

dir.create(
  dirname(results_file),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  pollobs_rwiring_results,
  results_file
)

# Useful objects to inspect interactively after sourcing this script:
# summary(prepared)
# comparison_counts
# mechanism_summary_by_garden
# node_summary_by_garden
# head(candidate_transitions, 10)

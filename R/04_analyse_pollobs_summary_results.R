library(dplyr)
library(tidyr)

results <- readRDS(
  "DATA/PollObs/processed/pollobs_rwiring_results.RDS"
)

# -------------------------------------------------------------------------
# 1. Persistent, lost, and gained links
# -------------------------------------------------------------------------

direction_by_transition <- results$link_changes %>%
  count(
    Site,
    Time_1,
    Time_2,
    Direction,
    name = "n_links"
  ) %>%
  complete(
    nesting(Site, Time_1, Time_2),
    Direction = c("persistent", "lost", "gained"),
    fill = list(n_links = 0)
  ) %>%
  pivot_wider(
    names_from = Direction,
    values_from = n_links,
    values_fill = 0
  ) %>%
  mutate(
    Changed = lost + gained,
    Total = persistent + Changed,
    Proportion_persistent = persistent / Total,
    Proportion_lost = lost / Total,
    Proportion_gained = gained / Total,
    Assembly_balance = (gained - lost) / Changed
  )

direction_summary <- direction_by_transition %>%
  group_by(Site) %>%
  summarise(
    Comparisons = n(),
    
    Mean_total_links = mean(Total),
    SD_total_links = sd(Total),
    Minimum_total_links = min(Total),
    Maximum_total_links = max(Total),
    
    Mean_changed_links = mean(Changed),
    SD_changed_links = sd(Changed),
    Minimum_changed_links = min(Changed),
    Maximum_changed_links = max(Changed),
    
    Mean_persistent = mean(Proportion_persistent),
    SD_persistent = sd(Proportion_persistent),
    
    Mean_lost = mean(Proportion_lost),
    SD_lost = sd(Proportion_lost),
    
    Mean_gained = mean(Proportion_gained),
    SD_gained = sd(Proportion_gained),
    
    Mean_balance = mean(Assembly_balance),
    SD_balance = sd(Assembly_balance),
    Proportion_gain_dominated =
      mean(Assembly_balance > 0),
    Proportion_loss_dominated =
      mean(Assembly_balance < 0),
    
    .groups = "drop"
  )

direction_summary

readr::write_csv(direction_summary,"direction_summary.csv")

# -------------------------------------------------------------------------
# 2. Relative contribution of the three mechanisms
# -------------------------------------------------------------------------

mechanism_summary <- results$mechanism_contribution %>%
  group_by(Site, Mechanism_label) %>%
  summarise(
    Comparisons = n(),
    Mean = mean(Proportion_interaction_change),
    SD = sd(Proportion_interaction_change),
    Median = median(Proportion_interaction_change),
    Minimum = min(Proportion_interaction_change),
    Maximum = max(Proportion_interaction_change),
    Q1 = quantile(Proportion_interaction_change, 0.25),
    Q3 = quantile(Proportion_interaction_change, 0.75),
    .groups = "drop"
  )

mechanism_summary

# -------------------------------------------------------------------------
# 3. Species persistence and assembly balance
# -------------------------------------------------------------------------

node_status_summary <- results$node_status_proportions %>%
  group_by(
    Site,
    Guild,
    Node_status_label
  ) %>%
  summarise(
    Comparisons = n(),
    Mean = mean(Proportion_species),
    SD = sd(Proportion_species),
    Median = median(Proportion_species),
    Minimum = min(Proportion_species),
    Maximum = max(Proportion_species),
    .groups = "drop"
  )

node_balance_summary <- results$node_assembly_balance %>%
  group_by(Site, Guild) %>%
  summarise(
    Comparisons = n(),
    Mean_balance = mean(
      Node_assembly_balance,
      na.rm = TRUE
    ),
    SD_balance = sd(
      Node_assembly_balance,
      na.rm = TRUE
    ),
    Mean_absolute_balance = mean(
      abs(Node_assembly_balance),
      na.rm = TRUE
    ),
    SD_absolute_balance = sd(
      abs(Node_assembly_balance),
      na.rm = TRUE
    ),
    Proportion_gain_dominated = mean(
      Node_assembly_balance > 0,
      na.rm = TRUE
    ),
    Proportion_loss_dominated = mean(
      Node_assembly_balance < 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

node_status_summary
node_balance_summary

halle_25_26 <- results$link_category_summary %>%
  filter(
    Site == "Halle",
    Week_1 == 25,
    Week_2 == 26
  ) %>%
  select(Category, n_links)

halle_25_26

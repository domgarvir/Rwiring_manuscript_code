# PollObs temporal-lag comparisons
#
# This script compares every non-empty network with later networks from the
# same botanical garden. Temporal lag is defined as the number of sampled
# networks separating the two observations. Actual elapsed time in days is
# also retained because sampling intervals were not perfectly regular.

# =============================================================================
# 1. Packages and settings
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(Rwiring)

results_file <- file.path(
  "DATA",
  "PollObs",
  "processed",
  "pollobs_rwiring_results.RDS"
)

figure_dir <- file.path(
  "FIGURES",
  "PollObs",
  "Supplementary"
)

# Maximum number of sampling steps displayed in the analysis.
max_lag <- 4L

garden_levels <- c(
  "Halle",
  "Jena",
  "Leipzig"
)

lag_colours <- c(
  "t+1" = "#00441B",
  "t+2" = "#238B45",
  "t+3" = "#66C2A4",
  "t+4" = "#C7E9C0"
)

metric_labels <- c(
  permanent_floral_visitors =
    "Permanent floral visitors",
  permanent_plants =
    "Permanent flowering plants",
  persistent_interactions =
    "Persistent interactions"
)

theme_pollobs <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(
        fill = "white",
        colour = "black"
      ),
      strip.text = element_text(
        face = "bold",
        size = base_size
      ),
      axis.title = element_text(size = base_size),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.position = "bottom"
    )
}


# =============================================================================
# 2. Load the previously prepared networks
# =============================================================================

results <- readRDS(results_file)

required_objects <- c(
  "interactions_rwiring",
  "network_index"
)

missing_objects <- setdiff(
  required_objects,
  names(results)
)

if (length(missing_objects) > 0) {
  stop(
    "The results file is missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

interactions_rwiring <- results$interactions_rwiring
network_index <- results$network_index

stopifnot(
  all(c(
    "Site",
    "Time",
    "Origin_node",
    "Destiny_node",
    "Origin_guild",
    "Destiny_guild"
  ) %in% names(interactions_rwiring)),
  all(c(
    "Site",
    "ISO_year",
    "Week",
    "Date",
    "Time"
  ) %in% names(network_index))
)

prepared <- Rwiring::prepare_interactions(
  interactions_rwiring,
  directed = FALSE,
  multipartite = TRUE
)


# =============================================================================
# 3. Construct all within-garden comparisons
# =============================================================================

network_index_ordered <- network_index %>%
  arrange(
    Site,
    Date
  ) %>%
  group_by(Site) %>%
  mutate(
    Network_order = row_number()
  ) %>%
  ungroup()

# Join the network index to itself and retain comparisons in which Network 2
# occurs after Network 1. Lag_steps is the number of sampling steps separating
# the two networks; Elapsed_days is their actual temporal separation.
comparisons_by_lag <- network_index_ordered %>%
  rename_with(
    ~ paste0(.x, "_1"),
    -Site
  ) %>%
  inner_join(
    network_index_ordered %>%
      rename_with(
        ~ paste0(.x, "_2"),
        -Site
      ),
    by = "Site"
  ) %>%
  filter(
    Network_order_2 > Network_order_1
  ) %>%
  mutate(
    Lag_steps =
      Network_order_2 - Network_order_1,
    Elapsed_days = as.numeric(
      Date_2 - Date_1
    ),
    Mid_date = Date_1 +
      (Date_2 - Date_1) / 2,
    Transition = paste0(
      Week_1,
      "\u2192",
      Week_2
    ),
    Lag_label = paste0(
      "t+",
      Lag_steps
    )
  ) %>%
  filter(
    Lag_steps <= max_lag
  ) %>%
  arrange(
    Site,
    Network_order_1,
    Network_order_2
  )

stopifnot(
  all(comparisons_by_lag$Elapsed_days > 0),
  all(comparisons_by_lag$Lag_steps >= 1),
  all(comparisons_by_lag$Lag_steps <= max_lag)
)


# =============================================================================
# 4. Quantify network changes for every comparison
# =============================================================================

comparisons_by_lag <- comparisons_by_lag %>%
  mutate(
    Changes = pmap(
      list(
        Site,
        Time_1,
        Time_2
      ),
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
# 5. Extract node and interaction classifications
# =============================================================================

# Time_1 and Time_2 are already included in the node and link tables returned
# by Rwiring. They are therefore excluded from the external metadata to avoid
# duplicated names during unnesting.
comparison_metadata <- c(
  "Site",
  "ISO_year_1",
  "ISO_year_2",
  "Week_1",
  "Week_2",
  "Date_1",
  "Date_2",
  "Network_order_1",
  "Network_order_2",
  "Lag_steps",
  "Lag_label",
  "Elapsed_days",
  "Mid_date",
  "Transition"
)

node_changes_by_lag <- comparisons_by_lag %>%
  transmute(
    across(
      all_of(comparison_metadata)
    ),
    Nodes = map(
      Changes,
      "nodes"
    )
  ) %>%
  unnest(Nodes)

link_changes_by_lag <- comparisons_by_lag %>%
  transmute(
    across(
      all_of(comparison_metadata)
    ),
    Links = map(
      Changes,
      "links"
    )
  ) %>%
  unnest(Links)

stopifnot(
  all(
    node_changes_by_lag$Site ==
      node_changes_by_lag$Site_1
  ),
  all(
    node_changes_by_lag$Site_1 ==
      node_changes_by_lag$Site_2
  ),
  all(
    link_changes_by_lag$Site ==
      link_changes_by_lag$Site_1
  ),
  all(
    link_changes_by_lag$Site_1 ==
      link_changes_by_lag$Site_2
  )
)

node_changes_by_lag <- node_changes_by_lag %>%
  select(
    -any_of(c(
      "Site_1",
      "Site_2",
      "Comparison_type"
    ))
  )

link_changes_by_lag <- link_changes_by_lag %>%
  select(
    -any_of(c(
      "Site_1",
      "Site_2",
      "Comparison_type"
    ))
  )


# =============================================================================
# 6. Calculate proportions of permanent nodes and persistent interactions
# =============================================================================

# For each guild, the denominator is the union of species represented across
# the two networks: permanent + lost + gained species.
permanent_node_proportions <- node_changes_by_lag %>%
  count(
    Site,
    Time_1,
    Time_2,
    Lag_steps,
    Lag_label,
    Elapsed_days,
    Guild,
    Node_status,
    name = "n_species"
  ) %>%
  complete(
    nesting(
      Site,
      Time_1,
      Time_2,
      Lag_steps,
      Lag_label,
      Elapsed_days,
      Guild
    ),
    Node_status = c(
      "permanent",
      "lost",
      "gained"
    ),
    fill = list(n_species = 0)
  ) %>%
  group_by(
    Site,
    Time_1,
    Time_2,
    Lag_steps,
    Lag_label,
    Elapsed_days,
    Guild
  ) %>%
  mutate(
    Total_species = sum(n_species),
    Proportion_permanent = if_else(
      Total_species > 0,
      n_species / Total_species,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(
    Node_status == "permanent"
  ) %>%
  transmute(
    Site,
    Time_1,
    Time_2,
    Lag_steps,
    Lag_label,
    Elapsed_days,
    Metric = if_else(
      Guild == "floral_visitor",
      "permanent_floral_visitors",
      "permanent_plants"
    ),
    Proportion = Proportion_permanent
  )

# The denominator is the union of interactions represented across the paired
# networks: persistent + lost + gained interactions.
persistent_link_proportions <- link_changes_by_lag %>%
  group_by(
    Site,
    Time_1,
    Time_2,
    Lag_steps,
    Lag_label,
    Elapsed_days
  ) %>%
  summarise(
    Total_interactions = n(),
    Persistent_interactions = sum(
      Direction == "persistent"
    ),
    Proportion = Persistent_interactions /
      Total_interactions,
    .groups = "drop"
  ) %>%
  mutate(
    Metric = "persistent_interactions"
  ) %>%
  select(
    Site,
    Time_1,
    Time_2,
    Lag_steps,
    Lag_label,
    Elapsed_days,
    Metric,
    Proportion
  )

lag_proportions <- bind_rows(
  permanent_node_proportions,
  persistent_link_proportions
) %>%
  mutate(
    Site = factor(
      Site,
      levels = garden_levels
    ),
    Lag_label = factor(
      Lag_label,
      levels = paste0(
        "t+",
        seq_len(max_lag)
      )
    ),
    Metric = factor(
      Metric,
      levels = names(metric_labels),
      labels = unname(metric_labels)
    )
  )


# =============================================================================
# 7. Summarise comparisons within each garden and temporal lag
# =============================================================================

lag_summary <- lag_proportions %>%
  group_by(
    Site,
    Lag_steps,
    Lag_label,
    Metric
  ) %>%
  summarise(
    Comparisons = sum(!is.na(Proportion)),
    Mean = mean(
      Proportion,
      na.rm = TRUE
    ),
    SD = sd(
      Proportion,
      na.rm = TRUE
    ),
    Median = median(
      Proportion,
      na.rm = TRUE
    ),
    Q1 = quantile(
      Proportion,
      0.25,
      na.rm = TRUE
    ),
    Q3 = quantile(
      Proportion,
      0.75,
      na.rm = TRUE
    ),
    Mean_elapsed_days = mean(
      Elapsed_days,
      na.rm = TRUE
    ),
    SD_elapsed_days = sd(
      Elapsed_days,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Lower_SD = pmax(
      0,
      Mean - SD
    ),
    Upper_SD = pmin(
      1,
      Mean + SD
    )
  )

lag_summary


# =============================================================================
# 8. Check the number and duration of comparisons
# =============================================================================

lag_comparison_summary <- comparisons_by_lag %>%
  group_by(
    Site,
    Lag_steps,
    Lag_label
  ) %>%
  summarise(
    Comparisons = n(),
    Minimum_days = min(Elapsed_days),
    Mean_days = mean(Elapsed_days),
    SD_days = sd(Elapsed_days),
    Median_days = median(Elapsed_days),
    Maximum_days = max(Elapsed_days),
    .groups = "drop"
  )

lag_comparison_summary

# Expected number of pairs at lag k is N - k for a garden with N networks.
lag_count_check <- network_index_ordered %>%
  count(
    Site,
    name = "Networks"
  ) %>%
  crossing(
    Lag_steps = seq_len(max_lag)
  ) %>%
  mutate(
    Expected_comparisons = pmax(
      Networks - Lag_steps,
      0L
    )
  ) %>%
  left_join(
    comparisons_by_lag %>%
      count(
        Site,
        Lag_steps,
        name = "Observed_comparisons"
      ),
    by = c(
      "Site",
      "Lag_steps"
    )
  ) %>%
  mutate(
    Observed_comparisons = replace_na(
      Observed_comparisons,
      0L
    )
  )

stopifnot(
  all(
    lag_count_check$Expected_comparisons ==
      lag_count_check$Observed_comparisons
  )
)


# =============================================================================
# 9. Plot temporal-lag patterns
# =============================================================================

position_lag <- position_dodge(
  width = 0.65
)

lag_labels <- c(
  "t+1" = "t+1 (~7 days)",
  "t+2" = "t+2 (~15 days)",
  "t+3" = "t+3 (~23 days)",
  "t+4" = "t+4 (~30 days)"
)

plot_temporal_lags <- ggplot(
  lag_summary,
  aes(
    x = Site,
    y = Mean,
    colour = Lag_label,
    group = Lag_label
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower_SD,
      ymax = Upper_SD
    ),
    width = 0,
    linewidth = 0.75,
    position = position_lag
  ) +
  geom_point(
    size = 2.8,
    position = position_lag
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1
  ) +
  scale_colour_manual(
    values = lag_colours,
    breaks = names(lag_labels),
    labels = lag_labels,
    drop = FALSE,
    name = "Temporal lag"
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      0.5,
      by = 0.1
    ),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  coord_cartesian(
    ylim = c(0, 0.5)
  ) +
  labs(
    x = "Botanical garden",
    y = "Mean proportion"
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_pollobs(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    panel.spacing = unit(
      1.2,
      "lines"
    )
  )

plot_temporal_lags


# =============================================================================
# 10. Export data and figure
# =============================================================================

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    comparisons_by_lag = comparisons_by_lag,
    node_changes_by_lag = node_changes_by_lag,
    link_changes_by_lag = link_changes_by_lag,
    lag_proportions = lag_proportions,
    lag_summary = lag_summary,
    lag_comparison_summary = lag_comparison_summary
  ),
  file.path(
    "DATA",
    "PollObs",
    "processed",
    "pollobs_temporal_lag_results.RDS"
  )
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_S5_temporal_lag_comparisons.pdf"
  ),
  plot = plot_temporal_lags,
  width = 13,
  height = 4.8,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_S5_temporal_lag_comparisons.png"
  ),
  plot = plot_temporal_lags,
  width = 13,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

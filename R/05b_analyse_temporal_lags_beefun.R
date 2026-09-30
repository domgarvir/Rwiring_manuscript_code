# =============================================================================
# BeeFun interannual comparisons
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(Rwiring)


# =============================================================================
# 1. Paths and settings
# =============================================================================

data_file <- file.path(
  "DATA",
  "Ilist_byyear.csv"
)

max_lag <- 4L
selected_site <- "aznalcazar"

metric_levels <- c(
  "permanent_pollinators",
  "permanent_plants",
  "persistent_interactions",
  "persistent_frequency"
)

metric_labels <- c(
  permanent_pollinators =
    "Permanent pollinators",
  permanent_plants =
    "Permanent flowering plants",
  persistent_interactions =
    "Permanent interactions",
  persistent_frequency =
    "Permanent frequency"
)


# =============================================================================
# 2. Read and prepare the interaction list
# =============================================================================

beefun_raw <- read_csv(
  data_file,
  show_col_types = FALSE
)

stopifnot(
  all(c(
    "Site_ID",
    "Year",
    "Pollinator_gen_sp",
    "Plant_gen_sp",
    "Frequency"
  ) %in% names(beefun_raw))
)

beefun_interactions <- beefun_raw %>%
  select(
    -any_of(c("...1", "X1"))
  ) %>%
  filter(
    !is.na(Site_ID),
    !is.na(Year),
    !is.na(Pollinator_gen_sp),
    !is.na(Plant_gen_sp),
    !is.na(Frequency),
    Frequency > 0
  ) %>%
  mutate(
    Year = as.integer(Year),
    across(
      c(
        Site_ID,
        Pollinator_gen_sp,
        Plant_gen_sp
      ),
      trimws
    )
  ) %>%
  # Sum duplicated observations of the same interaction.
  group_by(
    Site_ID,
    Year,
    Plant_gen_sp,
    Pollinator_gen_sp
  ) %>%
  summarise(
    Frequency = sum(Frequency),
    .groups = "drop"
  ) %>%
  transmute(
    Site = Site_ID,
    Time = as.character(Year),
    Year,
    Origin_node = Plant_gen_sp,
    Destiny_node = Pollinator_gen_sp,
    Origin_guild = "plant",
    Destiny_guild = "floral_visitor",
    Strength = Frequency
  )

stopifnot(
  all(beefun_interactions$Strength > 0)
)


# =============================================================================
# 3. Inspect network coverage
# =============================================================================

network_index <- beefun_interactions %>%
  distinct(
    Site,
    Year,
    Time
  ) %>%
  arrange(
    Site,
    Year
  )

network_coverage <- network_index %>%
  group_by(Site) %>%
  summarise(
    First_year = min(Year),
    Last_year = max(Year),
    Networks = n(),
    Years = paste(Year, collapse = ", "),
    .groups = "drop"
  )

network_coverage

year_sequence_check <- network_index %>%
  group_by(Site) %>%
  summarise(
    Consecutive_years = all(diff(Year) == 1),
    .groups = "drop"
  )

year_sequence_check

stopifnot(
  all(year_sequence_check$Consecutive_years)
)


# =============================================================================
# 4. Prepare data for Rwiring
# =============================================================================

prepared_beefun <- Rwiring::prepare_interactions(
  beefun_interactions %>%
    select(
      Site,
      Time,
      Origin_node,
      Destiny_node,
      Strength,
      Origin_guild,
      Destiny_guild
    ),
  directed = FALSE,
  multipartite = TRUE
)


# =============================================================================
# 5. Construct all within-site comparisons at lags L1-L4
# =============================================================================

comparisons_by_lag <- network_index %>%
  rename(
    Year_1 = Year,
    Time_1 = Time
  ) %>%
  inner_join(
    network_index %>%
      rename(
        Year_2 = Year,
        Time_2 = Time
      ),
    by = "Site"
  ) %>%
  filter(
    Year_2 > Year_1
  ) %>%
  mutate(
    Lag = Year_2 - Year_1,
    Lag_label = paste0(
      "L",
      Lag
    ),
    Comparison = paste0(
      Year_1,
      "\u2192",
      Year_2
    )
  ) %>%
  filter(
    Lag <= max_lag
  ) %>%
  arrange(
    Site,
    Year_1,
    Year_2
  )

stopifnot(
  all(comparisons_by_lag$Lag %in% seq_len(max_lag))
)

comparisons_by_lag %>%
  count(
    Lag,
    Lag_label,
    name = "Comparisons"
  )


# =============================================================================
# 6. Quantify network changes using Rwiring
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
          prepared_beefun,
          site_1 = site,
          time_1 = time_1,
          site_2 = site,
          time_2 = time_2
        )
      }
    )
  )


# =============================================================================
# 7. Extract node and interaction classifications
# =============================================================================

comparison_metadata <- c(
  "Site",
  "Year_1",
  "Year_2",
  "Lag",
  "Lag_label",
  "Comparison"
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
  unnest(Nodes) %>%
  select(
    -any_of(c(
      "Site_1",
      "Site_2",
      "Comparison_type"
    ))
  )

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
  unnest(Links) %>%
  select(
    -any_of(c(
      "Site_1",
      "Site_2",
      "Comparison_type"
    ))
  )


# =============================================================================
# 8. Permanent species proportions
# =============================================================================

permanent_species <- node_changes_by_lag %>%
  mutate(
    Node_status = tolower(
      as.character(Node_status)
    )
  ) %>%
  count(
    Site,
    Year_1,
    Year_2,
    Lag,
    Lag_label,
    Comparison,
    Guild,
    Node_status,
    name = "n_species"
  ) %>%
  complete(
    nesting(
      Site,
      Year_1,
      Year_2,
      Lag,
      Lag_label,
      Comparison,
      Guild
    ),
    Node_status = c(
      "permanent",
      "lost",
      "gained"
    ),
    fill = list(
      n_species = 0
    )
  ) %>%
  group_by(
    Site,
    Year_1,
    Year_2,
    Lag,
    Lag_label,
    Comparison,
    Guild
  ) %>%
  summarise(
    Total_species = sum(n_species),
    Permanent_species = sum(
      n_species[
        Node_status == "permanent"
      ]
    ),
    Proportion = if_else(
      Total_species > 0,
      Permanent_species / Total_species,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = case_when(
      Guild == "floral_visitor" ~
        "permanent_pollinators",
      Guild == "plant" ~
        "permanent_plants",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Metric)
  )


# =============================================================================
# 9. Persistent interaction proportions
# =============================================================================

persistent_interactions <- link_changes_by_lag %>%
  mutate(
    Direction = tolower(
      as.character(Direction)
    )
  ) %>%
  group_by(
    Site,
    Year_1,
    Year_2,
    Lag,
    Lag_label,
    Comparison
  ) %>%
  summarise(
    Total_interactions = n(),
    Persistent_interactions = sum(
      Direction == "persistent",
      na.rm = TRUE
    ),
    Proportion = if_else(
      Total_interactions > 0,
      Persistent_interactions /
        Total_interactions,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = "persistent_interactions"
  )

# =============================================================================
# 10. Persistent interaction frequency
# =============================================================================

calculate_persistent_frequency <- function(
    site,
    year_1,
    year_2,
    interactions,
    tolerance = 0) {
  
  network_1 <- interactions %>%
    filter(
      Site == site,
      Year == year_1
    ) %>%
    select(
      Origin_node,
      Destiny_node,
      Strength_1 = Strength
    )
  
  network_2 <- interactions %>%
    filter(
      Site == site,
      Year == year_2
    ) %>%
    select(
      Origin_node,
      Destiny_node,
      Strength_2 = Strength
    )
  
  aligned_interactions <- full_join(
    network_1,
    network_2,
    by = c(
      "Origin_node",
      "Destiny_node"
    )
  ) %>%
    mutate(
      Strength_1 = replace_na(
        Strength_1,
        0
      ),
      Strength_2 = replace_na(
        Strength_2,
        0
      ),
      Strength_change =
        Strength_2 - Strength_1,
      
      Direction = case_when(
        Strength_1 > 0 &
          Strength_2 > 0 ~
          "persistent",
        
        Strength_1 > 0 &
          Strength_2 == 0 ~
          "lost",
        
        Strength_1 == 0 &
          Strength_2 > 0 ~
          "gained",
        
        TRUE ~ NA_character_
      ),
      
      Strength_direction = case_when(
        Direction == "persistent" &
          abs(Strength_change) <= tolerance ~
          "unchanged",
        
        Direction == "persistent" &
          Strength_change < -tolerance ~
          "weakened",
        
        Direction == "persistent" &
          Strength_change > tolerance ~
          "strengthened",
        
        Direction == "lost" ~
          "weakened",
        
        Direction == "gained" ~
          "strengthened",
        
        TRUE ~ NA_character_
      ),
      
      Frequency_contribution = case_when(
        Direction == "persistent" &
          Strength_direction == "unchanged" ~
          Strength_1,
        
        Direction == "persistent" &
          Strength_direction %in%
          c(
            "weakened",
            "strengthened"
          ) ~
          abs(Strength_change),
        
        Direction == "lost" ~
          Strength_1,
        
        Direction == "gained" ~
          Strength_2,
        
        TRUE ~ 0
      )
    )
  
  permanent_frequency <- aligned_interactions %>%
    filter(
      Direction == "persistent"
    ) %>%
    summarise(
      Value = sum(
        Frequency_contribution
      )
    ) %>%
    pull(Value)
  
  total_frequency <- sum(
    aligned_interactions$Frequency_contribution
  )
  
  if (total_frequency == 0) {
    return(
      tibble(
        Permanent_frequency = NA_real_,
        Total_frequency = NA_real_,
        Proportion = NA_real_
      )
    )
  }
  
  tibble(
    Permanent_frequency = permanent_frequency,
    Total_frequency = total_frequency,
    Proportion =
      permanent_frequency / total_frequency
  )
}

persistent_frequency <- comparisons_by_lag %>%
  mutate(
    Frequency_components = pmap(
      list(
        Site,
        Year_1,
        Year_2
      ),
      calculate_persistent_frequency,
      interactions = beefun_interactions,
      tolerance = 0
    )
  ) %>%
  unnest(
    Frequency_components
  ) %>%
  mutate(
    Metric = "persistent_frequency"
  ) %>%
  select(
    Site,
    Year_1,
    Year_2,
    Lag,
    Lag_label,
    Comparison,
    Permanent_frequency,
    Total_frequency,
    Metric,
    Proportion
  )


persistent_frequency %>%
  filter(
    Site == selected_site,
    Lag == 1
  ) %>%
  summarise(
    Comparisons = n(),
    Mean_of_proportions =
      mean(Proportion),
    SD_of_proportions =
      sd(Proportion),
    Pooled_proportion =
      sum(Permanent_frequency) /
      sum(Total_frequency)
  ) %>%
  mutate(
    across(
      c(
        Mean_of_proportions,
        SD_of_proportions,
        Pooled_proportion
      ),
      ~ 100 * .x
    )
  )

# =============================================================================
# 11. Combine the four metrics
# =============================================================================

beefun_lag_proportions <- bind_rows(
  permanent_species %>%
    select(
      Site,
      Year_1,
      Year_2,
      Lag,
      Lag_label,
      Comparison,
      Metric,
      Proportion
    ),
  persistent_interactions %>%
    select(
      Site,
      Year_1,
      Year_2,
      Lag,
      Lag_label,
      Comparison,
      Metric,
      Proportion
    ),
  persistent_frequency
) %>%
  mutate(
    Lag_label = factor(
      Lag_label,
      levels = paste0(
        "L",
        seq_len(max_lag)
      )
    ),
    Metric = factor(
      Metric,
      levels = metric_levels,
      labels = unname(metric_labels)
    )
  )

beefun_lag_summary <- beefun_lag_proportions %>%
  filter(
    Site == selected_site
  ) %>%
  group_by(
    Lag,
    Lag_label,
    Metric
  ) %>%
  summarise(
    Comparisons = sum(
      !is.na(Proportion)
    ),
    Mean = mean(
      Proportion,
      na.rm = TRUE
    ),
    SD = sd(
      Proportion,
      na.rm = TRUE
    ),
    SE = SD / sqrt(Comparisons),
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

beefun_lag_summary


# No metric should fall outside 0-1.
stopifnot(
  all(
    beefun_lag_proportions$Proportion >= 0 &
      beefun_lag_proportions$Proportion <= 1,
    na.rm = TRUE
  )
)

stopifnot(
  selected_site %in%
    unique(beefun_lag_proportions$Site)
)

beefun_lag_summary %>%
  filter(
    Lag == 1
  ) %>%
  transmute(
    Metric,
    Comparisons,
    Mean_percent = 100 * Mean,
    SD_percent = 100 * SD
  )

# Number of comparisons contributing to every metric.
beefun_lag_proportions %>%
  count(
    Metric,
    Lag_label,
    name = "Comparisons"
  )

# Means and SD in percentage units.
beefun_lag_summary %>%
  mutate(
    across(
      c(
        Mean,
        SD,
        Median,
        Q1,
        Q3
      ),
      ~ 100 * .x
    )
  )

# =============================================================================
# 12. Plot BeeFun temporal-lag patterns
# =============================================================================

beefun_lag_colours <- c(
  "L1" = "#00441B",
  "L2" = "#238B45",
  "L3" = "#41AB5D",
  "L4" = "#A1D99B"
)

beefun_lag_labels <- c(
  "L1" = "L1: t+1",
  "L2" = "L2: t+2",
  "L3" = "L3: t+3",
  "L4" = "L4: t+4"
)

theme_beefun <- function(base_size = 15) {
  theme_classic(
    base_size = base_size
  ) +
    theme(
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.6
      ),
      strip.text = element_text(
        face = "plain",
        size = base_size + 1
      ),
      axis.title = element_text(
        size = base_size
      ),
      axis.text = element_text(
        size = base_size - 1
      ),
      legend.title = element_text(
        size = base_size
      ),
      legend.text = element_text(
        size = base_size - 1
      ),
      legend.position = "bottom"
    )
}

plot_beefun_temporal_lags <- ggplot(
  beefun_lag_summary,
  aes(
    x = Lag_label,
    y = Mean,
    colour = Lag_label
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower_SD,
      ymax = Upper_SD
    ),
    width = 0,
    linewidth = 0.9
  ) +
  geom_point(
    size = 3.2
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1
  ) +
  scale_colour_manual(
    values = beefun_lag_colours,
    breaks = paste0(
      "L",
      1:4
    ),
    labels = beefun_lag_labels,
    drop = FALSE,
    name = "Time lag"
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      0.7,
      by = 0.1
    ),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      0.70
    )
  ) +
  labs(
    x = "Time lag",
    y = "Mean proportion"
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_beefun(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      size = 13
    ),
    panel.spacing = grid::unit(
      1.2,
      "lines"
    )
  )

plot_beefun_temporal_lags

# =============================================================================
# 13. Save BeeFun temporal-lag figure
# =============================================================================

beefun_output_file <- file.path(
  "DATA",
  "BeeFun",
  "beefun_temporal_lag_results.RDS"
)

beefun_figure_dir <- file.path(
  "FIGURES",
  "BeeFun"
)

dir.create(
  dirname(beefun_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  beefun_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    comparisons_by_lag = comparisons_by_lag,
    node_changes_by_lag = node_changes_by_lag,
    link_changes_by_lag = link_changes_by_lag,
    lag_proportions = beefun_lag_proportions,
    lag_summary = beefun_lag_summary,
    network_coverage = network_coverage
  ),
  beefun_output_file
)

ggsave(
  filename = file.path(
    beefun_figure_dir,
    "beefun_temporal_reorganization.pdf"
  ),
  plot = plot_beefun_temporal_lags,
  width = 15,
  height = 4.8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    beefun_figure_dir,
    "beefun_temporal_reorganization.png"
  ),
  plot = plot_beefun_temporal_lags,
  width = 15,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)


beefun_interactions %>%
  summarise(
    Origin_taxa = n_distinct(
      Origin_node
    ),
    Destiny_taxa = n_distinct(
      Destiny_node
    )
  )

beefun_lag_proportions %>%
  count(
    Metric,
    Lag_label
  )

beefun_interactions %>%
  distinct(
    Origin_guild,
    Destiny_guild
  )

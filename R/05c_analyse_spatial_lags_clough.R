# =============================================================================
# Clough spatial network comparisons
#
# TAG is used as the reference network and is compared with each of the
# remaining networks. Comparisons are subsequently grouped into four
# distance classes containing approximately the same number of networks.
# =============================================================================


# =============================================================================
# 1. Packages and settings
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(geosphere)
library(Rwiring)

data_file <- file.path(
  "DATA",
  "Clough",
  "clough_data.rds"
)

output_file <- file.path(
  "DATA",
  "Clough",
  "clough_spatial_comparisons.rds"
)

figure_dir <- file.path(
  "FIGURES",
  "Clough"
)

reference_site <- "TAG"

distance_colours <- c(
  "D1" = "#00441B",
  "D2" = "#238B45",
  "D3" = "#41AB5D",
  "D4" = "#A1D99B"
)

metric_levels <- c(
  "permanent_floral_visitors",
  "permanent_flowering_plants",
  "persistent_interactions",
  "persistent_frequency"
)

metric_labels <- c(
  permanent_floral_visitors =
    "Permanent floral visitors",
  permanent_flowering_plants =
    "Permanent flowering plants",
  persistent_interactions =
    "Permanent interactions",
  persistent_frequency =
    "Permanent frequency"
)

theme_rwiring <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(
        fill = "white",
        colour = "black"
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


# =============================================================================
# 2. Load and validate the Clough data
# =============================================================================

clough_data <- readRDS(data_file)

required_columns <- c(
  "Network_id",
  "Latitude",
  "Longitude",
  "Interaction",
  "Plant_accepted_name",
  "Pollinator_accepted_name"
)

missing_columns <- setdiff(
  required_columns,
  names(clough_data)
)

if (length(missing_columns) > 0) {
  stop(
    "The Clough data are missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

stopifnot(
  reference_site %in%
    clough_data$Network_id
)


# =============================================================================
# 3. Construct one pooled weighted network per site
# =============================================================================

clough_interactions <- clough_data %>%
  filter(
    !is.na(Network_id),
    !is.na(Plant_accepted_name),
    !is.na(Pollinator_accepted_name),
    !is.na(Interaction),
    Interaction > 0
  ) %>%
  mutate(
    across(
      c(
        Network_id,
        Plant_accepted_name,
        Pollinator_accepted_name
      ),
      trimws
    )
  ) %>%
  filter(
    Plant_accepted_name != "",
    Pollinator_accepted_name != ""
  ) %>%
  group_by(
    Network_id,
    Plant_accepted_name,
    Pollinator_accepted_name
  ) %>%
  summarise(
    Strength = sum(
      Interaction,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  transmute(
    Site = Network_id,
    Time = "pooled",
    Origin_node = Pollinator_accepted_name,
    Destiny_node = Plant_accepted_name,
    Strength,
    Origin_guild = "floral_visitor",
    Destiny_guild = "plant"
  )

stopifnot(
  all(clough_interactions$Strength > 0)
)

network_summary <- clough_interactions %>%
  group_by(Site) %>%
  summarise(
    Flowering_plants = n_distinct(
      Destiny_node
    ),
    Floral_visitors = n_distinct(
      Origin_node
    ),
    Interactions = n(),
    Interaction_events = sum(Strength),
    .groups = "drop"
  )

network_summary

# =============================================================================
# 4. Extract coordinates
# =============================================================================

coordinates <- clough_data %>%
  distinct(
    Network_id,
    Latitude,
    Longitude
  ) %>%
  filter(
    !is.na(Latitude),
    !is.na(Longitude)
  )

coordinate_check <- coordinates %>%
  count(
    Network_id,
    name = "Coordinate_pairs"
  )

stopifnot(
  all(
    coordinate_check$Coordinate_pairs == 1
  ),
  nrow(coordinates) ==
    n_distinct(clough_interactions$Site)
)

# =============================================================================
# 5. Calculate distances from TAG
# =============================================================================

reference_coordinates <- coordinates %>%
  filter(
    Network_id == reference_site
  ) %>%
  select(
    Reference_latitude = Latitude,
    Reference_longitude = Longitude
  )

stopifnot(
  nrow(reference_coordinates) == 1
)

# Compare TAG with every other site.
tag_comparisons <- coordinates %>%
  filter(
    Network_id != reference_site
  ) %>%
  transmute(
    Reference_site = reference_site,
    Comparison_site = Network_id,
    Comparison_latitude = Latitude,
    Comparison_longitude = Longitude
  ) %>%
  crossing(
    reference_coordinates
  ) %>%
  rowwise() %>%
  mutate(
    Distance_km =
      geosphere::distVincentyEllipsoid(
        c(
          Reference_longitude,
          Reference_latitude
        ),
        c(
          Comparison_longitude,
          Comparison_latitude
        )
      ) / 1000
  ) %>%
  ungroup() %>%
  arrange(Distance_km)

stopifnot(
  nrow(tag_comparisons) ==
    n_distinct(clough_interactions$Site) - 1,
  !anyDuplicated(
    tag_comparisons$Comparison_site
  ),
  all(
    tag_comparisons$Distance_km > 0
  )
)


# =============================================================================
# 6. Divide comparisons into four distance groups
# =============================================================================

# Calculate the quartile boundaries from the 17 TAG comparisons.
distance_breaks <- quantile(
  tag_comparisons$Distance_km,
  probs = seq(
    0,
    1,
    by = 0.25
  ),
  na.rm = TRUE,
  names = FALSE
)

distance_breaks

stopifnot(
  length(unique(distance_breaks)) == 5
)

# Assign every comparison to its quartile.
tag_comparisons <- tag_comparisons %>%
  mutate(
    Distance_group = cut(
      Distance_km,
      breaks = distance_breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = paste0(
        "D",
        1:4
      )
    ),
    Distance_group = factor(
      Distance_group,
      levels = paste0(
        "D",
        1:4
      )
    )
  )

# Construct legend labels using the quartile boundaries, not the observed
# minimum and maximum within each group.
distance_label_table <- tibble(
  Distance_group = factor(
    paste0(
      "D",
      1:4
    ),
    levels = paste0(
      "D",
      1:4
    )
  ),
  Distance_label = paste0(
    "D",
    1:4,
    ": (",
    sprintf(
      "%.1f",
      distance_breaks[1:4]
    ),
    ", ",
    sprintf(
      "%.1f",
      distance_breaks[2:5]
    ),
    "] km"
  )
)

tag_comparisons <- tag_comparisons %>%
  left_join(
    distance_label_table,
    by = "Distance_group"
  ) %>%
  arrange(Distance_km)

distance_group_summary <- tag_comparisons %>%
  group_by(
    Distance_group,
    Distance_label
  ) %>%
  summarise(
    Comparisons = n(),
    Minimum_observed_km = min(
      Distance_km
    ),
    Mean_km = mean(
      Distance_km
    ),
    SD_km = sd(
      Distance_km
    ),
    Maximum_observed_km = max(
      Distance_km
    ),
    .groups = "drop"
  )

distance_group_summary

tag_comparisons %>%
  select(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label
  )

# =============================================================================
# 7. Prepare the networks for Rwiring
# =============================================================================

prepared_clough <- Rwiring::prepare_interactions(
  clough_interactions,
  directed = FALSE,
  multipartite = TRUE
)


# =============================================================================
# 8. Compare TAG with every other network
# =============================================================================

tag_comparisons <- tag_comparisons %>%
  mutate(
    Changes = map(
      Comparison_site,
      function(comparison_site) {
        Rwiring::quantify_network_changes(
          prepared_clough,
          site_1 = reference_site,
          time_1 = "pooled",
          site_2 = comparison_site,
          time_2 = "pooled"
        )
      }
    )
  )


# =============================================================================
# 9. Extract node and link classifications
# =============================================================================

tag_node_changes <- tag_comparisons %>%
  transmute(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label,
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

tag_link_changes <- tag_comparisons %>%
  transmute(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label,
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
# 10. Permanent node proportions
# =============================================================================

calculate_permanent_species <- function(
    comparison_site,
    interactions,
    reference_site = "TAG") {
  
  reference_data <- interactions %>%
    filter(
      Site == reference_site
    )
  
  comparison_data <- interactions %>%
    filter(
      Site == comparison_site
    )
  
  pollinators_1 <- unique(
    reference_data$Origin_node
  )
  
  pollinators_2 <- unique(
    comparison_data$Origin_node
  )
  
  plants_1 <- unique(
    reference_data$Destiny_node
  )
  
  plants_2 <- unique(
    comparison_data$Destiny_node
  )
  
  tibble(
    Metric = c(
      "permanent_floral_visitors",
      "permanent_flowering_plants"
    ),
    Proportion = c(
      length(
        intersect(
          pollinators_1,
          pollinators_2
        )
      ) /
        length(
          union(
            pollinators_1,
            pollinators_2
          )
        ),
      
      length(
        intersect(
          plants_1,
          plants_2
        )
      ) /
        length(
          union(
            plants_1,
            plants_2
          )
        )
    )
  )
}

permanent_node_proportions <- tag_comparisons %>%
  mutate(
    Species_metrics = map(
      Comparison_site,
      calculate_permanent_species,
      interactions = clough_interactions,
      reference_site = reference_site
    )
  ) %>%
  select(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label,
    Species_metrics
  ) %>%
  unnest(Species_metrics)


# =============================================================================
# 11. Persistent interaction proportions
# =============================================================================

persistent_interaction_proportions <- tag_link_changes %>%
  mutate(
    Direction = tolower(
      as.character(Direction)
    )
  ) %>%
  group_by(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label
  ) %>%
  summarise(
    Total_interactions = n(),
    Persistent_interactions = sum(
      Direction == "persistent",
      na.rm = TRUE
    ),
    Proportion =
      Persistent_interactions /
      Total_interactions,
    .groups = "drop"
  ) %>%
  mutate(
    Metric = "persistent_interactions"
  )

# =============================================================================
# 12. Persistent interaction frequency
# =============================================================================

calculate_persistent_frequency <- function(
    comparison_site,
    interactions,
    reference_site = "TAG",
    tolerance = 0) {
  
  reference_network <- interactions %>%
    filter(
      Site == reference_site
    ) %>%
    select(
      Origin_node,
      Destiny_node,
      Strength_1 = Strength
    )
  
  comparison_network <- interactions %>%
    filter(
      Site == comparison_site
    ) %>%
    select(
      Origin_node,
      Destiny_node,
      Strength_2 = Strength
    )
  
  aligned_interactions <- full_join(
    reference_network,
    comparison_network,
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
      
      # Frequency contribution defined in rewiring_freq_py().
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
    return(NA_real_)
  }
  
  permanent_frequency / total_frequency
}

persistent_frequency_proportions <- tag_comparisons %>%
  mutate(
    Proportion = map_dbl(
      Comparison_site,
      calculate_persistent_frequency,
      interactions = clough_interactions,
      reference_site = reference_site,
      tolerance = 0
    ),
    Metric = "persistent_frequency"
  ) %>%
  select(
    Reference_site,
    Comparison_site,
    Distance_km,
    Distance_group,
    Distance_label,
    Metric,
    Proportion
  )


# =============================================================================
# 13. Combine and summarise the four metrics
# =============================================================================

clough_spatial_proportions <- bind_rows(
  permanent_node_proportions %>%
    select(
      Reference_site,
      Comparison_site,
      Distance_km,
      Distance_group,
      Distance_label,
      Metric,
      Proportion
    ),
  persistent_interaction_proportions %>%
    select(
      Reference_site,
      Comparison_site,
      Distance_km,
      Distance_group,
      Distance_label,
      Metric,
      Proportion
    ),
  persistent_frequency_proportions
) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = metric_levels,
      labels = unname(metric_labels)
    )
  )

stopifnot(
  all(
    clough_spatial_proportions$Proportion >= 0 &
      clough_spatial_proportions$Proportion <= 1,
    na.rm = TRUE
  )
)

clough_spatial_summary <- clough_spatial_proportions %>%
  group_by(
    Distance_group,
    Distance_label,
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

clough_spatial_summary

# =============================================================================
# 14. Plot spatial-distance patterns
# =============================================================================

distance_legend_labels <- setNames(
  distance_label_table$Distance_label,
  as.character(
    distance_label_table$Distance_group
  )
)

names(distance_legend_labels) <- paste0(
  "D",
  1:4
)

plot_clough_spatial <- ggplot(
  clough_spatial_summary,
  aes(
    x = Distance_group,
    y = Mean,
    colour = Distance_group
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
    values = distance_colours,
    breaks = paste0(
      "D",
      1:4
    ),
    labels = distance_legend_labels,
    drop = FALSE,
    name = "Spatial distance"
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
      mult = c(
        0.02,
        0.06
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      0.45
    )
  ) +
  labs(
    x = "Spatial distance",
    y = "Mean proportion"
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_rwiring(
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

plot_clough_spatial

# =============================================================================
# 15. Save results and figure
# =============================================================================

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    clough_interactions =
      clough_interactions,
    network_summary =
      network_summary,
    tag_comparisons =
      tag_comparisons,
    distance_group_summary =
      distance_group_summary,
    tag_node_changes =
      tag_node_changes,
    tag_link_changes =
      tag_link_changes,
    spatial_proportions =
      clough_spatial_proportions,
    spatial_summary =
      clough_spatial_summary
  ),
  output_file
)

ggsave(
  filename = file.path(
    figure_dir,
    "clough_spatial_reorganization.pdf"
  ),
  plot = plot_clough_spatial,
  width = 15,
  height = 4.8,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "clough_spatial_reorganization.png"
  ),
  plot = plot_clough_spatial,
  width = 15,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)



tag_comparisons %>%
  count(
    Distance_group,
    Distance_label
  )

round(
  distance_breaks,
  1
)

clough_interactions %>%
  summarise(
    Floral_visitors = n_distinct(
      Origin_node
    ),
    Flowering_plants = n_distinct(
      Destiny_node
    )
  )

clough_interactions %>%
  summarise(
    Origin_taxa = n_distinct(
      Origin_node
    ),
    Destiny_taxa = n_distinct(
      Destiny_node
    )
  )

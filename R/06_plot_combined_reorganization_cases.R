# =============================================================================
# Combined reorganization figure: within-season, interannual, and spatial cases
#
# This script uses only summaries saved by the three preceding analysis scripts:
#   05a_analyse_temporal_lags_pollobs.R
#   05b_analyse_temporal_lags_beefun.R
#   05c_analyse_spatial_lags_clough.R
#
# It does not repeat network preparation or Rwiring calculations.
# =============================================================================


# =============================================================================
# 1. Packages, paths, and graphical settings
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

pollobs_file <- file.path(
  "DATA",
  "PollObs",
  "processed",
  "pollobs_temporal_lag_results.RDS"
)

beefun_file <- file.path(
  "DATA",
  "BeeFun",
  "beefun_temporal_lag_results.RDS"
)

clough_file <- file.path(
  "DATA",
  "Clough",
  "clough_spatial_comparisons.rds"
)

figure_dir <- file.path(
  "FIGURES",
  "Combined"
)

output_name <- "reorganization_across_temporal_and_spatial_lags"

lag_colours_pollobs <- c(
  "L1" = "#00441B",
  "L2" = "#238B45",
  "L3" = "#66C2A4",
  "L4" = "#C7E9C0"
)

lag_colours_common <- c(
  "L1" = "#00441B",
  "L2" = "#238B45",
  "L3" = "#41AB5D",
  "L4" = "#A1D99B"
)

metric_order <- c(
  "Permanent floral visitors",
  "Permanent flowering plants",
  "Permanent interactions",
  "Permanent frequency"
)

theme_combined <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(
        fill = "white",
        colour = "grey35",
        linewidth = 0.65
      ),
      strip.text = element_text(
        face = "plain",
        size = 16
      ),
      axis.title.x = element_text(
        size = 16,
        margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        face = "plain",
        size = 18,
        margin = margin(r = 12)
      ),
      axis.text.x = element_text(
        size = 15,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 15,
        colour = "black"
      ),
      axis.line = element_line(
        colour = "grey45",
        linewidth = 0.6
      ),
      axis.ticks = element_line(
        colour = "grey45",
        linewidth = 0.6
      ),
      axis.ticks.length = grid::unit(
        0.12,
        "cm"
      ),
      panel.spacing.x = grid::unit(
        0.75,
        "lines"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(
        size = 15
      ),
      legend.text = element_text(
        size = 14
      ),
      legend.key.width = grid::unit(
        0.72,
        "cm"
      ),
      legend.key.height = grid::unit(
        0.50,
        "cm"
      ),
      legend.spacing.x = grid::unit(
        0.20,
        "cm"
      ),
      plot.margin = margin(
        6,
        6,
        6,
        6
      )
    )
}

# =============================================================================
# 2. Load and validate the saved summaries
# =============================================================================

input_files <- c(
  pollobs_file,
  beefun_file,
  clough_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Run scripts 05a, 05b, and 05c before creating the combined figure. ",
    "Missing files: ",
    paste(missing_files, collapse = ", ")
  )
}

pollobs_results <- readRDS(pollobs_file)
beefun_results <- readRDS(beefun_file)
clough_results <- readRDS(clough_file)

stopifnot(
  "lag_summary" %in% names(pollobs_results),
  "lag_summary" %in% names(beefun_results),
  "spatial_summary" %in% names(clough_results)
)


# =============================================================================
# 3. Harmonise the three summary tables
# =============================================================================

# Case 1 uses Halle only and contains three unweighted metrics.
case_1_data <- pollobs_results$lag_summary %>%
  filter(Site == "Halle") %>%
  mutate(
    Lag_axis = factor(
      paste0("L", Lag_steps),
      levels = paste0("L", 1:4)
    ),
    Metric = recode(
      as.character(Metric),
      "Permanent pollinators" =
        "Permanent floral visitors",
      "Persistent interactions" =
        "Permanent interactions"
    ),
    Metric = factor(
      Metric,
      levels = metric_order[1:3]
    )
  )

# Case 2 contains four temporal metrics pooled across BeeFun sites.
case_2_data <- beefun_results$lag_summary %>%
  mutate(
    Lag_axis = factor(
      as.character(Lag_label),
      levels = paste0("L", 1:4)
    ),
    Metric = recode(
      as.character(Metric),
      "Permanent pollinators" =
        "Permanent floral visitors"
    ),
    Metric = factor(
      Metric,
      levels = metric_order
    )
  )

# Case 3 contains the four spatial metrics from TAG comparisons.
case_3_data <- clough_results$spatial_summary %>%
  mutate(
    Distance_axis = factor(
      as.character(Distance_group),
      levels = paste0("D", 1:4)
    ),
    Metric = factor(
      as.character(Metric),
      levels = metric_order
    )
  )

stopifnot(
  nrow(case_1_data) == 12L,
  nrow(case_2_data) == 16L,
  nrow(case_3_data) == 16L,
  !anyNA(case_1_data$Metric),
  !anyNA(case_2_data$Metric),
  !anyNA(case_3_data$Metric)
)


# =============================================================================
# 4. Legend labels
# =============================================================================

case_1_legend <- c(
  "L1" = "L1: t+1 (~7 days)",
  "L2" = "L2: t+2 (~15 days)",
  "L3" = "L3: t+3 (~23 days)",
  "L4" = "L4: t+4 (~30 days)"
)

case_2_legend <- c(
  "L1" = "L1: t+1",
  "L2" = "L2: t+2",
  "L3" = "L3: t+3",
  "L4" = "L4: t+4"
)

distance_label_table <- clough_results$distance_group_summary %>%
  distinct(
    Distance_group,
    Distance_label
  ) %>%
  arrange(Distance_group)

case_3_legend <- setNames(
  distance_label_table$Distance_label,
  as.character(distance_label_table$Distance_group)
)


# =============================================================================
# 5. Construct the three case-study rows
# =============================================================================

case_1_plot <- ggplot(
  case_1_data,
  aes(
    x = Lag_axis,
    y = Mean,
    colour = Lag_axis
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower_SD,
      ymax = Upper_SD
    ),
    width = 0,
    linewidth = 0.85
  ) +
  geom_point(size = 3.1) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    drop = FALSE
  ) +
  scale_colour_manual(
    values = lag_colours_pollobs,
    breaks = names(case_1_legend),
    labels = case_1_legend,
    drop = FALSE,
    name = "Time lag"
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.4, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(
    x = NULL,
    y = "CASE 1"
  ) +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 15,
        size = 6,
        linewidth = 0
      )
    )
  ) +
  theme_combined()

# Preserve an empty fourth-panel space in Case 1 so that columns align with
# the four-metric rows below.
case_1_row <- (
  case_1_plot +
    plot_spacer()
) +
  plot_layout(
    widths = c(3, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center"
  )

case_2_plot <- ggplot(
  case_2_data,
  aes(
    x = Lag_axis,
    y = Mean,
    colour = Lag_axis
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower_SD,
      ymax = Upper_SD
    ),
    width = 0,
    linewidth = 0.85
  ) +
  geom_point(size = 3.1) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    drop = FALSE
  ) +
  scale_colour_manual(
    values = lag_colours_common,
    breaks = names(case_2_legend),
    labels = case_2_legend,
    drop = FALSE,
    name = "Time lag"
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.6, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(ylim = c(0, 0.70)) +
  labs(
    x = NULL,
    y = "CASE 2"
  ) +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 15,
        size = 6,
        linewidth = 0
      )
    )
  ) +
  theme_combined()

case_3_plot <- ggplot(
  case_3_data,
  aes(
    x = Distance_axis,
    y = Mean,
    colour = Distance_axis
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower_SD,
      ymax = Upper_SD
    ),
    width = 0,
    linewidth = 0.85
  ) +
  geom_point(size = 3.1) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    drop = FALSE
  ) +
  scale_colour_manual(
    values = setNames(
      unname(lag_colours_common),
      paste0("D", 1:4)
    ),
    breaks = names(case_3_legend),
    labels = case_3_legend,
    drop = FALSE,
    name = "Spatial distance"
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.4, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(
    x = NULL,
    y = "CASE 3"
  ) +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 15,
        size = 6,
        linewidth = 0
      )
    )
  ) +
  theme_combined()


# =============================================================================
# 6. Assemble and export the complete figure
# =============================================================================

plot_reorganization_cases <-
  case_1_row /
  case_2_plot /
  case_3_plot +
  plot_layout(
    heights = c(1, 1, 1)
  )

plot_reorganization_cases

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    paste0(output_name, ".pdf")
  ),
  plot = plot_reorganization_cases,
  width = 15,
  height = 11.5,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    paste0(output_name, ".png")
  ),
  plot = plot_reorganization_cases,
  width = 15,
  height = 11.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

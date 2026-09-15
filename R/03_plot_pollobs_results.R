# PollObs within-season figures
#
# This script produces:
#   1. the three-panel main-text figure for Halle; and
#   2. supplementary figures showing all three botanical gardens.
#
# No network calculations are repeated here. Keeping analysis and plotting in
# separate scripts makes the numerical results easier to validate and reuse.


# =============================================================================
# 1. Packages, paths, and graphical settings
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

results_file <- file.path(
  "DATA",
  "PollObs",
  "processed",
  "pollobs_rwiring_results.RDS"
)

figure_dir <- file.path("FIGURES", "PollObs")
supplementary_dir <- file.path(figure_dir, "Supplementary")

selected_site <- "Halle"

category_levels <- c(
  "Persistent",
  "Lost: rewiring",
  "Lost: one transient species",
  "Lost: two transient species",
  "Gained: rewiring",
  "Gained: one transient species",
  "Gained: two transient species"
)

category_colours <- c(
  "Persistent" = "#009455",
  
  "Lost: rewiring" = "#0072B2",
  "Lost: one transient species" = "#72B1D4",
  "Lost: two transient species" = "#C6DFEE",
  
  "Gained: rewiring" = "#E69F00",
  "Gained: one transient species" = "#F1CA72",
  "Gained: two transient species" = "#F9E9C6"
)

mechanism_colours <- c(
  "Rewiring" = "#009455",
  "Turnover: one transient species" = "#E69F00",
  "Turnover: two transient species" = "#F1CA72"
)

node_status_colours <- c(
  "Permanent" = "#009455",
  "Gained" = "#E69F00",
  "Lost" = "#0072B2"
)


garden_colours <- c(
  "Halle" = "#0072B2",
  "Jena" = "#D55E00",
  "Leipzig" = "#009E73"
)

change_levels <- c(
  "Permanent / persistent",
  "Lost",
  "Gained"
)

change_colours <- c(
  "Permanent / persistent" = "#009455",
  "Lost" = "#0072B2",
  "Gained" = "#E69F00"
)

guild_labels <- c(
  floral_visitor = "Floral visitors",
  plant = "Plants"
)

# Apply a consistent manuscript theme to all plots.
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

# Save each figure in vector and high-resolution raster formats.
save_plot_both <- function(
    plot,
    filename,
    directory,
    width,
    height,
    dpi = 600) {
  ggsave(
    filename = file.path(directory, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(directory, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
}


# =============================================================================
# 2. Load and validate analytical results
# =============================================================================

results <- readRDS(results_file)

required_objects <- c(
  "comparisons",
  "link_category_summary",
  "mechanism_contribution",
  "node_status_summary",
  "assembly_balance",
  "node_assembly_balance"
)

missing_objects <- setdiff(required_objects, names(results))

if (length(missing_objects) > 0) {
  stop(
    "The results file is missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

comparisons <- results$comparisons
link_category_summary <- results$link_category_summary
mechanism_contribution <- results$mechanism_contribution
node_status_summary <- results$node_status_summary
assembly_balance <- results$assembly_balance
node_assembly_balance <- results$node_assembly_balance

# Restore the intended plotting order even if factor attributes were lost.
link_category_summary <- link_category_summary %>%
  mutate(Category = factor(Category, levels = category_levels))

node_status_summary <- node_status_summary %>%
  mutate(
    Node_status_label = factor(
      Node_status_label,
      levels = c(
        "Permanent",
        "Lost",
        "Gained"
      )
    )
  )

stopifnot(
  all(category_levels %in% levels(link_category_summary$Category)),
  all(comparisons$Elapsed_days > 0)
)


## =============================================================================
# 3. Main-text figure, panel A: link-change categories through time in Halle
# =============================================================================

transition_axis_halle <- comparisons %>%
  filter(Site == selected_site) %>%
  distinct(Mid_day, Transition) %>%
  arrange(Mid_day)

stopifnot(!anyDuplicated(transition_axis_halle$Mid_day))

main_panel_A <- ggplot(
  link_category_summary %>%
    filter(Site == selected_site),
  aes(
    x = Mid_day,
    y = n_links,
    fill = Category
  )
) +
  geom_col(width = 4.5) +
  scale_fill_manual(
    values = category_colours,
    breaks = category_levels,
    drop = FALSE,
    name = "Interaction-change category"
  ) +
  scale_x_continuous(
    breaks = transition_axis_halle$Mid_day,
    labels = transition_axis_halle$Transition
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    tag = "A",
    x = "Successive network comparison",
    y = "Number of interactions"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme_pollobs(base_size = 15) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.text.x = element_text(
      size = 12,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.margin = margin(t = 2, b = 2)
  )

# =============================================================================
# 4. Prepare data for the species and interaction summary panels
# =============================================================================

species_plot_data <- node_status_summary %>%
  mutate(
    Guild_label = recode(
      Guild,
      plant = "Flowering plants",
      floral_visitor = "Floral visitors"
    ),
    Guild_label = factor(
      Guild_label,
      levels = c(
        "Flowering plants",
        "Floral visitors"
      )
    ),
    Change_status = recode(
      as.character(Node_status_label),
      "Permanent" = "Permanent / persistent",
      "Gained" = "Gained",
      "Lost" = "Lost"
    ),
    Change_status = factor(
      Change_status,
      levels = change_levels
    ),
    Plot_value = if_else(
      Change_status == "Lost",
      -n_nodes,
      n_nodes
    )
  )

interaction_plot_data <- link_category_summary %>%
  mutate(
    Category_character = as.character(Category),
    Mechanism_group = case_when(
      Category_character == "Persistent" ~
        "Persistent",
      grepl("rewiring", Category_character) ~
        "Rewiring",
      grepl("one transient species", Category_character) ~
        "One transient\nspecies",
      grepl("two transient species", Category_character) ~
        "Two transient\nspecies",
      TRUE ~ NA_character_
    ),
    Mechanism_group = factor(
      Mechanism_group,
      levels = c(
        "Persistent",
        "Rewiring",
        "One transient\nspecies",
        "Two transient\nspecies"
      )
    ),
    Change_status = case_when(
      Category_character == "Persistent" ~
        "Permanent / persistent",
      grepl("^Gained:", Category_character) ~
        "Gained",
      grepl("^Lost:", Category_character) ~
        "Lost",
      TRUE ~ NA_character_
    ),
    Change_status = factor(
      Change_status,
      levels = change_levels
    ),
    Plot_value = if_else(
      Change_status == "Lost",
      -n_links,
      n_links
    )
  ) %>%
  select(-Category_character)

halle_species_data <- species_plot_data %>%
  filter(Site == selected_site) %>%
  mutate(
    Change_status = factor(
      Change_status,
      levels = change_levels
    )
  )

halle_species_summary <- halle_species_data %>%
  group_by(
    Guild_label,
    Change_status
  ) %>%
  summarise(
    Mean = mean(Plot_value, na.rm = TRUE),
    SD = sd(Plot_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )

halle_interaction_data <- interaction_plot_data %>%
  filter(Site == selected_site) %>%
  mutate(
    Category = factor(
      Category,
      levels = category_levels
    )
  )

halle_interaction_summary <- halle_interaction_data %>%
  group_by(
    Mechanism_group,
    Change_status,
    Category
  ) %>%
  summarise(
    Mean = mean(Plot_value, na.rm = TRUE),
    SD = sd(Plot_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )
  
# =============================================================================
# 5. Main-text figure, panel B: mean species reorganization
# =============================================================================

position_species <- position_dodge(width = 0.72)

main_panel_B <- ggplot(
  halle_species_summary,
  aes(
    x = Guild_label,
    y = Mean,
    fill = Change_status,
    group = Change_status
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_col(
    width = 0.65,
    position = position_species
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65,
    position = position_species
  ) +
  geom_point(
    data = halle_species_data,
    aes(
      x = Guild_label,
      y = Plot_value,
      group = Change_status
    ),
    inherit.aes = FALSE,
    shape = 1,
    fill = "grey85",
    colour = "grey25",
    size = 1.5,
    stroke = 0.4,
    alpha = 0.65,
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.72,
      seed = 24
    )
  ) +
  scale_fill_manual(
    values = change_colours,
    breaks = change_levels,
    labels = c(
      "Permanent / persistent" = "Permanent",
      "Lost" = "Lost",
      "Gained" = "Gained"
    ),
    drop = FALSE,
    name = "Species status"
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    tag = "B",
    title = "Species",
    x = NULL,
    y = "Mean number of species"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_pollobs(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.margin = margin(t = 2, b = 2)
  )

# =============================================================================
# 6. Main-text figure, panel C: mean interaction reorganization
# =============================================================================

position_interactions <- position_dodge(width = 0.72)

main_panel_C <- ggplot(
  halle_interaction_summary,
  aes(
    x = Mechanism_group,
    y = Mean,
    fill = Category,
    group = Change_status
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_col(
    width = 0.65,
    position = position_interactions
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65,
    position = position_interactions
  ) +
  geom_point(
    data = halle_interaction_data,
    aes(
      x = Mechanism_group,
      y = Plot_value,
      group = Change_status
    ),
    inherit.aes = FALSE,
    shape = 1,
    fill = "grey85",
    colour = "grey25",
    size = 1.5,
    stroke = 0.4,
    alpha = 0.65,
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.72,
      seed = 24
    )
  ) +
  scale_fill_manual(
    values = category_colours,
    breaks = category_levels,
    drop = FALSE,
    guide = "none"
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    tag = "C",
    title = "Interactions",
    x = NULL,
    y = "Mean number of interactions"
  ) +
  theme_pollobs(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "none"
  )

# =============================================================================
# 7. Assemble the main-text figure
# =============================================================================

main_bottom_row <- (
  main_panel_B |
    main_panel_C
) +
  plot_layout(
    widths = c(0.85, 1.35)
  )

pollobs_main_figure <- patchwork::wrap_plots(
  main_panel_A,
  main_bottom_row,
  ncol = 1,
  heights = c(1.10, 0.90)
)

# =============================================================================
# 8. Supplementary figure S1: link-change categories in all gardens
# =============================================================================

plot_link_categories_all <- ggplot(
  link_category_summary,
  aes(
    x = Mid_date,
    y = n_links,
    fill = Category
  )
) +
  geom_col(width = 4.5) +
  facet_wrap(~ Site, ncol = 1) +
  scale_fill_manual(
    values = category_colours,
    breaks = category_levels,
    drop = FALSE
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%d %b"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Comparison midpoint",
    y = "Number of interactions",
    fill = "Interaction-change category"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pollobs() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# =============================================================================
# 9. Supplementary figure S2: mechanism proportions through time
# =============================================================================

plot_mechanism_proportions_all <- ggplot(
  mechanism_contribution,
  aes(
    x = Mid_date,
    y = Proportion_interaction_change,
    fill = Mechanism_label
  )
) +
  geom_col(width = 4.5) +
  facet_wrap(~ Site, ncol = 1) +
  scale_fill_manual(values = mechanism_colours, drop = FALSE) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%d %b"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Comparison midpoint",
    y = "Proportion of interaction change",
    fill = "Mechanism"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_pollobs() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# =============================================================================
# 10. Supplementary figure S3:
#    species assembly and disassembly through time in all gardens
# =============================================================================

# Summarise permanent, lost, and gained species for every guild, garden,
# and comparison.
node_dynamics_all <- node_status_summary %>%
  select(
    Site,
    Guild,
    Week_1,
    Week_2,
    Mid_day,
    Mid_date,
    Transition,
    Node_status_label,
    n_nodes
  ) %>%
  pivot_wider(
    names_from = Node_status_label,
    values_from = n_nodes,
    values_fill = 0
  ) %>%
  mutate(
    # Total number of distinct species represented across the paired networks.
    Total = Permanent + Lost + Gained,
    Guild = factor(
      Guild,
      levels = c(
        "floral_visitor",
        "plant"
      )
    ),
    Site = factor(
      Site,
      levels = c(
        "Halle",
        "Jena",
        "Leipzig"
      )
    )
  ) %>%
  arrange(
    Guild,
    Site,
    Mid_date
  )


# Convert gained and lost species to long format. Losses are shown below zero
# to distinguish disassembly from assembly.
node_dynamics_bars <- node_dynamics_all %>%
  select(
    Site,
    Guild,
    Mid_day,
    Mid_date,
    Transition,
    Gained,
    Lost
  ) %>%
  pivot_longer(
    cols = c(Lost, Gained),
    names_to = "Species_change",
    values_to = "n_species"
  ) %>%
  mutate(
    Species_change = factor(
      Species_change,
      levels = c(
        "Lost",
        "Gained"
      )
    ),
    Plot_value = if_else(
      Species_change == "Lost",
      -n_species,
      n_species
    )
  )


# Convert total and permanent species to long format for the two lines.
node_dynamics_lines <- node_dynamics_all %>%
  select(
    Site,
    Guild,
    Mid_day,
    Mid_date,
    Transition,
    Total,
    Permanent
  ) %>%
  pivot_longer(
    cols = c(Total, Permanent),
    names_to = "Species_pool",
    values_to = "n_species"
  ) %>%
  mutate(
    Species_pool = factor(
      Species_pool,
      levels = c(
        "Total",
        "Permanent"
      )
    )
  )


plot_node_dynamics_all <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    data = node_dynamics_bars,
    aes(
      x = Mid_date,
      y = Plot_value,
      fill = Species_change
    ),
    width = 4.5
  ) +
  geom_line(
    data = node_dynamics_lines,
    aes(
      x = Mid_date,
      y = n_species,
      colour = Species_pool,
      group = interaction(
        Site,
        Guild,
        Species_pool
      )
    ),
    linewidth = 1.1
  ) +
  geom_point(
    data = node_dynamics_lines,
    aes(
      x = Mid_date,
      y = n_species,
      colour = Species_pool
    ),
    size = 2.8
  ) +
  facet_grid(
    Guild ~ Site,
    scales = "free_y",
    labeller = labeller(
      Guild = c(
        floral_visitor = "Floral visitors",
        plant = "Flowering plants"
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      Lost = node_status_colours[["Lost"]],
      Gained = node_status_colours[["Gained"]]
    ),
    breaks = c(
      "Lost",
      "Gained"
    ),
    name = "Species change"
  ) +
  scale_colour_manual(
    values = c(
      Total = "grey45",
      Permanent = node_status_colours[["Permanent"]]
    ),
    breaks = c(
      "Total",
      "Permanent"
    ),
    name = "Species represented"
  ) +
  scale_x_date(
    date_breaks = "3 weeks",
    date_labels = "%d %b",
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  labs(
    x = "Comparison midpoint",
    y = "Number of species"
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      override.aes = list(
        linewidth = 1.1,
        size = 2.8
      )
    ),
    fill = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  theme_pollobs(base_size = 15) +
  theme(
    axis.title.x = element_text(
      size = 16,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 16,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(
      size = 12
    ),
    strip.text = element_text(
      face = "bold",
      size = 15
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 13
    ),
    legend.key.width = grid::unit(
      0.65,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.50,
      "cm"
    ),
    legend.spacing.x = grid::unit(
      0.25,
      "cm"
    ),
    panel.spacing = grid::unit(
      1,
      "lines"
    )
  )

# =============================================================================
# 12. Supplementary figure S5: net species assembly balance
# =============================================================================

plot_node_assembly_balance_all <- ggplot(
  node_assembly_balance,
  aes(
    x = Mid_date,
    y = Node_assembly_balance,
    colour = Site
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_line(
    aes(group = interaction(Site, Guild)),
    linewidth = 0.6
  ) +
  geom_point(aes(size = Total_node_change)) +
  facet_grid(
    Guild ~ Site,
    labeller = labeller(Guild = guild_labels)
  ) +
  scale_colour_manual(values = garden_colours, guide = "none") +
  scale_x_date(
    date_breaks = "3 weeks",
    date_labels = "%d %b"
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  scale_size_continuous(
    name = expression("Species changes, " ~ G[S] + L[S]),
    range = c(1.5, 5)
  ) +
  labs(
    x = "Comparison midpoint",
    y = expression(
      "Net species assembly balance, " ~
        frac(G[S] - L[S], G[S] + L[S])
    )
  ) +
  guides(size = guide_legend(nrow = 1)) +
  theme_pollobs() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# =============================================================================
# 13. Supplementary figure S6:
#     species assembly and disassembly through time in Halle
# =============================================================================

# Summarise permanent, lost, and gained species for each guild and comparison.
halle_node_dynamics <- node_status_summary %>%
  filter(Site == selected_site) %>%
  select(
    Site,
    Guild,
    Week_1,
    Week_2,
    Mid_day,
    Mid_date,
    Transition,
    Node_status_label,
    n_nodes
  ) %>%
  pivot_wider(
    names_from = Node_status_label,
    values_from = n_nodes,
    values_fill = 0
  ) %>%
  mutate(
    Total = Permanent + Lost + Gained
  ) %>%
  arrange(
    Guild,
    Mid_day
  )

# Convert gained and lost species to a long format. Losses are plotted below
# zero to distinguish species assembly from species disassembly.
halle_node_bars <- halle_node_dynamics %>%
  select(
    Site,
    Guild,
    Mid_day,
    Mid_date,
    Transition,
    Gained,
    Lost
  ) %>%
  pivot_longer(
    cols = c(Gained, Lost),
    names_to = "Species_change",
    values_to = "n_species"
  ) %>%
  mutate(
    Plot_value = if_else(
      Species_change == "Lost",
      -n_species,
      n_species
    ),
    Species_change = factor(
      Species_change,
      levels = c("Gained", "Lost")
    )
  )

# Convert total and permanent species to a long format for the two lines.
halle_node_lines <- halle_node_dynamics %>%
  select(
    Site,
    Guild,
    Mid_day,
    Mid_date,
    Transition,
    Total,
    Permanent
  ) %>%
  pivot_longer(
    cols = c(Total, Permanent),
    names_to = "Species_pool",
    values_to = "n_species"
  ) %>%
  mutate(
    Species_pool = factor(
      Species_pool,
      levels = c("Total", "Permanent")
    )
  )

# Use symmetric y-axis limits to make gains and losses directly comparable.
node_axis_limit <- max(
  abs(halle_node_bars$Plot_value),
  halle_node_lines$n_species,
  na.rm = TRUE
)

node_axis_limit <- ceiling(
  node_axis_limit / 10
) * 10

plot_node_dynamics_halle <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    data = halle_node_bars,
    aes(
      x = Mid_day,
      y = Plot_value,
      fill = Species_change
    ),
    width = 4.5,
    alpha = 0.80
  ) +
  geom_line(
    data = halle_node_lines,
    aes(
      x = Mid_day,
      y = n_species,
      colour = Species_pool,
      group = Species_pool
    ),
    linewidth = 1.05
  ) +
  geom_point(
    data = halle_node_lines,
    aes(
      x = Mid_day,
      y = n_species,
      colour = Species_pool
    ),
    size = 2.8
  ) +
  facet_wrap(
    ~ Guild,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(
      Guild = guild_labels
    )
  ) +
  scale_fill_manual(
    values = c(
      Gained = node_status_colours[["Gained"]],
      Lost = node_status_colours[["Lost"]]
    ),
    breaks = c("Gained", "Lost"),
    name = "Species change"
  ) +
  scale_colour_manual(
    values = c(
      Total = "grey45",
      Permanent = "#007A4D"
    ),
    breaks = c("Total", "Permanent"),
    name = "Species pool"
  ) +
  scale_x_continuous(
    breaks = transition_axis_halle$Mid_day,
    labels = transition_axis_halle$Transition
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  labs(
    x = "Successive network comparison",
    y = "Number of species"
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(
        linewidth = 1.05,
        size = 2.8
      )
    ),
    fill = guide_legend(
      order = 2,
      nrow = 1
    )
  ) +
  theme_pollobs(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.box = "vertical",
    legend.title.position = "top"
  )

plot_node_dynamics_halle

# =============================================================================
# 15. Supplementary figure S7: species and interaction reorganization by garden
# =============================================================================

change_colours <- c(
  "Permanent / persistent" = "#009455",
  "Gained" = "#E69F00",
  "Lost" = "#0072B2"
)

change_levels <- c(
  "Permanent / persistent",
  "Lost",
  "Gained"
)

# -----------------------------------------------------------------------------
# 14A. Species reorganization
# -----------------------------------------------------------------------------

species_plot_data <- node_status_summary %>%
  mutate(
    Guild_label = recode(
      Guild,
      plant = "Flowering plants",
      floral_visitor = "Floral visitors"
    ),
    Guild_label = factor(
      Guild_label,
      levels = c(
        "Flowering plants",
        "Floral visitors"
      )
    ),
    Change_status = recode(
      as.character(Node_status_label),
      Permanent = "Permanent / persistent",
      Gained = "Gained",
      Lost = "Lost"
    ),
    Change_status = factor(
      Change_status,
      levels = change_levels
    ),
    # Losses are represented below zero.
    Plot_value = if_else(
      Change_status == "Lost",
      -n_nodes,
      n_nodes
    )
  )

species_plot_summary <- species_plot_data %>%
  group_by(
    Site,
    Guild_label,
    Change_status
  ) %>%
  summarise(
    Mean = mean(Plot_value, na.rm = TRUE),
    SD = sd(Plot_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )

position_species <- position_dodge(
  width = 0.72
)

plot_species_reorganization <- ggplot(
  species_plot_summary,
  aes(
    x = Guild_label,
    y = Mean,
    fill = Change_status,
    group = Change_status
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    width = 0.65,
    position = position_species
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65,
    position = position_species
  ) +
  geom_point(
    data = species_plot_data,
    aes(
      x = Guild_label,
      y = Plot_value,
      group = Change_status
    ),
    inherit.aes = FALSE,
    shape = 1,
    size = 1.5,
    stroke = 0.45,
    colour = "grey20",
    fill = NA,
    alpha = 0.65,
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.72,
      seed = 24
    )
  ) +
  facet_wrap(
    ~ Site,
    ncol = 1
  ) +
  scale_fill_manual(
    values = change_colours,
    breaks = change_levels,
    labels = c(
      "Permanent / persistent" = "Permanent",
      "Gained" = "Gained",
      "Lost" = "Lost"
    ),
    drop = FALSE,
    name = "Species status"
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  labs(
    title = "Species",
    x = NULL,
    y = "Mean number of species"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_pollobs(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.text.x = element_text(
      size = 11
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    legend.key.width = grid::unit(
      0.45,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.45,
      "cm"
    )
  )


# -----------------------------------------------------------------------------
# 15B. Interaction reorganization
# -----------------------------------------------------------------------------

interaction_plot_data <- link_category_summary %>%
  mutate(
    Mechanism_group = case_when(
      Category == "Persistent" ~
        "Persistent",
      grepl("rewiring", as.character(Category)) ~
        "Rewiring",
      grepl(
        "one transient species",
        as.character(Category)
      ) ~
        "One transient\nspecies",
      grepl(
        "two transient species",
        as.character(Category)
      ) ~
        "Two transient\nspecies",
      TRUE ~ NA_character_
    ),
    Mechanism_group = factor(
      Mechanism_group,
      levels = c(
        "Persistent",
        "Rewiring",
        "One transient\nspecies",
        "Two transient\nspecies"
      )
    ),
    Change_status = case_when(
      Category == "Persistent" ~
        "Permanent / persistent",
      grepl("^Gained:", as.character(Category)) ~
        "Gained",
      grepl("^Lost:", as.character(Category)) ~
        "Lost",
      TRUE ~ NA_character_
    ),
    Change_status = factor(
      Change_status,
      levels = change_levels
    ),
    Plot_value = if_else(
      Change_status == "Lost",
      -n_links,
      n_links
    )
  )

stopifnot(
  !anyNA(interaction_plot_data$Mechanism_group),
  !anyNA(interaction_plot_data$Change_status)
)

interaction_plot_summary <- interaction_plot_data %>%
  group_by(
    Site,
    Mechanism_group,
    Change_status,
    Category
  ) %>%
  summarise(
    Mean = mean(Plot_value, na.rm = TRUE),
    SD = sd(Plot_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )

position_interactions <- position_dodge(
  width = 0.72
)

plot_interaction_reorganization <- ggplot(
  interaction_plot_summary,
  aes(
    x = Mechanism_group,
    y = Mean,
    fill = Category,
    group = Change_status
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    width = 0.65,
    position = position_interactions
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65,
    position = position_interactions
  ) +
  geom_point(
    data = interaction_plot_data,
    aes(
      x = Mechanism_group,
      y = Plot_value,
      group = Change_status
    ),
    inherit.aes = FALSE,
    shape = 1,
    size = 1.5,
    stroke = 0.45,
    colour = "grey20",
    fill = NA,
    alpha = 0.65,
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.72,
      seed = 24
    )
  ) +
  facet_wrap(
    ~ Site,
    ncol = 1
  ) +
  scale_fill_manual(
    values = category_colours,
    breaks = category_levels,
    drop = FALSE,
    name = "Interaction-change category"
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  labs(
    title = "Interactions",
    x = NULL,
    y = "Mean number of interactions"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      ncol = 2,
      byrow = FALSE
    )
  ) +
  theme_pollobs(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.text.x = element_text(
      size = 11
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(
      size = 12
    ),
    legend.text = element_text(
      size = 10.5
    ),
    legend.key.width = grid::unit(
      0.45,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.45,
      "cm"
    ),
    legend.spacing.x = grid::unit(
      0.15,
      "cm"
    )
  )

# -----------------------------------------------------------------------------
# 15C. Assemble the two columns
# -----------------------------------------------------------------------------

plot_reorganization_all_gardens <- (
  plot_species_reorganization |
    plot_interaction_reorganization
) +
  plot_layout(
    widths = c(0.85, 1.35)
  )

plot_reorganization_all_gardens

# =============================================================================
# 16. Export figures
# =============================================================================

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_both(
  pollobs_main_figure,
  "pollobs_main_figure_Halle",
  directory = figure_dir,
  width = 12,
  height = 10.5,
  dpi = 600
)

save_plot_both(
  plot_link_categories_all,
  "Figure_S1_link_change_categories_all_gardens",
  directory = supplementary_dir,
  width = 13,
  height = 8
)

save_plot_both(
  plot_mechanism_proportions_all,
  "Figure_S2_mechanism_proportions_all_gardens",
  directory = supplementary_dir,
  width = 13,
  height = 8
)

save_plot_both(
  plot_node_dynamics_all,
  "Figure_S3_node_dynamics_all_gardens",
  directory = supplementary_dir,
  width = 15,
  height = 8.5,
  dpi = 600
)

save_plot_both(
  plot_link_assembly_balance_all,
  "Figure_S4_link_assembly_balance_all_gardens",
  directory = supplementary_dir,
  width = 8,
  height = 9
)

save_plot_both(
  plot_node_assembly_balance_all,
  "Figure_S5_node_assembly_balance_all_gardens",
  directory = supplementary_dir,
  width = 13,
  height = 7
)

save_plot_both(
  plot_node_dynamics_halle,
  "Figure_S6_node_dynamics_Halle",
  directory = supplementary_dir,
  width = 12,
  height = 8
)

save_plot_both(
  plot_reorganization_all_gardens,
  "Figure_S8_species_and_interaction_reorganization",
  directory = supplementary_dir,
  width = 14,
  height = 11,
  dpi = 600
)

# Useful objects to inspect after sourcing this script:
# pollobs_main_figure
# panel_A
# panel_B
# panel_C
# plot_link_categories_all
# plot_mechanism_proportions_all
# plot_node_status_all
# plot_link_assembly_balance_all
# plot_node_assembly_balance_all
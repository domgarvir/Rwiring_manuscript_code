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

# Horizontal separation between the centres of the lost and gained bars.
# A value of 0 places both bars on the same vertical axis; 0.20 produces the
# slight overlap used in the Python reference figure.
loss_gain_separation <- 0.16

# Width of the bars in panels B and C and in Supplementary Figure S4. Because
# this is slightly larger than loss_gain_separation, lost and gained bars
# overlap horizontally while remaining visually distinguishable.
reorganization_bar_width <- 0.28

# In the species panels, keep the permanent bar to the left and place the
# overlapping lost/gained pair to its right, as in the reference figure.
species_permanent_offset <- -0.27
species_change_pair_centre <- 0.13

guild_labels <- c(
  floral_visitor = "Pollinators",
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

legend_breaks <- c(
  "Persistent",
  "Lost: rewiring",
  "Lost: one transient species",
  "Lost: two transient species",
  "Spacer",
  "Gained: rewiring",
  "Gained: one transient species",
  "Gained: two transient species"
)

legend_labels <- c(
  "Persistent",
  "Lost: rewiring",
  "Lost: one transient species",
  "Lost: two transient species",
  "",
  "Gained: rewiring",
  "Gained: one transient species",
  "Gained: two transient species"
)

theme_panel_frame <- theme(
  panel.border = element_rect(
    colour = "grey40",
    fill = NA,
    linewidth = 0.55
  ),
  axis.line = element_blank()
)

# Each comparison has two vertically aligned panels: persistent links above,
# and gained (positive) / lost (negative) links below. Both use the same x axis.
plot_link_change_pair <- function(
    data,
    x_var,
    x_scale,
    site_label = NULL,
    panel_tag = NULL,
    x_title = "Weekly network comparison",
    base_size = 13) {
  
  persistent <- data %>% filter(Category == "Persistent")
  changes <- data %>%
    filter(Category != "Persistent") %>%
    mutate(signed_links = if_else(
      grepl("^Lost:", as.character(Category)), -n_links, n_links
    ))
  
  changes_for_plot <- bind_rows(
    changes,
    persistent %>%
      slice(1) %>%
      mutate(
        n_links = 0,
        signed_links = 0
      )
  )

  common <- list(
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35),
    x_scale,
    theme_pollobs(base_size),
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(2, 8, 2, 8)
    )
  )

  top <- ggplot(persistent, aes(x = .data[[x_var]], y = n_links)) +
    geom_col(width = 4.5, fill = unname(category_colours["Persistent"])) +
    common +
    scale_y_continuous(
      breaks = seq(0, 100, by = 50),
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = site_label, tag = panel_tag,
         x = NULL, y = "Persistent\ninteractions") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.margin = margin(
        t = 2,
        r = 8,
        b = 6,
        l = 8
      )
    ) +
    theme_panel_frame

  bottom <- ggplot(
    changes_for_plot,
    aes(
      x = .data[[x_var]],
      y = signed_links,
      fill = Category
    )
  ) +
    geom_col(width = 4.5) +
    common +
    scale_fill_manual(
      values = c(
        category_colours,
        "Spacer" = "transparent"
      ),
      breaks = legend_breaks,
      labels = legend_labels,
      name = "Interaction-change category",
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.06, 0.06)
      )
    ) +
    labs(
      x = x_title,
      y = "Gained / lost\ninteractions"
    ) +
    guides(
      fill = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 2,
        byrow = TRUE
      )
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      plot.margin = margin(
        t = 6,
        r = 8,
        b = 2,
        l = 8
      )
    ) +
    theme_panel_frame
  
  (top / bottom) +
    plot_layout(
      heights = c(0.8, 1.6)
    )
}

main_panel_A <- plot_link_change_pair(
  link_category_summary %>%
    filter(Site == selected_site),
  x_var = "Mid_day",
  x_scale = scale_x_continuous(
    breaks = transition_axis_halle$Mid_day,
    labels = transition_axis_halle$Transition
  ),
  site_label = "Interaction reorganization through time",
  panel_tag = "A",
  base_size = 15
)

# =============================================================================
# 4. Prepare data for the species and interaction summary panels
# =============================================================================

species_plot_data <- node_status_summary %>%
  mutate(
    Guild_label = recode(
      Guild,
      plant = "Flowering plants",
      floral_visitor = "Pollinators"
    ),
    Guild_label = factor(
      Guild_label,
      levels = c(
        "Flowering plants",
        "Pollinators"
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
    ),
    X_base = as.numeric(Guild_label),
    X_plot = X_base + case_when(
      Change_status == "Permanent / persistent" ~
        species_permanent_offset,
      Change_status == "Lost" ~
        species_change_pair_centre - loss_gain_separation / 2,
      Change_status == "Gained" ~
        species_change_pair_centre + loss_gain_separation / 2,
      TRUE ~ 0
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
    X_plot = first(X_plot),
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
    ),
    X_base = as.numeric(Mechanism_group),
    X_plot = X_base + case_when(
      Change_status == "Lost" ~
        -loss_gain_separation / 2,
      Change_status == "Gained" ~
        loss_gain_separation / 2,
      TRUE ~ 0
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
    X_plot = first(X_plot),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )
  
# =============================================================================
# 5. Main-text figure, panel B: mean species reorganization
# =============================================================================

main_panel_B <- ggplot(
  halle_species_summary,
  aes(
    x = X_plot,
    y = Mean,
    fill = Change_status
  )
) +
  geom_col(
    width = reorganization_bar_width
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65
  ) +
  geom_point(
    data = halle_species_data,
    aes(
      x = X_plot,
      y = Plot_value
    ),
    inherit.aes = FALSE,
    shape = 1,
    fill = "grey85",
    colour = "grey25",
    size = 1.5,
    stroke = 0.4,
    alpha = 0.65,
    position = position_jitter(
      width = 0.045,
      height = 0,
      seed = 24
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(halle_species_data$Guild_label)),
    labels = levels(halle_species_data$Guild_label),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_fill_manual(
    values = change_colours,
    breaks = change_levels,
    labels = c(
      "Permanent / persistent" = "Persistent",
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
  ) +
  theme_panel_frame

# =============================================================================
# 6. Main-text figure, panel C: mean interaction reorganization
# =============================================================================

main_panel_C <- ggplot(
  halle_interaction_summary,
  aes(
    x = X_plot,
    y = Mean,
    fill = Category
  )
) +
  geom_col(
    width = reorganization_bar_width
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65
  ) +
  geom_point(
    data = halle_interaction_data,
    aes(
      x = X_plot,
      y = Plot_value
    ),
    inherit.aes = FALSE,
    shape = 1,
    fill = "grey85",
    colour = "grey25",
    size = 1.5,
    stroke = 0.4,
    alpha = 0.65,
    position = position_jitter(
      width = 0.045,
      height = 0,
      seed = 24
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(halle_interaction_data$Mechanism_group)),
    labels = levels(halle_interaction_data$Mechanism_group),
    expand = expansion(mult = c(0.08, 0.08))
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
  ) +
  theme_panel_frame

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
# 8. Supplementary figure S1:
#    interaction reorganization through time in all gardens
# =============================================================================

# Set the order in which gardens will appear.
s1_site_levels <- c(
  "Halle",
  "Jena",
  "Leipzig"
)

# Retain only gardens represented in the analytical results.
s1_site_levels <- s1_site_levels[
  s1_site_levels %in%
    unique(as.character(link_category_summary$Site))
]

# Use a common temporal range so that the seasonal position of the
# comparisons is directly comparable among gardens.
s1_midday_limits <- range(
  link_category_summary$Mid_day,
  na.rm = TRUE
)

# Create one vertically divided panel for each garden. Break positions and
# labels are garden-specific because the available weekly networks differ
# among gardens.
s1_panels <- lapply(
  s1_site_levels,
  function(site_i) {
    
    site_transition_axis <- link_category_summary %>%
      filter(Site == site_i) %>%
      distinct(
        Mid_day,
        Transition
      ) %>%
      arrange(Mid_day)
    
    stopifnot(
      !anyDuplicated(
        site_transition_axis$Mid_day
      )
    )
    
    plot_link_change_pair(
      data = link_category_summary %>%
        filter(Site == site_i),
      x_var = "Mid_day",
      x_scale = scale_x_continuous(
        limits = s1_midday_limits,
        breaks = site_transition_axis$Mid_day,
        labels = site_transition_axis$Transition,
        expand = expansion(
          mult = c(0.01, 0.01)
        )
      ),
      site_label = site_i,
      panel_tag = NULL,
      x_title = "Weekly network comparison",
      base_size = 14
    )
  }
)

# Stack the three gardens and collect their identical legends into a single
# legend at the bottom of the figure.
plot_link_categories_all <- patchwork::wrap_plots(
  plotlist = s1_panels,
  ncol = 1,
  guides = "collect"
) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title.position = "top",
    legend.title = element_text(
      size = 13
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.key.width = grid::unit(
      0.50,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.45,
      "cm"
    )
  )

plot_link_categories_all
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
        floral_visitor = "Pollinators",
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
# 14. Supplementary figure S4: species and interaction reorganization by garden
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
      floral_visitor = "Pollinators"
    ),
    Guild_label = factor(
      Guild_label,
      levels = c(
        "Flowering plants",
        "Pollinators"
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
    ),
    X_base = as.numeric(Guild_label),
    X_plot = X_base + case_when(
      Change_status == "Permanent / persistent" ~
        species_permanent_offset,
      Change_status == "Lost" ~
        species_change_pair_centre - loss_gain_separation / 2,
      Change_status == "Gained" ~
        species_change_pair_centre + loss_gain_separation / 2,
      TRUE ~ 0
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
    X_plot = first(X_plot),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )

plot_species_reorganization <- ggplot(
  species_plot_summary,
  aes(
    x = X_plot,
    y = Mean,
    fill = Change_status
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    width = reorganization_bar_width
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65
  ) +
  geom_point(
    data = species_plot_data,
    aes(
      x = X_plot,
      y = Plot_value
    ),
    inherit.aes = FALSE,
    shape = 1,
    size = 1.5,
    stroke = 0.45,
    colour = "grey20",
    fill = NA,
    alpha = 0.65,
    position = position_jitter(
      width = 0.045,
      height = 0,
      seed = 24
    )
  ) +
  facet_wrap(
    ~ Site,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(species_plot_data$Guild_label)),
    labels = levels(species_plot_data$Guild_label),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_fill_manual(
    values = change_colours,
    breaks = change_levels,
    labels = c(
      "Permanent / persistent" = "Persistent",
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
# 14B. Interaction reorganization
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
    ),
    X_base = as.numeric(Mechanism_group),
    X_plot = X_base + case_when(
      Change_status == "Lost" ~
        -loss_gain_separation / 2,
      Change_status == "Gained" ~
        loss_gain_separation / 2,
      TRUE ~ 0
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
    X_plot = first(X_plot),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = Mean - SD,
    Upper = Mean + SD
  )

plot_interaction_reorganization <- ggplot(
  interaction_plot_summary,
  aes(
    x = X_plot,
    y = Mean,
    fill = Category
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    colour = "grey25"
  ) +
  geom_col(
    width = reorganization_bar_width
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.12,
    linewidth = 0.65
  ) +
  geom_point(
    data = interaction_plot_data,
    aes(
      x = X_plot,
      y = Plot_value
    ),
    inherit.aes = FALSE,
    shape = 1,
    size = 1.5,
    stroke = 0.45,
    colour = "grey20",
    fill = NA,
    alpha = 0.65,
    position = position_jitter(
      width = 0.045,
      height = 0,
      seed = 24
    )
  ) +
  facet_wrap(
    ~ Site,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(interaction_plot_data$Mechanism_group)),
    labels = levels(interaction_plot_data$Mechanism_group),
    expand = expansion(mult = c(0.08, 0.08))
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
# 14C. Assemble the two columns
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
  plot = plot_link_categories_all,
  filename = "Figure_S1_link_change_categories_all_gardens",
  directory = supplementary_dir,
  width = 14,
  height = 15,
  dpi = 600
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
  plot_reorganization_all_gardens,
  "Figure_S4_species_and_interaction_reorganization",
  directory = supplementary_dir,
  width = 14,
  height = 11,
  dpi = 600
)

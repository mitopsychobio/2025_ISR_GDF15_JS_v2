# ============================================================================
# Figure 2E: Dose Response — with individual replicate data points
# Journal requirement: show individual data points behind mean ± SEM
# ============================================================================

library(here)

# Source the original script to get all data objects
source(here("Main Figure Scripts", "Helper Scripts", "Figure_2E_Base.R"), local = TRUE)

# Re-define output directory AFTER sourcing (source script clears environment)
out_dir <- here("Results", "Figures", "Figure_2E")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# merged_df has raw replicates: Hours, value, Agent, fa_vs_GDF15, Rep, Kinase
# df_fa_all / df_gdf_all have summarized: Hours, mean, sem, Agent (with baseline)

# Prepare raw replicate data for the same agents
agents_to_plot <- c("Dim01", "Dim195", "Dim3")

# The Dim agents don't have their own hour 0 data in merged_df.
# Hour 0 comes from UT/UT2 baselines (same as the summary data).
# Use baseline_map to pull UT/UT2 raw replicates at hour 0 and relabel them.
raw_baselines <- purrr::map_dfr(agents_to_plot, function(agent_name) {
  base_agent <- baseline_map[agent_name]
  if (!is.na(base_agent)) {
    merged_df %>%
      filter(Agent == base_agent, Hours == 0) %>%
      mutate(Agent = agent_name)
  }
})

# Dim agents' own data (hours > 0) + relabeled baseline data (hour 0)
raw_dim <- merged_df %>%
  filter(Agent %in% agents_to_plot)

raw_all <- bind_rows(raw_baselines, raw_dim)

raw_fa <- raw_all %>%
  filter(fa_vs_GDF15 == "Factor1")

raw_gdf <- raw_all %>%
  filter(fa_vs_GDF15 == "GDF15")

# ---- Factor1 plot with individual replicates ----
p_Factor1_ind <- ggplot() +
  # Error bars
  geom_errorbar(
    data = df_fa_all,
    aes(x = Hours, ymin = ymin, ymax = ymax, color = Agent),
    width = 0.6, linewidth = 0.5
  ) +
  # Mean line
  geom_line(
    data = df_fa_all,
    aes(x = Hours, y = mean, group = Agent, color = Agent),
    linewidth = 1, alpha = 0.9
  ) +
  # Mean points — translucent fill
  geom_point(
    data = df_fa_all,
    aes(x = Hours, y = mean, color = Agent),
    size = 3, alpha = 0.3
  ) +
  # Mean points — solid outline
  geom_point(
    data = df_fa_all,
    aes(x = Hours, y = mean, color = Agent),
    shape = 1, size = 3, stroke = 1
  ) +
  # Individual replicates ON TOP so they are always visible
  geom_point(
    data = raw_fa,
    aes(x = Hours, y = value, color = Agent),
    size = 2.8, alpha = 0.4,
    position = position_identity()
  ) +
  scale_color_manual(name = "Agent", values = agent_cols_Factor1) +
  scale_x_continuous(name = "Time (hours)", breaks = seq(0, 24, 4), expand = c(0, 0)) +
  ylab("Factor1 value") +
  labs(title = "Timecourse: Factor1 — Individual Replicates Shown") +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
  theme_minimal(base_size = 19) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text  = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
  )

print(p_Factor1_ind)
ggsave(file.path(out_dir, "Figure_2E_Factor1_with_individual_points.png"),
       plot = p_Factor1_ind, width = 8, height = 5.5, dpi = 300)

# ---- GDF15 plot with individual replicates ----
p_GDF15_ind <- ggplot() +
  # Error bars
  geom_errorbar(
    data = df_gdf_all,
    aes(x = Hours, ymin = ymin, ymax = ymax, color = Agent),
    width = 0.6, linewidth = 0.5
  ) +
  # Mean line
  geom_line(
    data = df_gdf_all,
    aes(x = Hours, y = mean, group = Agent, color = Agent),
    linewidth = 1, alpha = 0.9
  ) +
  # Mean points — translucent fill
  geom_point(
    data = df_gdf_all,
    aes(x = Hours, y = mean, color = Agent),
    size = 3, alpha = 0.3
  ) +
  # Mean points — solid outline
  geom_point(
    data = df_gdf_all,
    aes(x = Hours, y = mean, color = Agent),
    shape = 1, size = 3, stroke = 1
  ) +
  # Individual replicates ON TOP so they are always visible
  geom_point(
    data = raw_gdf,
    aes(x = Hours, y = value, color = Agent),
    size = 2.8, alpha = 0.4,
    position = position_identity()
  ) +
  scale_color_manual(name = "Agent", values = agent_cols_GDF15) +
  scale_x_continuous(name = "Time (hours)", breaks = seq(0, 24, 4), expand = c(0, 0)) +
  ylab("GDF15 expression (log2)") +
  labs(title = "Timecourse: GDF15 — Individual Replicates Shown") +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
  theme_minimal(base_size = 19) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text  = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
  )

print(p_GDF15_ind)
ggsave(file.path(out_dir, "Figure_2E_GDF15_with_individual_points.png"),
       plot = p_GDF15_ind, width = 8, height = 5.5, dpi = 300)

# ============================================================================
# EXPORT SOURCE DATA CSV — individual datapoints
# ============================================================================
# Combine Factor1 and GDF15 raw replicates into one tidy CSV
# Join summary stats (mean, sem) from df_fa_all / df_gdf_all
summ_fa_join <- df_fa_all %>% select(Agent, Hours, mean, sem) %>% mutate(Measure = "Factor1")
summ_gdf_join <- df_gdf_all %>% select(Agent, Hours, mean, sem) %>% mutate(Measure = "GDF15")
summ_join <- bind_rows(summ_fa_join, summ_gdf_join)

source_data_individual <- bind_rows(
  raw_fa  %>% mutate(Measure = "Factor1"),
  raw_gdf %>% mutate(Measure = "GDF15")
) %>%
  select(Agent, Hours, Rep, Measure, value) %>%
  left_join(summ_join, by = c("Agent", "Hours", "Measure")) %>%
  arrange(Measure, Agent, Hours, Rep)

write.csv(source_data_individual,
          file.path(out_dir, "Figure_2E_source_data_individual.csv"),
          row.names = FALSE)
message("  Saved: Figure_2E_source_data_individual.csv (",
        nrow(source_data_individual), " rows)")

message("\n=== Figure 2E journal-compliant plots saved ===")
message("  ", file.path(out_dir, "Figure_2E_Factor1_with_individual_points.png"))
message("  ", file.path(out_dir, "Figure_2E_GDF15_with_individual_points.png"))

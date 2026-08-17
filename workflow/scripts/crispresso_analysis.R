## =====================================================================
## VLP Base Editor Benchmarking - Analysis & Plotting
##
##
## Expected input structure (relative to this script):
##   ../../results/crispresso_tables/modifications/*.txt      - CRISPRessoPooled modification tables
##   ../../results/crispresso_tables/substitutions/*.txt      - CRISPRessoPooled substitution tables
##   ../../results/crispresso_tables/editing_summary.tsv      - per-sample aligned-read / edit counts
## Output:
##   ../../results/crispresso_tables/plots/*.pdf                   - all figures generated below
##
## =====================================================================


## 0. Setup ===========================================================

### 0a. Load packages ----

library(tidyverse)
library(tools)
library(ggh4x)
library(patchwork)
library(ggtext)


### 0b. Constants ----

# --- Dosage ---
dosage_levels        <- c("D0", "D1", "D2", "D3", "D4", "D5")
dosage_labels         <- c(D0 = "0 µL",   D1 = "0.1 µL", D2 = "0.5 µL",
                           D3 = "1 µL",   D4 = "5 µL",   D5 = "10 µL")
dosage_labels_x_axis  <- c(D0 = "0",      D1 = "0.1",    D2 = "0.5",
                           D3 = "1",      D4 = "5",      D5 = "10")
dosage_colors         <- setNames(RColorBrewer::brewer.pal(6, "Blues"), dosage_levels)

# --- Editor panel (13 benchmarked editors) ---
A_editors   <- c("ABE8e", "ABE8e-SpRY", "ABEmax-NG", "AXBE", "AYBE")
C_editors   <- c("CBE6d", "TadCBEd", "CGBE1", "CGBE1-NG", "miniCGBE1-SpRY")
AC_editors  <- c("SPACE", "SPACE-DeltaUGI", "SPACE-NG")

# --- Genomic targets ---
targets <- c("PDCD1", "B2M", "HEK3")

# --- Plotting ---
error_bar_type <- "sem"
dodge_width    <- 0.9

# Shared per-editor color palette used by every A-/C-editor comparison
# plot (dosage-across-editor scatter/ratio/line plots). Consolidated
# here from three identical copies in the original script.
editor_palette <- c(
  # A-editors: distinct, saturated blues/teals -- not shades of one hue
  "ABE8e"           = "#08306B",  # deep navy
  "ABE8e-SpRY"      = "#1F78B4",  # strong blue
  "ABEmax-NG"       = "#33A02C",  # teal-green (deliberately distinct from the blues)
  "AXBE"            = "#00B4D8",  # bright cyan-blue
  "AYBE"            = "#6A3D9A",  # blue-violet
  # C-editors: distinct, saturated reds/warm hues
  "CBE6d"           = "#67001F",  # deep maroon
  "TadCBEd"         = "#E31A1C",  # strong red
  "CGBE1"           = "#FF7F00",  # orange
  "CGBE1-NG"        = "#FB9A99",  # salmon pink
  "miniCGBE1-SpRY"  = "#B15928"   # brown-red
)

# Fixed target facet color scheme, reused by both spacer line plots.
target_colors <- c(HEK3 = "#1b9e77", B2M = "#d95f02", PDCD1 = "#7570b3")


### 0c. Paths ----

data_dir  <- "./Data"
plots_dir <- "./plots"


### 0d. Data loading / wrangling helper functions ----

#' Load a folder of CRISPRessoPooled per-sample TSVs into a named list of
#' data frames, parsing the "Folder" column into base_editor/dosage/
#' replicate/target.
load_crispresso_folder <- function(path) {
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
  names(files) <- file_path_sans_ext(basename(files))
  
  dfs <- lapply(files, read_tsv)
  
  # Strip the CRISPRessoPooled run-name prefix from every character column
  dfs <- lapply(dfs, function(df) {
    df[] <- lapply(df, function(x) {
      if (is.character(x)) str_remove(x, ".*CRISPRessoPooled_on_") else x
    })
    df
  })
  
  map(dfs, ~ .x %>%
        separate(
          col = Folder,
          into = c("base_editor", "dosage", "replicate", "tmp1", "tmp2", "target"),
          sep = "_"
        ) %>%
        select(-tmp1, -tmp2))
}

#' Split a list of per-sample data frames into one bound data frame per
#' genomic target.
split_by_target <- function(df_list, target_names = targets) {
  set_names(target_names) %>%
    map(~ bind_rows(Filter(\(x) unique(x$target) == .x, df_list)))
}

#' Load and reshape CRISPResso nucleotide-frequency distribution tables,
#' converting raw counts (fq) to a per-sample percentage.
load_distribution_data <- function(path, suffix_pattern) {
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
  names(files) <- file_path_sans_ext(basename(files))
  
  bind_rows(lapply(files, read_tsv), .id = "source") %>%
    mutate(source = str_remove(source, suffix_pattern)) %>%
    separate(source, into = c("base_editor", "dosage", "target", "replicate"), sep = "_") %>%
    group_by(base_editor, dosage, target, replicate) %>%
    mutate(total_fq = sum(fq), percent = fq / total_fq * 100) %>%
    ungroup()
}

#' Rename per-nucleotide sgRNA columns (cols 6:ncol) from their raw base
#' letter to "<letter><position>" (e.g. "A1", "C2", ...).
rename_sgrna <- function(df) {
  nucs <- substr(names(df)[6:ncol(df)], 1, 1)
  names(df)[6:ncol(df)] <- paste0(nucs, seq_along(nucs))
  df
}

#' Reverse-complement a target's per-nucleotide columns and Modification
#' labels (used for B2M, which is targeted on the minus strand).
flip_reverse_complement <- function(df) {
  meta_cols <- 1:5
  nuc_cols  <- 6:ncol(df)
  rc <- function(x) chartr("ATCG", "TAGC", x)
  
  df <- df %>%
    mutate(Modification = case_when(
      Modification == "A" ~ "T",
      Modification == "T" ~ "A",
      Modification == "C" ~ "G",
      Modification == "G" ~ "C",
      TRUE ~ Modification
    ))
  
  nuc_data <- df[, nuc_cols]
  nuc_data <- nuc_data[, ncol(nuc_data):1]
  colnames(nuc_data) <- sapply(seq_along(colnames(nuc_data)), function(i) {
    paste0(rc(substr(colnames(nuc_data)[i], 1, 1)), i)
  })
  
  cbind(df[, meta_cols], nuc_data)
}

#' Restrict rows to the base-editor-appropriate base letter (A-editors ->
#' "A", C-editors -> "C", dual A/C editors -> either), optionally
#' excluding specific "target__position" combinations.
get_qualifying_positions <- function(df, exclude = character(0)) {
  qp <- df %>%
    filter(
      (base_editor %in% A_editors  & base_letter == "A") |
        (base_editor %in% C_editors  & base_letter == "C") |
        (base_editor %in% AC_editors & base_letter %in% c("A", "C"))
    )
  if (length(exclude) > 0) {
    qp <- qp %>% filter(!paste(target, position, sep = "__") %in% exclude)
  }
  qp
}

#' Compute grouped mean +/- error (SEM or SD) for a numeric value column.
summarise_stats <- function(df, group_cols, value_col = "value", error_type = error_bar_type) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      sd_val   = sd(.data[[value_col]], na.rm = TRUE),
      n        = n(),
      sem_val  = sd_val / sqrt(n),
      .groups  = "drop"
    ) %>%
    mutate(err = if (error_type == "sd") sd_val else sem_val)
}

#' Build a per-dosage white-to-hue color ramp for each named base hue
#' (used for the dosage-graded fill legends).
make_dosage_fill_lookup <- function(hues, levels = dosage_levels, start_index = NULL) {
  imap(hues, function(hue, key) {
    ramp <- if (is.null(start_index)) {
      colorRampPalette(c("white", hue))(length(levels) + 1)[-1]
    } else {
      start_col <- colorRampPalette(c("white", hue))(10)[start_index]
      colorRampPalette(c(start_col, hue))(length(levels) + 1)[-1]
    }
    set_names(ramp, paste0(key, "_", levels))
  }) %>% flatten_chr()
}

#' Build a simple tile-and-label legend for a set of named fill colors.
make_tile_legend <- function(items, hues) {
  legend_df <- data.frame(
    item = factor(items, levels = items),
    y = seq_along(items)
  )
  ggplot(legend_df, aes(x = 1, y = y, fill = item)) +
    geom_tile(width = 0.7, height = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = item), x = 1.5, hjust = 0, size = 4, color = "black") +
    scale_fill_manual(values = hues, guide = "none") +
    coord_cartesian(xlim = c(0.5, 3), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5))
}

#' Write a list of ggplot objects to a single multi-page PDF.
save_plots_pdf <- function(plots, filename, width, height) {
  pdf(filename, width = width, height = height)
  walk(plots, print)
  dev.off()
}


## 1. Load & tidy modification / substitution / summary data =========

### 1a. Load raw CRISPResso data ----
mods <- load_crispresso_folder(file.path(data_dir, "modifications"))
subs <- load_crispresso_folder(file.path(data_dir, "substitutions"))

# Per-sample editing summary table (aligned read counts, total
# substitutions/insertions/deletions). Loaded once here and reused by
# every section that needs it (Sections 2c, 4), instead of being
# re-read from disk in each section as in the original script.
editing_summary <- read_tsv(
  file.path(data_dir, "editing_summary.tsv"),   # adjust path if not under data_dir
  show_col_types = FALSE
)

### 1b. Format data ----

#### i. Split into per-target data frames ----
mods_by_target <- split_by_target(mods)
subs_by_target <- split_by_target(subs) %>%
  map(~ .x %>% filter(Nucleotide != "-") %>% rename(Modification = Nucleotide))

#### ii. Combine modifications & substitutions ----
target_data <- map2(mods_by_target, subs_by_target, bind_rows)
list2env(target_data, envir = environment())

rm(mods, subs, mods_by_target, subs_by_target, target_data)

#### iii. Isolate & format sgRNA sequence columns ----
target_data <- list(PDCD1 = PDCD1, B2M = B2M, HEK3 = HEK3) %>%
  map(rename_sgrna)

# B2M is targeted on the minus strand -> reverse complement
target_data$B2M <- flip_reverse_complement(target_data$B2M)

list2env(target_data, envir = environment())
rm(target_data)

#### iv. Build long-format data ----

dat_long <- bind_rows(
  B2M   %>% pivot_longer(cols = matches("^[A-Z][0-9]+$"), names_to = "position", values_to = "value"),
  HEK3  %>% pivot_longer(cols = matches("^[A-Z][0-9]+$"), names_to = "position", values_to = "value"),
  PDCD1 %>% pivot_longer(cols = matches("^[A-Z][0-9]+$"), names_to = "position", values_to = "value")
)

dat_long <- dat_long %>%
  mutate(
    dosage      = factor(dosage, levels = dosage_levels),
    base_letter = str_extract(position, "^[A-Za-z]"),
    pos_num     = as.numeric(str_extract(position, "\\d+")),
    position    = fct_reorder(position, pos_num),
    target_position = factor(
      paste(target, position, sep = "__"),
      levels = paste(target, position, sep = "__")[order(target, pos_num)] %>% unique()
    )
  )

rm(B2M, HEK3, PDCD1)


## 2. Substitution plots ==============================================

### 2a. Total substitutions ----

#### Bar plot ----

main_data <- dat_long %>% filter(Modification == "Substitutions")

adenine_editors   <- A_editors
cytosine_editors  <- C_editors

# Same threshold logic as "D5 Editing - Adenine vs. Cytosine Heatmaps": within a
# family, a target/position qualifies if ANY member editor of that family reaches
# >=1% mean editing at D5. Applied per-family so A editors only ever see
# base_letter == "A" rows and C editors only ever see base_letter == "C" rows.
# No pos_num range restriction -- the full protospacer is considered.
build_family_data <- function(editors, letter, family_name) {
  
  d5_raw <- main_data %>%
    filter(dosage == "D5", base_editor %in% editors, base_letter == letter)
  
  d5_summary <- summarise_stats(d5_raw, c("base_editor", "target", "position", "pos_num"))
  
  qualifying_tp <- d5_summary %>%
    group_by(target, position) %>%
    summarise(max_mean = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
    filter(max_mean >= 0.01) %>%
    select(target, position)
  
  main_data %>%
    filter(base_editor %in% editors, base_letter == letter) %>%
    semi_join(qualifying_tp, by = c("target", "position")) %>%
    mutate(family = family_name)
}

plot_data_raw <- bind_rows(
  build_family_data(adenine_editors,  "A", "Adenine editors"),
  build_family_data(cytosine_editors, "C", "Cytosine editors")
)

summary_data <- summarise_stats(
  plot_data_raw,
  c("base_editor", "target", "position", "pos_num", "dosage", "family")
)

make_family_plot <- function(family_name, show_x_title = TRUE) {
  
  sd  <- summary_data  %>% filter(family == family_name) %>% mutate(position = droplevels(position))
  pdr <- plot_data_raw %>% filter(family == family_name) %>% mutate(position = droplevels(position))
  
  ggplot(sd, aes(x = position, y = mean_val, fill = dosage)) +
    geom_col(position = position_dodge(width = dodge_width), color = "black", linewidth = 0.15) +
    geom_errorbar(
      aes(ymin = pmax(mean_val - err, 0), ymax = mean_val + err),
      position = position_dodge(width = dodge_width), width = 0.25
    ) +
    geom_point(
      data = pdr,
      aes(x = position, y = value, group = dosage),
      position = position_dodge(width = dodge_width),
      size = 0.7, alpha = 0.7, color = "black"
    ) +
    facet_grid(
      rows = vars(base_editor), cols = vars(target),
      scales = "free",
      space  = "free_x"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_fill_manual(values = dosage_colors, labels = dosage_labels, drop = FALSE, name = "Dosage") +
    labs(
      title = family_name,
      x = if (show_x_title) "Nucleotide position" else NULL,
      y = "% Substitutions"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5.5),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 6),
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.25, "lines")
    )
}

plot_A <- make_family_plot("Adenine editors")
plot_C <- make_family_plot("Cytosine editors")

n_A <- length(intersect(adenine_editors,  unique(summary_data$base_editor)))
n_C <- length(intersect(cytosine_editors, unique(summary_data$base_editor)))

combined_plot <- (plot_A / plot_C) +
  theme(legend.position = "right")

combined_plot <- combined_plot +
  plot_annotation(
    title = "Total Substitutions Across Base Editors"
  ) 

ggsave(
  file.path(plots_dir, "substitutions_total_compressed.pdf"),
  plot = combined_plot, width = 10, height = 12, units = "in", limitsize = FALSE
)

rm(main_data, adenine_editors, cytosine_editors,
   build_family_data, plot_data_raw, summary_data, make_family_plot,
   plot_A, plot_C, combined_plot, n_A, n_C)

#### Heatmap ----

main_data_d5 <- dat_long %>%
  filter(Modification == "Substitutions", dosage == "D5")

adenine_editors  <- c(A_editors, AC_editors)
cytosine_editors <- c(C_editors, AC_editors)

label_sf <- function(x) {
  pct <- x * 100
  formatC(pct, digits = 2, format = "fg", flag = "#")
}

build_panel_data <- function(editors, letter) {
  
  raw <- main_data_d5 %>%
    filter(base_editor %in% editors, base_letter == letter)
  
  summary_data <- summarise_stats(raw, c("base_editor", "target", "position", "pos_num"))
  
  # Keep target/position only if >=1 base editor hits >=1% mean editing at D5
  qualifying_tp <- summary_data %>%
    group_by(target, position) %>%
    summarise(max_mean = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
    filter(max_mean >= 0.01) %>%
    select(target, position)
  
  summary_data %>%
    semi_join(qualifying_tp, by = c("target", "position")) %>%
    mutate(
      position    = fct_drop(position),
      base_editor = factor(base_editor, levels = rev(editors)),   # rev() so first editor plots on top
      label       = label_sf(mean_val)   # every cell labeled, 3 sig figs, no >=1% restriction
    )
}

adenine_data  <- build_panel_data(adenine_editors,  "A")
cytosine_data <- build_panel_data(cytosine_editors, "C")

make_heatmap_panel <- function(panel_data, title, hue, show_x_title = TRUE) {
  
  ggplot(panel_data, aes(x = position, y = base_editor, fill = mean_val)) +
    geom_tile(color = "grey40", linewidth = 0.2) +
    geom_text(aes(label = label), size = 5, color = "black") +
    facet_grid(~ target, scales = "free_x", space = "free_x") +
    scale_fill_gradient(
      low = "white", high = hue, limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      name = "% Substitutions\n(D5)"
    ) +
    labs(
      title = title,
      x = if (show_x_title) "Nucleotide position" else NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title        = element_text(size = 12, face = "bold"),
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y       = element_text(size = 14),
      strip.background  = element_rect(fill = "grey90"),
      strip.text        = element_text(size = 8),
      panel.grid         = element_blank(),
      panel.spacing      = unit(0.2, "lines"),
      legend.key.height  = unit(0.7, "cm")
    )
}

adenine_plot  <- make_heatmap_panel(adenine_data,  "Adenine editing",  hue = "#7f007f", show_x_title = FALSE)
cytosine_plot <- make_heatmap_panel(cytosine_data, "Cytosine editing", hue = "#00771d", show_x_title = TRUE)

d5_ac_heatmap <- (adenine_plot / cytosine_plot) +
  plot_annotation(
    title = "D5 (10 \u00b5L) Editing \u2014 Adenine vs. Cytosine Base Editors"
  )

ggsave(
  file.path(plots_dir, "D5_adenine_cytosine_heatmap.pdf"),
  plot = d5_ac_heatmap, width = 12, height = 8, units = "in", limitsize = FALSE
)

rm(main_data_d5, adenine_editors, cytosine_editors, label_sf, build_panel_data,
   adenine_data, cytosine_data, make_heatmap_panel,
   adenine_plot, cytosine_plot, d5_ac_heatmap)

### 2b. Substitution composition proportions ----

#### Substitution composition plot ----

composition_editors <- list(
  A = c("ABE8e", "ABEmax-NG", "ABE8e-SpRY", "AYBE", "AXBE"),
  C = c("CBE6d", "TadCBEd", "CGBE1","CGBE1-NG", "miniCGBE1-SpRY")
)

base_hues <- c(A = "darkgreen", T = "firebrick", C = "dodgerblue4", G = "#d8a822")

get_d5_qualifying_positions <- function(editors, letter) {
  dat_long %>%
    filter(Modification == "Substitutions", dosage == "D5",
           base_editor %in% editors, base_letter == letter) %>%
    summarise_stats(c("base_editor", "target", "position")) %>%
    group_by(target, position) %>%
    summarise(max_mean = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
    filter(max_mean >= 0.01) %>%
    select(target, position)
}

build_composition_data <- function(editors, letter) {
  
  qualifying_positions <- get_d5_qualifying_positions(editors, letter)
  
  main_data <- dat_long %>%
    filter(Modification %in% c("A", "T", "C", "G")) %>%
    filter(Modification != base_letter) %>%
    filter(dosage == "D5") %>%
    filter(base_editor %in% editors) %>%
    semi_join(qualifying_positions, by = c("target", "position"))
  
  summary_data <- main_data %>%
    group_by(base_editor, target, position, pos_num, Modification) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
    rename(subst_base = Modification) %>%
    mutate(subst_base = factor(subst_base, levels = c("A", "T", "C", "G")))
  
  # Renormalize: express each substitution base as % of TOTAL substitutions
  summary_data <- summary_data %>%
    group_by(base_editor, target, position, pos_num) %>%
    mutate(
      total_sub_val = sum(mean_val, na.rm = TRUE),
      mean_val      = if_else(total_sub_val > 0, mean_val / total_sub_val, 0)
    ) %>%
    ungroup() %>%
    mutate(
      base_editor = factor(base_editor, levels = editors),
      position    = droplevels(position)
    )
  
  summary_data
}

adenine_comp  <- build_composition_data(composition_editors$A, "A")
cytosine_comp <- build_composition_data(composition_editors$C, "C")

build_total_label_data <- function(comp_data) {
  comp_data %>%
    distinct(base_editor, target, position, pos_num, total_sub_val) %>%
    mutate(label = paste0(signif(total_sub_val * 100, 2), "%"))
}

adenine_labels  <- build_total_label_data(adenine_comp)
cytosine_labels <- build_total_label_data(cytosine_comp)

make_composition_panel <- function(comp_data, label_data, title) {
  
  ggplot(comp_data, aes(x = position, y = mean_val, fill = subst_base)) +
    geom_col(color = "black", linewidth = 0.15, position = "stack") +
    geom_text(
      data = label_data,
      aes(x = position, y = 1.08, label = label),
      inherit.aes = FALSE, size = 2.8, vjust = 0
    ) +
    facet_grid(
      rows = vars(base_editor), cols = vars(target),
      scales = "free_x", space = "free_x"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      expand = expansion(mult = c(0, 0.25)),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0)
    ) +
    scale_fill_manual(values = base_hues, name = "Substituted to", drop = FALSE) +
    labs(title = title, x = "Nucleotide position", y = "% Substitutions") +
    theme_bw(base_size = 9) +
    theme(
      plot.title        = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y       = element_text(size = 6),
      strip.background  = element_rect(fill = "grey90"),
      strip.text        = element_text(size = 6),
      strip.placement   = "outside",
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(0.15, "lines")
    )
}

adenine_panel  <- make_composition_panel(adenine_comp,  adenine_labels,  "Adenine editors")
cytosine_panel <- make_composition_panel(cytosine_comp, cytosine_labels, "Cytosine editors")

composition_plot <- (adenine_panel / cytosine_panel)

ggsave(
  file.path(plots_dir, "abe_cbe_composition.pdf"),
  plot = composition_plot, width = 8, height = 7, units = "in", limitsize = FALSE
)

rm(composition_editors, base_hues, get_d5_qualifying_positions,
   build_composition_data, adenine_comp, cytosine_comp,
   build_total_label_data, adenine_labels, cytosine_labels,
   make_composition_panel, adenine_panel, cytosine_panel, composition_plot)


### 2c. Substitution frequency across spacer ----

main_data <- dat_long %>% filter(Modification == "Substitutions", dosage == "D5")

# Order rows by editor family (A-editors -> C-editors -> SPACE) so the panel
# reads as categorized by BE type, top to bottom
editor_order <- c(A_editors, C_editors, AC_editors)
main_data <- main_data %>%
  mutate(base_editor = factor(base_editor, levels = editor_order))

# Mean + SEM across replicates, per base_editor/target/position (D5 only)
summary_data <- summarise_stats(main_data, c("base_editor", "target", "position", "pos_num"))

compressed_lineplot <- ggplot(
  summary_data,
  aes(x = pos_num, y = mean_val, color = target, group = target)
) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = pmax(mean_val - sem_val, 0), ymax = mean_val + sem_val),
    width = 0.25, linewidth = 0.3
  ) +
  geom_point(
    data = main_data,
    aes(x = pos_num, y = value, color = target, group = target),
    size = 0.7, alpha = 0.6,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  facet_grid(
    rows   = vars(base_editor),
    scales = "free",      # free x (positions differ per target) AND free y
    space  = "free_x"     # panel widths scale with number of positions
  ) +
  scale_color_manual(values = target_colors, name = "Target") +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "D5 Substitution Frequency Across Spacer - All Base Editors (all targets, one panel)",
    x = "Nucleotide position",
    y = "% Substitutions"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y       = element_text(size = 5.5),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(size = 6),
    strip.placement   = "outside",
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(0.15, "lines")
  )

ggsave(
  file.path(plots_dir, "substitutions_spacers_lineplot_compressed_D5.pdf"),
  plot = compressed_lineplot, width = 15, height = 24, units = "in", limitsize = FALSE
)

rm(main_data, summary_data, editor_order, compressed_lineplot)

### Substitution percentage across dosage (A- and C-editors only) ----

sub_data <- editing_summary %>%
  rename(base_editor = Base_Editor, dosage = Dose, target = Name, replicate = Replicate) %>%
  filter(!base_editor %in% AC_editors) %>%
  mutate(
    Substitutions_pct = Substitutions / Reads_aligned * 100,
    dosage = factor(dosage, levels = dosage_levels),
    target = factor(target, levels = c("HEK3", "B2M", "PDCD1"))
  )

editor_order <- c(A_editors, C_editors)
editor_family <- bind_rows(
  tibble(base_editor = A_editors, family = "A-editors"),
  tibble(base_editor = C_editors, family = "C-editors")
)

sub_data <- sub_data %>%
  left_join(editor_family, by = "base_editor") %>%
  mutate(
    base_editor = factor(base_editor, levels = editor_order),
    family      = factor(family, levels = c("A-editors", "C-editors"))
  )

sub_summary <- sub_data %>%
  group_by(base_editor, family, target, dosage) %>%
  summarise(
    mean_sub_pct = mean(Substitutions_pct, na.rm = TRUE),
    sem_sub_pct  = sd(Substitutions_pct, na.rm = TRUE) / sqrt(n()),
    n_rep        = n(),
    .groups = "drop"
  )

substitutions_pct_plot <- ggplot(
  sub_summary,
  aes(x = dosage, y = mean_sub_pct, color = base_editor, group = base_editor)
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = pmax(mean_sub_pct - sem_sub_pct, 0), ymax = mean_sub_pct + sem_sub_pct),
    width = 0.15
  ) +
  geom_point(
    data = sub_data,
    aes(x = dosage, y = Substitutions_pct, color = base_editor),
    size = 0.8, alpha = 0.35, inherit.aes = FALSE
  ) +
  facet_grid(rows = vars(family), cols = vars(target), scales = "free_y") +
  scale_x_discrete(labels = dosage_labels_x_axis) +   # actual VLP volumes (µL) instead of D0-D5 codes
  scale_color_manual(
    values = editor_palette, name = "Base editor",
    guide = guide_legend(nrow = 2, byrow = TRUE)   # 2 rows x 5 columns
  ) +
  labs(
    title = "Substitution Percentage Across Dosage",
    x = "VLP dosage (\u00b5L)", y = "% Substitutions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(size = 8),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(0.3, "lines"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8)
  )

ggsave(
  file.path(plots_dir, "substitutions_percentage_dosage.pdf"),
  plot = substitutions_pct_plot, width = 10, height = 8, units = "in", limitsize = FALSE
)

rm(sub_data, editor_order, editor_family, sub_summary, substitutions_pct_plot)


### 2d. Substitution composition table ----
#
# JUDGMENT CALL: this block needs a `qualifying_positions` object with the
# individual A/T/C/G substitution-outcome rows intact (Modification, value,
# dosage, replicate, etc.) so it can compute per-position substitution
# composition. That's different from the `qualifying_positions` built in
# Section 3 (which collapses to a bare distinct target/position list with no
# Modification/value columns, and is unused there beyond being rm()'d) -
# reusing that version would fail here for lack of columns to sum/average.
# Built to mirror the "qualifying position" logic used for the composition
# plot in Section 2b: restrict to the A/T/C/G outcome rows, then keep only
# base_editor x base_letter combinations appropriate to each editor's panel
# (A-editors -> A positions, C-editors -> C positions, AC-editors -> either).

qualifying_positions <- dat_long %>%
  filter(Modification %in% c("A", "T", "C", "G")) %>%
  get_qualifying_positions()

substitution_composition <- qualifying_positions %>%
  filter(Modification != base_letter) %>%
  group_by(base_editor, dosage, replicate, target, position, pos_num) %>%
  mutate(total_subs = sum(value, na.rm = TRUE)) %>%
  mutate(substitution_percent = value / total_subs) %>%
  ungroup() %>%
  select(
    base_editor,
    dosage,
    replicate,
    target,
    position,
    pos_num,
    subst_base = Modification,
    substitution_percent,
    total_subs               # <-- kept: needed to label total edit % on the plot
  )

# Add average across replicates:
substitution_composition <- substitution_composition %>%
  group_by(base_editor, dosage, target, position, subst_base) %>%
  mutate(avg_substitution_percent = mean(substitution_percent, na.rm = TRUE)) %>%
  ungroup()

# Simplified data w/o replicate level
substitution_composition_avg <- substitution_composition %>%
  select(base_editor, dosage, target, position, subst_base, pos_num, avg_substitution_percent) %>%
  distinct()

# Average total edit % (across replicates) - one value per base_editor/dosage/target/position,
# independent of subst_base
total_edit_avg <- substitution_composition %>%
  distinct(base_editor, dosage, replicate, target, position, pos_num, total_subs) %>%
  group_by(base_editor, dosage, target, position, pos_num) %>%
  summarise(avg_total_subs = mean(total_subs, na.rm = TRUE), .groups = "drop")

# CBE edits
sub_composition_summary <- substitution_composition_avg %>%
  filter(base_editor %in% c("CBE6d", "TadCBEd")) %>%
  group_by(base_editor, dosage, target, position, pos_num) %>%
  summarise(
    C_to_T_percent = sum(avg_substitution_percent[subst_base == "T"], na.rm = TRUE),
    C_to_nonT_percent = sum(avg_substitution_percent[subst_base != "T"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(total_edit_avg, by = c("base_editor", "dosage", "target", "position", "pos_num"))

rm(qualifying_positions, substitution_composition, substitution_composition_avg, total_edit_avg)
# NOTE: sub_composition_summary is intentionally NOT removed here (unlike
# the rest of this section) - it's the actual canonical/non-canonical CBE
# editing summary and isn't saved to a file anywhere in this script, so
# clearing it would discard the only output this block produces.


## 3. Deletion plots ===================================================

qualifying_positions <- dat_long %>%
  distinct(base_editor, target, position, pos_num, base_letter, target_position) %>%
  get_qualifying_positions()

main_data <- dat_long %>% filter(Modification == "Deletions")

mean_data <- main_data %>%
  group_by(base_editor, target, position, pos_num, target_position, dosage) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_type = "Mean", grp = paste(dosage, "mean", sep = "_"))

replicate_data <- main_data %>%
  mutate(group_type = "Replicate", grp = paste(dosage, replicate, sep = "_")) %>%
  select(base_editor, target, position, pos_num, target_position, dosage, value, group_type, grp)

line_data <- bind_rows(replicate_data, mean_data)

y_axis_min <- 0.005

main_data <- dat_long %>%
  filter(Modification == "Deletions", dosage == "D5", !base_editor %in% AC_editors)

# Order rows by editor family (A-editors -> C-editors) so the panel
# reads as categorized by BE type, top to bottom
editor_order <- c(A_editors, C_editors)
main_data <- main_data %>%
  mutate(base_editor = factor(base_editor, levels = editor_order))

# Mean + SEM across replicates, per base_editor/target/position (D5 only)
summary_data <- summarise_stats(main_data, c("base_editor", "target", "position", "pos_num"))

compressed_deletion_lineplot <- ggplot(
  summary_data,
  aes(x = pos_num, y = mean_val, color = target, group = target)
) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = pmax(mean_val - sem_val, 0), ymax = mean_val + sem_val),
    width = 0.25, linewidth = 0.3
  ) +
  geom_point(
    data = main_data,
    aes(x = pos_num, y = value, color = target, group = target),
    size = 0.7, alpha = 0.6,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  facet_grid(
    rows   = vars(base_editor),
    scales = "free", space = "free_x"
  ) +
  scale_color_manual(values = target_colors, name = "Target") +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(
    labels = function(x) paste0(formatC(signif(x * 100, 2), format = "fg"), "%"),
    limits = c(0, NA)
  ) +
  labs(
    title = "D5 Deletion Frequency Across Spacer - All Base Editors (all targets, one panel)",
    x = "Nucleotide position", y = "% Deletions"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y       = element_text(size = 5.5),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(size = 6),
    strip.placement   = "outside",
    panel.spacing     = unit(0.15, "lines"),
    panel.grid.minor  = element_blank()
  )

ggsave(
  file.path(plots_dir, "deletions_spacers_lineplot_compressed_D5.pdf"),
  plot = compressed_deletion_lineplot, width = 15, height = 20, units = "in", limitsize = FALSE
)

rm(main_data, mean_data, replicate_data, line_data, editor_order, y_axis_min,
   qualifying_positions, summary_data, compressed_deletion_lineplot)


## 4. Cross-metric comparisons (editing_summary-derived) ===============

### Substitutions vs. % Deletions (scatter, A- and C-editors only) ----

scatter_data <- editing_summary %>%
  dplyr::filter(Dose != "D0") %>%
  rename(base_editor = Base_Editor, dosage = Dose, target = Name, replicate = Replicate) %>%
  filter(!base_editor %in% AC_editors) %>%
  mutate(
    Substitutions_pct = Substitutions / Reads_aligned * 100,
    Deletions_pct        = Deletions / Reads_aligned * 100,
    Indels_pct        = (Insertions + Deletions) / Reads_aligned * 100,
    dosage = factor(dosage, levels = dosage_levels),
    target = factor(target, levels = c("HEK3", "B2M", "PDCD1"))
  )

editor_order <- c(A_editors, C_editors)
editor_family <- bind_rows(
  tibble(base_editor = A_editors, family = "A-editors"),
  tibble(base_editor = C_editors, family = "C-editors")
)

scatter_data <- scatter_data %>%
  left_join(editor_family, by = "base_editor") %>%
  mutate(
    base_editor = factor(base_editor, levels = editor_order),
    family      = factor(family, levels = c("A-editors", "C-editors"))
  )

scatter_palette <- editor_palette[names(editor_palette) %in% c(A_editors, C_editors)]

subs_vs_indels_plot <- ggplot(scatter_data, aes(x = Substitutions_pct, y = Deletions_pct, color = base_editor)) +
  geom_point(size = 2.5, alpha = 0.6) +
  ggh4x::facet_grid2(rows = vars(family), cols = vars(target), scales = "free") +
  scale_color_manual(
    values = scatter_palette, name = "Base editor",
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    title = "% Substitutions vs. % Deletions",
    x = "% Substitutions", y = "% Deletions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(size = 8),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(0.3, "lines"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8)
  )

ggsave(
  file.path(plots_dir, "substitutions_vs_dels_scatter.pdf"),
  plot = subs_vs_indels_plot, width = 12, height = 8, units = "in", limitsize = FALSE
)

rm(scatter_data, editor_order, editor_family, scatter_palette, subs_vs_indels_plot)

### Deletion : substitution ratio ----

ratio_data <- editing_summary %>%
  rename(base_editor = Base_Editor, dosage = Dose, target = Name, replicate = Replicate) %>%
  filter(
    !base_editor %in% AC_editors,              # excluded: SPACE dual-base editors underperform
    dosage == "D5"                             # 10 uL only
  ) %>%
  mutate(
    Substitutions_pct = Substitutions / Reads_aligned * 100,
    Insertions_pct    = Insertions    / Reads_aligned * 100,
    Deletions_pct     = Deletions     / Reads_aligned * 100,
    Indels_pct        = Insertions_pct + Deletions_pct,
    ratio  = if_else(Substitutions_pct > 0, Deletions_pct / Substitutions_pct, NA_real_),
    target = factor(target, levels = c("HEK3", "B2M", "PDCD1"))
  ) %>%
  filter(!is.na(ratio))

n_dropped <- nrow(editing_summary %>% filter(Dose == "D5", !Base_Editor %in% AC_editors)) - nrow(ratio_data)
if (n_dropped > 0) {
  message(n_dropped, " D5 sample(s) dropped: zero substitutions (ratio undefined).")
}

editor_order <- c(A_editors, C_editors)
editor_family <- bind_rows(
  tibble(base_editor = A_editors,  family = "Adenine editors"),
  tibble(base_editor = C_editors,  family = "Cytosine editors")
)

ratio_data <- ratio_data %>%
  left_join(editor_family, by = "base_editor") %>%
  mutate(
    base_editor = factor(base_editor, levels = editor_order),
    family      = factor(family, levels = c("Adenine editors", "Cytosine editors"))
  )

ratio_summary <- ratio_data %>%
  group_by(base_editor, family, target) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    sem_ratio  = sd(ratio, na.rm = TRUE) / sqrt(n()),
    n_rep      = n(),
    .groups = "drop"
  )

ratio_palette <- editor_palette[names(editor_palette) %in% c(A_editors, C_editors)]

indel_subs_ratio_plot <- ggplot(
  ratio_summary,
  aes(x = base_editor, y = mean_ratio, fill = base_editor)
) +
  geom_col(width = 0.7, linewidth = 0.3, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = pmax(mean_ratio - sem_ratio, 0), ymax = mean_ratio + sem_ratio),
    width = 0.25, linewidth = 0.4
  ) +
  geom_point(
    data = ratio_data,
    aes(x = base_editor, y = ratio),
    size = 0.8, alpha = 0.35, position = position_jitter(width = 0.15), inherit.aes = FALSE
  ) +
  ggh4x::facet_nested(
    cols = vars(target, family),
    scales = "free", space = "free_x"
  ) +
  scale_fill_manual(values = ratio_palette, name = "Base editor") +
  labs(
    title = "Deletion-to-Substitution Ratio at 10 \u00b5L (D5)",
    x = "Base editor", y = "Deletions : Substitutions ratio"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(size = 8),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(0.3, "lines"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position    = "none"
  )

ggsave(
  file.path(plots_dir, "del_substitution_ratio_D5.pdf"),
  plot = indel_subs_ratio_plot, width = 10, height = 6, units = "in", limitsize = FALSE
)

rm(editing_summary, ratio_data, n_dropped, editor_order, editor_family,
   ratio_summary, ratio_palette, indel_subs_ratio_plot)

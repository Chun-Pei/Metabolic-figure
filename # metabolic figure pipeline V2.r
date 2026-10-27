############################################################
# Nucleotide isotope plotting pipeline — Prism-style
############################################################

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

############################################################
# USER INPUT — edit this section only
############################################################

INPUT_FILE <- file.choose()
OUT_DIR <- dirname(INPUT_FILE)
OUT_STEM <- file.path(OUT_DIR, tools::file_path_sans_ext(basename(INPUT_FILE)))

conditions <- tribble(
    ~code     , ~label                , ~order ,
    "UL"      , "Unlabelled"          ,      1 ,
    "N"       , "Normoxia"            ,      2 ,
    "N MTX"   , "Normoxia + MTX"      ,      3 ,
    "H4"      , "Hypoxia 4h"          ,      4 ,
    "H4 MTX"  , "Hypoxia 4h + MTX"    ,      5 ,
    "H6"      , "Hypoxia 6h"          ,      6 ,
    "H6 MTX"  , "Hypoxia 6h + MTX"    ,      7 ,
    "H48"     , "Hypoxia 48h"         ,      9 ,
    "H48 MTX" , "Hypoxia 48h + MTX"   ,     11 ,
    "0.4 mM"  , "0.4mM serine"        ,     12 ,
    "0.8 mM"  , "0.8mM serine"        ,     13 ,
    "UL KO"   , "Unlabelled HIF1B KO" ,     14 ,
    "N KO"    , "Normoxia KO"         ,     15 ,
    "H48 KO"  , "Hypoxia 48h KO"      ,     16
)

METABOLITES <- c(
    "AMP",
    "ADP",
    "ATP",
    "Serine",
    "dTMP",
    "SAM",
    "SAH",
    "Methionine"
)

ISOTOPES <- c("M+0", "M+1", "M+2", "M+3")

iso_colors <- c(
    "M+0" = "#FFFFFF",
    "M+1" = "#F4ABAB",
    "M+2" = "#E03131",
    "M+3" = "#7A0000"
)

# CONDITIONS TO SHOW — edit this line to change which appear in the figure
KEEP_CONDITIONS <- c(
    "Unlabelled",
    "Normoxia",
    "Hypoxia 6h",
    "Hypoxia 48h",
    "Hypoxia 6h + MTX",
    "Hypoxia 48h + MTX",
    "Unlabelled HIF1B KO",
    "Normoxia KO",
    "Hypoxia 48h KO"
)

############################################################
# Auto-built — no changes needed below
############################################################

CONDITION_MAP <- setNames(conditions$label, conditions$code)
condition_order <- conditions$label[order(conditions$order)]

############################################################
# PIPELINE
############################################################

# 1  Load data
data <- read_excel(INPUT_FILE)
colnames(data)[1] <- "Sample"

# 2  Coerce isotope columns to numeric; treat "NF" (not found) as NA
iso_cols <- grep("^%", colnames(data), value = TRUE)
data[iso_cols] <- lapply(data[iso_cols], function(x) {
    x <- as.character(x)
    x[x == "NF"] <- NA
    as.numeric(x)
})

# 3  Extract condition label
data$Condition <- sub("^[0-9]+\\. ", "", data$Sample)
data$Condition <- sub(" [0-9]+$", "", data$Condition)

# 4  Tidy isotope data
iso_data <- data %>%
    pivot_longer(
        cols = all_of(iso_cols),
        names_to = "variable",
        values_to = "fraction"
    ) %>%
    mutate(
        metabolite = str_extract(variable, paste(METABOLITES, collapse = "|")),
        isotope = str_extract(variable, "M\\+[0-9]+")
    ) %>%
    filter(!is.na(metabolite), !is.na(isotope)) %>%
    select(Sample, Condition, metabolite, isotope, fraction) %>%
    # Renormalise: within each sample x metabolite, rescale detected isotopes
    # so they sum to 100. NF (NA) isotopes remain NA and contribute nothing
    # to the stack, but the detected ones always reach 100% together.
    group_by(Sample, metabolite) %>%
    mutate(fraction = fraction / sum(fraction, na.rm = TRUE) * 100) %>%
    ungroup()

# 5  Sort
iso_sorted <- iso_data %>% arrange(metabolite, Condition, isotope)

# 6  Summarise (mean ± SEM per condition × metabolite × isotope)
iso_summary <- iso_sorted %>%
    group_by(Condition, metabolite, isotope) %>%
    summarise(
        mean_fraction = mean(fraction, na.rm = TRUE),
        sem = sd(fraction, na.rm = TRUE) / sqrt(sum(!is.na(fraction))),
        .groups = "drop"
    )

iso_summary$isotope <- factor(iso_summary$isotope, levels = ISOTOPES)
iso_summary$metabolite <- factor(iso_summary$metabolite, levels = METABOLITES)

# 7  Cumulative sums
iso_summary <- iso_summary %>%
    arrange(metabolite, Condition, isotope) %>%
    group_by(metabolite, Condition) %>%
    mutate(
        cumsum_top = cumsum(replace_na(mean_fraction, 0)),
        bar_mid = cumsum_top - replace_na(mean_fraction, 0) / 2
    ) %>%
    ungroup()

# 8  Recode & factor conditions
iso_summary$Condition <- recode(iso_summary$Condition, !!!CONDITION_MAP)
iso_summary <- iso_summary %>%
    mutate(Condition = factor(Condition, levels = condition_order))

# 8b  Filter to selected conditions only
iso_summary <- iso_summary %>%
    filter(Condition %in% KEEP_CONDITIONS) %>%
    mutate(Condition = droplevels(Condition))

############################################################
# 9  Stacked bar plot — error bar at every isotope's cumsum_top (Script 1 style)
############################################################

p_iso <- ggplot(
    iso_summary,
    aes(x = Condition, y = mean_fraction, fill = isotope)
) +

    geom_bar(
        stat = "identity",
        color = "black",
        width = 0.7,
        position = position_stack(reverse = TRUE),
        na.rm = TRUE
    ) +

    geom_errorbar(
        aes(y = cumsum_top, ymin = cumsum_top - sem, ymax = cumsum_top + sem),
        width = 0.18,
        linewidth = 0.55,
        colour = "black"
    ) +

    facet_wrap(~metabolite, nrow = 1) +

    scale_fill_manual(values = iso_colors) +

    labs(y = "Fraction of total", x = "") +

    theme_classic(base_size = 14) +
    theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))

############################################################
# 10  Total-abundance panel
############################################################

total_data <- data %>%
    select(Sample, Condition, starts_with("Total")) %>%
    pivot_longer(
        cols = starts_with("Total"),
        names_to = "metabolite",
        values_to = "abundance"
    ) %>%
    mutate(metabolite = str_remove(metabolite, "Total ")) %>%
    filter(
        !is.na(abundance),
        metabolite %in% METABOLITES
    )

total_data$metabolite <- factor(total_data$metabolite, levels = METABOLITES)
total_data$Condition <- recode(total_data$Condition, !!!CONDITION_MAP)
total_data$Condition <- factor(total_data$Condition, levels = condition_order)

total_data <- total_data %>%
    filter(Condition %in% KEEP_CONDITIONS) %>%
    mutate(Condition = droplevels(Condition))

p_total <- ggplot(total_data, aes(x = Condition, y = abundance)) +
    geom_boxplot(
        width = 0.55,
        fill = "grey88",
        outlier.shape = NA,
        linewidth = 0.45
    ) +
    geom_jitter(width = 0.12, size = 1.6, colour = "black") +
    facet_wrap(~metabolite, nrow = 1, scales = "free_y") +
    labs(y = "Total abundance", x = NULL) +
    theme_classic(base_size = 13) +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.line = element_line(colour = "black")
    )

############################################################
# 11  Combine and save
############################################################

final_plot <- p_iso / p_total + plot_layout(heights = c(1.2, 1))

ggsave(paste0(OUT_STEM, "_figure.pdf"), final_plot, width = 10, height = 8)
ggsave(
    paste0(OUT_STEM, "_figure.png"),
    final_plot,
    width = 10,
    height = 8,
    dpi = 600
)

message("Done — files saved:")
message("  ", OUT_STEM, "_figure.pdf")
message("  ", OUT_STEM, "_figure.png")

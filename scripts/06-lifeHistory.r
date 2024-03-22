source("scripts/00-setup.R")

metadata <- readRDS("./processed/metadata.rds")
ps_norare <- readRDS("./processed/ps_pruned_samplefiltered.rds")
ps <- readRDS("./processed/ps_pruned_rared10k.rds")

length(unique(sample_data(ps)$Isoline))
# 22

library(GGally)

# Function to plot lower panel of ggpairs with custom colors
lower_panelfn <- function(data, mapping, method = "lm", ...) {
    p <- ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_point(colour = "gray25") +
        ggplot2::geom_smooth(
            method = method, color = "darkred",
            fill = "#aa6868", ...
        )
    p
}

# Correlations between life history traits
lh_correlations <- data.frame(sample_data(ps)) %>%
    select(Isoline, AvLF, EarlyRS, CSearly, CSSlope, Rsen) %>%
    na.omit() %>%
    unique() %>%
    select(-Isoline) %>%
    rename(
        `Average Lifespan` = AvLF,
        `Early Rep Success` = EarlyRS,
        `Early Climb Speed` = CSearly,
        `Functional Aging` = CSSlope,
        `Rep Senescence` = Rsen
    ) %>%
    ggpairs(
        lower = list(continuous = wrap(lower_panelfn)),
        upper = list(continuous = wrap("cor", method = "spearman")),
        diag = list(continuous = wrap("densityDiag", fill = "gray75"))
    ) +
    theme_bw()

# Save plot
plotsave(lh_correlations, "./results/life_history/correlations",
    width = 16,
    height = 16
)

# Check the summary statistics for life history traits
data.frame(sample_data(ps)) %>%
    select(Isoline, AvLF, EarlyRS, CSearly, CSSlope, Rsen) %>%
    na.omit() %>%
    unique() %>%
    select(-Isoline) %>%
    summary()

#      AvLF          EarlyRS          CSearly          CSSlope             Rsen
# Min.   :25.14   Min.   :0.1275   Min.   :0.2233   Min.   :-0.6467   Min.   :-0.10001
# 1st Qu.:32.41   1st Qu.:0.3971   1st Qu.:0.3444   1st Qu.:-0.5360   1st Qu.: 0.09491
# Median :35.96   Median :0.4738   Median :0.5067   Median :-0.4286   Median : 0.19027
# Mean   :37.77   Mean   :0.4932   Mean   :0.5243   Mean   :-0.4102   Mean   : 0.16416
# 3rd Qu.:41.39   3rd Qu.:0.6120   3rd Qu.:0.7075   3rd Qu.:-0.2888   3rd Qu.: 0.24310
# Max.   :58.61   Max.   :0.7852   Max.   :0.8167   Max.   :-0.1238   Max.   : 0.38166

# Calculate the summary statistics for life history traits
summary_stats <- data.frame(sample_data(ps)) %>%
    select(Isoline, AvLF, EarlyRS, CSearly, CSSlope, Rsen) %>%
    na.omit() %>%
    unique() %>%
    select(-Isoline) %>%
    gather(key = "Life History", value = value) %>%
    group_by(`Life History`) %>%
    summarise(
        Minimum = min(value),
        Maximum = max(value),
        Median = median(value),
        Q1 = quantile(value, 0.25),
        Q3 = quantile(value, 0.75),
        IQR = IQR(value),
        Mean = mean(value),
        StDev = sd(value)
    ) %>%
    mutate(`Max/Min` = Maximum / Minimum) %>%
    gather(key = "stat", value = "value", -`Life History`) %>%
    spread(key = `Life History`, value) %>%
    mutate(stat = factor(stat, levels = c(
        "Minimum", "Maximum", "Max/Min", "Median",
        "Q1", "Q3", "IQR", "Mean", "StDev"
    ))) %>%
    arrange(stat) %>%
    rename(
        ` ` = stat,
        `Average Lifespan` = AvLF,
        `Early Rep Success` = EarlyRS,
        `Early Climb Speed` = CSearly,
        `Functional Aging` = CSSlope,
        `Rep Senescence` = Rsen
    )

# Save summary statistics
tablesave(summary_stats, "./results/life_history/summary_stats")

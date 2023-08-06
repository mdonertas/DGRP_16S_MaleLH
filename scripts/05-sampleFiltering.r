source("scripts/00-setup.R")
library(phyloseq)

ps_before <- readRDS("./processed/ps.rds")
ps <- readRDS("./processed/ps_length_phyla_treepruned_prev.rds")

# Plot the library size across samples
libsize_across_samples <- data.frame(sums = sample_sums(ps)) %>%
    ggplot(aes(x = sums)) +
    geom_histogram(color = "gray80", bins = 50, linewidth = 0.25) +
    xlab("Library size (log10 scale)") +
    ylab("Frequency") +
    geom_vline(
        xintercept = 10000, linetype = "dotted", color = "darkred",
        linewidth = 0.5
    ) +
    scale_x_continuous(labels = scales::comma, trans = "log10") +
    theme_pubr(base_size = 12)

plotsave(libsize_across_samples, "./results/EDA/libsize_across_samples",
    width = 8, height = 8
)

table(sample_sums(ps) <= 10000)

# Filter samples with less than 10,000 reads
ps_sample10kpruned <- prune_samples(names(which(sample_sums(ps) > 10000)), ps)

objsave(ps_sample10kpruned, "./processed/ps_pruned_sample10kfiltered")

# only keep samples where we can analyse both early and late with 3 samples
samples2keep <- (data.frame(sample_data(ps_sample10kpruned)) %>%
    group_by(Isoline, age) %>%
    summarise(n = length(unique(sampleID))) %>%
    spread(age, n) %>%
    filter(early >= 3 & late >= 3) %>%
    left_join(sample_data(ps_sample10kpruned)))$sampleID
ps_samplepruned <- subset_samples(
    ps_sample10kpruned,
    sample_names(ps_sample10kpruned) %in% samples2keep
)

objsave(ps_samplepruned, "./processed/ps_pruned_samplefiltered")

# Create 100 rarefied tables with 10,000 reads per sample
ps10k <- lapply(1:100, function(i) {
    rarefy_even_depth(ps_samplepruned, 10000)
})

# Save one of the rarefied tables for use in the manuscript
objsave(ps10k[[1]], "./processed/ps_pruned_rared10k")

# Save all rarefied tables
objsave(ps10k, "./processed/ps_pruned_rared10k_100tables")

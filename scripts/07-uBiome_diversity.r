source("scripts/00-setup.R")

library(phyloseq)

ps <- readRDS("./processed/ps_pruned_rared10k.rds")
ps_all <- readRDS("./processed/ps_pruned_rared10k_100tables.rds")

# which phyla are present and most dominant?
all(sapply(ps_all, function(psx) {
    sort(unique(tax_table(psx)[, "Phylum"]))
}) == c("Actinobacteriota", "Cyanobacteria", "Firmicutes", "Proteobacteria"))
# in this filtered set, we detect only 4 phyla, present in all rarefied tables

ps_df <- reshape2::melt(otu_table(ps)) %>%
    set_names(c("sampleID", "ASV", "count")) %>%
    left_join(sample_data(ps)) %>%
    left_join(mutate(data.frame(tax_table(ps)), ASV = taxa_names(ps)))

numtaxa <- ps_df %>%
    filter(count > 0) %>%
    group_by(sampleID, age) %>%
    summarise(
        Kingdom = length(unique(setdiff(Kingdom, NA))),
        Phylum = length(unique(setdiff(Phylum, NA))),
        Class = length(unique(setdiff(Class, NA))),
        Order = length(unique(setdiff(Order, NA))),
        Family = length(unique(setdiff(Family, NA))),
        Genus = length(unique(setdiff(Genus, NA))),
        `Missing Sp` = sum(is.na(Species)),
        Species = length(unique(setdiff(Species, NA))),
        ASV = length(unique(setdiff(ASV, NA)))
    ) %>%
    select(-Kingdom, -Class, -Order, -Family) %>%
    gather(key = "taxon", value = "count", -sampleID, -age) %>%
    mutate(taxon = factor(taxon, levels = c(
        "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV",
        "Missing Sp"
    )))

numtaxa %>%
    group_by(taxon) %>%
    summarise(
        min = min(count), max = max(count), mean = mean(count),
        median = median(count)
    )

#  taxon        min   max   mean median
#  <fct>      <int> <int>  <dbl>  <int>
# 1 Phylum         1     3   1.88      2
# 2 Genus          2    11   4.88      5
# 3 Species        3    16   7.38      7
# 4 ASV           53   300 187.      184
# 5 Missing Sp     0    93  19.8      14

numtaxa %>%
    group_by(taxon, age) %>%
    summarise(
        min = min(count), max = max(count), mean = mean(count),
        median = median(count)
    )

#   taxon      age     min   max   mean median
#   <fct>      <chr> <int> <int>  <dbl>  <dbl>
# 1 Phylum     early     1     3   1.77      2
# 2 Phylum     late      1     3   1.98      2
# 3 Genus      early     2    11   5.34      6
# 4 Genus      late      2     9   4.45      4
# 5 Species    early     4    16   7.99      8
# 6 Species    late      3    10   6.82      7
# 7 ASV        early    53   300 181.      168
# 8 ASV        late    115   296 193.      187
# 9 Missing Sp early     0    93  15.7      13
# 10 Missing Sp late      1    66  23.5      16


numtaxa_p <- numtaxa %>%
    ggplot(aes(x = count, fill = age)) +
    scale_fill_manual(values = agecol) +
    geom_density(alpha = 0.5, linewidth = 0.1) +
    facet_wrap(. ~ taxon, scales = "free", ncol = 5, nrow = 1) +
    theme_pubr(base_size = 5) +
    theme(
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    scale_x_continuous(n.breaks = 4) +
    xlab("Number of taxa")

plotsave(numtaxa_p, "./results/16S/EDA/numtaxa", width = 8, height = 3)

isoline_numtaxa_summary <- ps_df %>%
    filter(count > 0) %>%
    group_by(sampleID, Isoline) %>%
    summarise(
        Kingdom = length(unique(setdiff(Kingdom, NA))),
        Phylum = length(unique(setdiff(Phylum, NA))),
        Class = length(unique(setdiff(Class, NA))),
        Order = length(unique(setdiff(Order, NA))),
        Family = length(unique(setdiff(Family, NA))),
        Genus = length(unique(setdiff(Genus, NA))),
        `Missing Sp` = sum(is.na(Species)),
        Species = length(unique(setdiff(Species, NA))),
        ASV = length(unique(setdiff(ASV, NA)))
    ) %>%
    select(-Kingdom, -Class, -Order, -Family) %>%
    gather(key = "taxon", value = "count", -sampleID, -Isoline) %>%
    mutate(taxon = factor(taxon, levels = c(
        "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV",
        "Missing Sp"
    ))) %>%
    group_by(taxon, Isoline) %>%
    summarise(
        min = min(count), max = max(count), mean = mean(count),
        median = median(count)
    ) %>%
    select(-min, -max, -median) %>%
    group_by(taxon) %>%
    summarise(
        minx = min(mean),
        maxx = max(mean),
        meanx = mean(mean),
        medianx = median(mean),
        sdx = sd(mean)
    ) %>%
    mutate(`max/min` = maxx / minx) %>%
    set_names(c("Taxa Level", "Minimum", "Maximum", "Mean", "Median", "SD", "Max/Min"))

tablesave(isoline_numtaxa_summary, "./results/16S/EDA/isoline_numtaxa_summary")

#   `Taxa Level` Minimum Maximum   Mean Median     SD `Max/Min`
#   <fct>          <dbl>   <dbl>  <dbl>  <dbl>  <dbl>     <dbl>
# 1 Phylum          1.25     2.5   1.88   1.88  0.306      2
# 2 Genus           3.25     7.5   4.90   4.88  1.06       2.31
# 3 Species         5.67     9.5   7.40   7.39  1.05       1.68
# 4 ASV           139.     246.  187.   188.   25.0        1.78
# 5 Missing Sp      8       50.3  20.1   16.9   9.81       6.29


taxaprev <- ps_df %>%
    filter(count > 0) %>%
    select(-Kingdom, -Class, -Order, -Family) %>%
    gather(key = "taxon", value = "value", Species, Genus, Phylum, ASV) %>%
    group_by(taxon, age, value) %>%
    summarise(n = length(unique(sampleID))) %>%
    mutate(taxon = factor(taxon, levels = c(
        "Phylum", "Genus", "Species", "ASV"
    ))) %>%
    ggplot(aes(x = n, fill = age)) +
    scale_fill_manual(values = agecol) +
    geom_density(alpha = 0.5, linewidth = 0.1) +
    facet_wrap(. ~ taxon, scales = "free_y", ncol = 4, nrow = 1) +
    theme_pubr(base_size = 5) +
    theme(
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    xlab("Number of samples")

plotsave(taxaprev, "./results/16S/EDA/taxaprev", width = 8, height = 3)

ps_df %>%
    filter(count > 0) %>%
    select(-Kingdom, -Class, -Order, -Family) %>%
    gather(key = "taxon", value = "value", Species, Genus, Phylum, ASV) %>%
    group_by(taxon, value) %>%
    summarise(n = length(unique(sampleID))) %>%
    mutate(taxon = factor(taxon, levels = c(
        "Phylum", "Genus", "Species", "ASV"
    ))) %>%
    group_by(taxon) %>%
    summarise(
        min = min(n), max = max(n), mean = mean(n), median = median(n)
    )
#  taxon     min   max  mean median
#  <fct>   <int> <int> <dbl>  <dbl>
# 1 Phylum      5   161  75.5     68
# 2 Genus       1   161  30.4      9
# 3 Species     1   161  48       16
# 4 ASV         1   160  44.5     21

numtaxa_prev <- ggarrange(numtaxa_p, taxaprev,
    common.legend = TRUE, legend = "top",
    labels = "auto", font.label = list(size = 8)
)

plotsave(numtaxa_prev, "./results/16S/EDA/numtaxa_and_prev", width = 16, height = 3)

taxaprev_isolines <- ps_df %>%
    group_by(Isoline, ASV, Species, Genus, Phylum) %>%
    summarise(count = sum(count)) %>%
    filter(count > 0) %>%
    gather(key = "taxon", value = "value", Species, Genus, Phylum, ASV) %>%
    group_by(taxon, value) %>%
    summarise(n = length(unique(Isoline))) %>%
    mutate(taxon = factor(taxon, levels = c(
        "Phylum", "Genus", "Species", "ASV"
    ))) %>%
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 2) +
    facet_wrap(~taxon, scales = "free_y", ncol = 4, nrow = 2) +
    theme_pubr(base_size = 5) +
    theme(
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    xlab("Number of isolines") +
    ylab("Number of taxa")

plotsave(taxaprev_isolines, "./results/16S/EDA/taxaprev_isolines", width = 8, height = 3)

ps_genus <- tax_glom(ps, taxrank = "Genus")
genusmat <- otu_table(ps_genus)
genus_taxatable <- mutate(data.frame(tax_table(ps_genus)),
    ASV = colnames(genusmat)
)

data.frame(t(genusmat), ASV = colnames(genusmat)) %>%
    gather(key = "sampleID", value = "count", -ASV) %>%
    left_join(genus_taxatable) %>%
    left_join(sample_data(ps_genus)) %>%
    group_by(age, Isoline, Genus) %>%
    summarise(abd = mean(count)) %>%
    ungroup() %>%
    mutate(Genus = fct_lump_n(
        f = Genus, n = 7, w = abd,
        other_level = "Other"
    )) %>%
    mutate(Isoline = factor(Isoline, levels = intersect(1:100, Isoline))) %>%
    filter(age == "early") %>%
    head()

isolinemap <- read_delim("./data/isoline_map.csv", delim = ";") %>%
    set_names(c("IsolineID", "Isoline")) %>%
    mutate(
        Isoline = as.character(Isoline),
        IsolineID = as.character(IsolineID)
    )


genus_bar <- data.frame(t(genusmat), ASV = colnames(genusmat)) %>%
    gather(key = "sampleID", value = "count", -ASV) %>%
    left_join(genus_taxatable) %>%
    left_join(sample_data(ps_genus)) %>%
    group_by(age, Isoline, Genus) %>%
    summarise(abd = mean(count)) %>%
    ungroup() %>%
    mutate(Genus = fct_lump_n(
        f = Genus, n = 7, w = abd,
        other_level = "Other"
    )) %>%
    left_join(isolinemap) %>%
    mutate(IsolineID = factor(IsolineID, levels = intersect(1:1000, IsolineID))) %>%
    ggplot(aes(x = IsolineID, y = abd, fill = Genus)) +
    geom_bar(
        stat = "identity", width = 1, color = "gray25",
        position = "fill", linewidth = 0.1,
    ) +
    scale_fill_brewer(type = "qual", palette = 1) +
    facet_wrap(~age) +
    ylab("Mean Relative Abundance") +
    theme_minimal(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

plotsave(genus_bar, "./results/16S/EDA/genus_bar", width = 16, height = 6)

ps_clr <- microbiome::transform(ps, transform = "clr")
asvmat <- otu_table(ps_clr)
pcax <- prcomp(asvmat)
pcx <- data.frame(pcax$x, sampleID = rownames(pcax$x)) %>%
    left_join(sample_data(ps_clr))

sig_pcs <- which(summary(pcax)$imp[2, ] * 100 >= 5)
sigvals <- paste("PC", seq_len(ncol(pcax$x)), " (",
    round(summary(pcax)$imp[2, ] * 100, 2), "%)",
    sep = ""
)

pc12_mean <- select(pcx, age, PC1, PC2) %>%
    group_by(age) %>%
    summarise(PC1mean = mean(PC1), PC2mean = mean(PC2))

pcaplot <- full_join(pcx, pc12_mean) %>%
    ggplot(aes(x = PC1, y = PC2, color = age)) +
    geom_point(size = 0.5) +
    geom_point(aes(x = PC1mean, y = PC2mean), size = 2, shape = 18) +
    geom_segment(aes(x = PC1mean, y = PC2mean, xend = PC1, yend = PC2),
        alpha = 0.5, arrow = arrow(length = unit(0.05, "cm")), linewidth = 0.1
    ) +
    scale_color_manual(values = agecol) +
    theme_bw(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    xlab(sigvals[1]) +
    ylab(sigvals[2])

plotsave(pcaplot, "./results/16S/EDA/pcaplot", width = 8, height = 8)

pslog <- microbiome::transform(ps, transform = "log10p")

ordx <- ordinate(pslog, "PCoA", "bray")
p1 <- plot_ordination(pslog, ordx, type = "taxa", color = "Phylum") +
    theme_bw(base_size = 6) +
    theme(
        legend.position = "right",
        legend.key.size = unit(0.3, "cm")
    ) + guides(color = guide_legend(nrow = 2))
p2 <- plot_ordination(pslog, ordx, type = "samples", color = "age") +
    scale_color_manual(values = agecol) +
    theme_bw(base_size = 6) +
    theme(
        legend.position = "right",
        legend.key.size = unit(0.3, "cm")
    )
p3 <- plot_ordination(pslog, ordx, type = "samples", color = "run") +
    theme_bw(base_size = 6) +
    theme(
        legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
    )
pcoa <- ggarrange(p1, p2, p3,
    ncol = 1, nrow = 3, align = "hv", labels = "auto",
    font.label = list(size = 8)
)

plotsave(pcoa, "./results/16S/EDA/pcoa", width = 16, height = 21)

alphadiv <- estimate_richness(ps) %>%
    mutate(sampleID = sample_names(ps)) %>%
    select(sampleID, Observed, Shannon, Simpson, Chao1, InvSimpson) %>%
    gather(key = "metric", value = "diversity", -sampleID) %>%
    left_join(sample_data(ps))

alphadiv_plot <- alphadiv %>%
    ggplot(aes(x = age, y = diversity, fill = age)) +
    geom_violin(scale = "width", linewidth = 0.1) +
    geom_boxplot(
        width = 0.1, fill = "white", outlier.size = 0.05,
        linewidth = 0.25
    ) +
    scale_fill_manual(values = agecol) +
    facet_wrap(~metric, scales = "free_y", ncol = 5, nrow = 1) +
    stat_compare_means(
        aes(label = paste0("p = ", after_stat(p.format))),
        method = "wilcox", paired = FALSE, size = 5 / pntnorm,
        label.y.npc = 0.95
    ) +
    theme_pubr(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    ylab("Alpha diversity") +
    xlab("Age group")

plotsave(alphadiv_plot, "./results/16S/EDA/alphadiv_plot", width = 16, height = 4)

alphadiv_genus <- estimate_richness(tax_glom(ps, "Genus"),
    measures = c("Observed", "Shannon", "Simpson", "Chao1", "InvSimpson")
) %>%
    select(-se.chao1) %>%
    mutate(sampleID = sample_names(ps)) %>%
    gather(key = "metric", value = "diversity", -sampleID) %>%
    left_join(sample_data(ps))

alphadiv_genus_plot <- alphadiv_genus %>%
    ggplot(aes(x = age, y = diversity, fill = age)) +
    geom_violin(scale = "width", linewidth = 0.1) +
    geom_boxplot(
        width = 0.1, fill = "white", outlier.size = 0.05,
        linewidth = 0.25
    ) +
    scale_fill_manual(values = agecol) +
    facet_wrap(~metric, scales = "free_y", ncol = 5, nrow = 1) +
    stat_compare_means(
        aes(label = paste0("p = ", after_stat(p.format))),
        method = "wilcox", paired = FALSE, size = 5 / pntnorm,
        label.y.npc = 0.95
    ) +
    theme_pubr(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    ylab("Alpha diversity") +
    xlab("Age group")

plotsave(alphadiv_genus_plot, "./results/16S/EDA/alphadiv_genus_plot",
    width = 16, height = 4
)

alphadivplots <- ggarrange(alphadiv_plot + ggtitle("ASV level"),
    alphadiv_genus_plot + ggtitle("Genus level"),
    ncol = 1, nrow = 2, align = "hv",
    labels = "auto", font.label = list(size = 8), common.legend = TRUE,
    legend = "right"
)

plotsave(alphadivplots, "./results/16S/EDA/alphadiv_plots",
    width = 16, height = 8
)

# beta diversity
betadiv <- distance(ps, method = "bray")

betadiv_df <- reshape2::melt(as.matrix(betadiv)) %>%
    set_names(c("sample1", "sample2", "beta")) %>%
    mutate(sampleID = sample1) %>%
    left_join(select(data.frame(sample_data(ps)), sampleID, age)) %>%
    rename(age1 = age) %>%
    mutate(sampleID = sample2) %>%
    left_join(select(data.frame(sample_data(ps)), sampleID, age)) %>%
    rename(age2 = age) %>%
    mutate(
        sameage = age1 == age2,
        agecomp = paste(age1, age2, sep = "-")
    ) %>%
    filter(sample1 != sample2) %>%
    rowwise() %>%
    mutate(comppair = paste(sort(c(sample1, sample2)), collapse = "-")) %>%
    select(-sample1, -sample2, -sampleID) %>%
    unique() %>%
    mutate(agecomp = ifelse(agecomp == "late-early", "early-late", agecomp))

betadiv_plot <- betadiv_df %>%
    ggplot(aes(x = agecomp, y = beta, fill = agecomp)) +
    geom_violin(scale = "width", linewidth = 0.1) +
    geom_boxplot(
        width = 0.1, fill = "white", outlier.size = 0.05,
        linewidth = 0.25
    ) +
    scale_fill_manual(values = setNames(
        c(agecol[1], "#56316d", agecol[2]),
        c("early-early", "early-late", "late-late")
    )) +
    xlab("Age comparison") +
    ylab("Bray-Curtis dissimilarity") +
    stat_compare_means(
        comparisons = list(
            c("early-early", "early-late"),
            c("early-early", "late-late"),
            c("early-late", "late-late")
        ),
        size = 5 / pntnorm,
        label.y.npc = 0.95, method = "wilcox"
    ) +
    theme_pubr(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    ) +
    guides(fill = guide_legend(title = NULL))

plotsave(betadiv_plot, "./results/16S/EDA/betadiv_plot",
    width = 4, height = 4
)


ubiomedivplot <- ggarrange(
    ggarrange(numtaxa_p, taxaprev,
        common.legend = TRUE, legend = "top",
        labels = "auto", font.label = list(size = 8)
    ), ggarrange(
        genus_bar, pcaplot, betadiv_plot,
        ncol = 3, nrow = 1, widths = c(1.8, 1.2, 0.9),
        labels = c("c.", "d.", "e."), font.label = list(size = 8)
    ), ggarrange(alphadiv_plot + ggtitle("ASV level"),
        alphadiv_genus_plot + ggtitle("Genus level"),
        ncol = 2, nrow = 1, align = "hv",
        labels = c("f.", "g."), font.label = list(size = 8),
        common.legend = TRUE,
        legend = "none"
    ),
    ncol = 1, nrow = 3,
    heights = c(0.8, 1.5, 0.8)
)

plotsave(ubiomedivplot, "./results/16S/EDA/ubiomedivplot",
    width = 16,
    height = 10
)



x <- readRDS("processed/dada2_run2_trackreads.rds")
x


genus_bar$data %>%
    filter(Genus == "Wolbachia") %>%
    filter(abd != 0)



#####

# to test the statement "The observed increase in ASV diversity may result from
# the dominance of certain species in the gut of older flies, accompanied by
# greater variability among strains within these dominant species."


source("scripts/00-setup.R")

library(phyloseq)

ps <- readRDS("./processed/ps_pruned_rared10k.rds")


ps_sp <- tax_glom(ps, taxrank = "Species")


spnames <- apply(tax_table(ps), 1, function(x) paste(x, collapse = ";"))
spasvcnts <- sort(table(spnames), dec = T)

datmat <- as.matrix(otu_table(ps_sp))

e_l_means <- t(apply(datmat, 2, function(x) {
    tapply(x, as.factor(substr(rownames(datmat), 1, 1)), mean)
}))

rownames(e_l_means) <- unname(spnames[rownames(e_l_means)])


data.frame(e_l_means) %>%
    mutate(spname = rownames(.)) %>%
    left_join(data.frame(spname = names(spasvcnts), counts = spasvcnts)) %>%
    ggplot(aes(x = counts.Freq, y = log2(L + 1) - log2(E + 1))) +
    geom_point()


# check instead if the NA species are the reason

naasvs <- is.na(tax_table(ps)[, "Species"])
table(naasvs)

datmat <- as.matrix(otu_table(ps))
e_l_means <- t(apply(datmat, 2, function(x) {
    tapply(x, as.factor(substr(rownames(datmat), 1, 1)), mean)
}))

data.frame(e_l_means) %>%
    mutate(spname = rownames(.)) %>%
    left_join(data.frame(isna = unname(naasvs[, 1]), spname = rownames(naasvs))) %>%
    ggplot(aes(x = isna, y = log2(L + 1) - log2(E + 1))) +
    geom_boxplot()

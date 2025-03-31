source("scripts/00-setup.R")

library(phyloseq)

ps_all <- readRDS("./processed/ps_pruned_rared10k_100tables.rds")

diversity_vals <- lapply(ps_all, function(ps) {
    alphadiv <- estimate_richness(ps) %>%
        mutate(sampleID = sample_names(ps)) %>%
        select(sampleID, Observed, Shannon, Simpson, Chao1, InvSimpson) %>%
        gather(key = "metric", value = "diversity", -sampleID) %>%
        left_join(sample_data(ps))

    alphadiv_genus <- estimate_richness(tax_glom(ps, "Genus"),
        measures = c("Observed", "Shannon", "Simpson", "Chao1", "InvSimpson")
    ) %>%
        select(-se.chao1) %>%
        mutate(sampleID = sample_names(ps)) %>%
        gather(key = "metric", value = "diversity", -sampleID) %>%
        left_join(sample_data(ps))

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

    divmeasures <- list(
        alpha_asv = alphadiv, alpha_genus = alphadiv_genus,
        beta = betadiv_df
    )
})


div_asv_plot <- reshape2::melt(lapply(diversity_vals, function(x) {
    x$alpha_asv %>%
        group_by(age, metric) %>%
        summarise(
            mean_div = mean(diversity)
        ) %>%
        spread(key = age, value = mean_div) %>%
        mutate(dif = late - early)
}), id.vars = c("metric", "early", "late", "dif")) %>%
    rename(perm = L1) %>%
    ggplot(aes(x = metric, y = dif)) +
    geom_boxplot() +
    facet_wrap(~metric, scales = "free", ncol = 5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    xlab(NULL) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    theme_pubr(base_size = 5)

objsave(diversity_vals, "./data/processed/rarefaction_checks/diversity_vals")
plotsave(div_asv_plot, "./results/rarefaction_checks/div_asv_plot", width = 12, height = 4)

div_genus_plot <- reshape2::melt(lapply(diversity_vals, function(x) {
    x$alpha_genus %>%
        group_by(age, metric) %>%
        summarise(
            mean_div = mean(diversity)
        ) %>%
        spread(key = age, value = mean_div) %>%
        mutate(dif = late - early)
}), id.vars = c("metric", "early", "late", "dif")) %>%
    rename(perm = L1) %>%
    ggplot(aes(x = metric, y = dif)) +
    geom_boxplot() +
    facet_wrap(~metric, scales = "free", ncol = 5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    xlab(NULL) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    theme_pubr(base_size = 5)

plotsave(div_genus_plot, "./results/rarefaction_checks/div_genus_plot", width = 12, height = 4)
div_plot <- ggarrange(div_asv_plot, div_genus_plot, ncol = 1, nrow = 2)
plotsave(div_plot, "./results/rarefaction_checks/div_plot", width = 12, height = 8)

###########
library(DESeq2)
da_all <- lapply(ps_all, function(ps) {
    metadata <- sample_data(ps)
    metadata$DNAbatch <- as.numeric(as.character(metadata$DNAbatch))
    metadata$DNAbatch[is.na(metadata$DNAbatch)] <- 19
    metadata$DNAbatch <- factor(metadata$DNAbatch)
    sample_data(ps) <- metadata
    relab_early <- otu_table(microbiome::transform(subset_samples(ps, sample_data(ps)$age == "early"), "compositional"))
    relab_late <- otu_table(microbiome::transform(subset_samples(ps, sample_data(ps)$age == "late"), "compositional"))
    taxatokeep <- names(which((colSums(relab_early > 0) >= 5) & (colSums(relab_late > 0) >= 5) & (colMeans(relab_early) >= 1 / 10000) & (colMeans(relab_late) >= 1 / 10000)))
    ps <- prune_taxa(taxatokeep, ps)
    otu_table(ps) <- otu_table(ps) + 1
    psds <- phyloseq_to_deseq2(ps, ~ run + DNAbatch + age + Isoline)
    dseqres_isoline <- DESeq(psds, "LRT", reduced = ~ run + DNAbatch + age)
    dseqres_age <- DESeq(psds, "LRT", reduced = ~ run + DNAbatch + Isoline)
    res_age <- results(dseqres_age, cooksCutoff = FALSE)
    alpha <- 0.1
    sigtab_age <- res_age[which(res_age$padj < alpha), ]
    sigtab_age <- cbind(as(sigtab_age, "data.frame"), as(tax_table(ps)[rownames(sigtab_age), ], "matrix"))
    sigtab_age$ASV <- rownames(sigtab_age)
    rownames(sigtab_age) <- NULL
    res_iso <- results(dseqres_isoline, cooksCutoff = FALSE)
    alpha <- 0.1
    sigtab_iso <- res_iso[which(res_iso$padj < alpha), ]
    sigtab_iso <- cbind(as(sigtab_iso, "data.frame"), as(tax_table(ps)[rownames(sigtab_iso), ], "matrix"))
    sigtab_iso$ASV <- rownames(sigtab_iso)
    IsolineSignif <- prune_taxa(rownames(sigtab_iso), ps)
    IsolineSignif_species <- tax_glom(IsolineSignif, "Species")
    asvmat <- otu_table(IsolineSignif_species)
    allres <- reshape2::melt(asvmat) %>%
        set_names(c("sampleID", "ASV", "count")) %>%
        left_join(sample_data(IsolineSignif_species))
    allres <- select(allres, -F1quality, -AvCS, -CSlate, -CSexp)
    allres2 <- allres %>% gather(
        key = "LH", value = "LHvalue", -sampleID, -ASV,
        -count, -run, -DNAbatch, -Isoline, -age
    )
    allres2 <- allres2 %>%
        left_join(data.frame(tax_table(IsolineSignif_species),
            ASV = rownames(tax_table(IsolineSignif_species))
        )) %>%
        mutate(Species = ifelse(is.na(Species), "", Species)) %>%
        mutate(species = paste(Genus, Species, sep = " "))
    xx <- group_by(allres2, species, LH, age, Isoline) %>%
        na.omit() %>%
        summarise(counts = mean(count), LHvalue = unique(LHvalue)) %>%
        summarise(
            co_p = cor.test(LHvalue, counts, method = "spearman")$p.val,
            co = cor.test(LHvalue, counts, method = "spearman")$est
        ) %>%
        left_join(allres2) %>%
        filter(!(LH == "CSearly" & age == "late")) %>%
        filter(!(LH == "EarlyRS" & age == "late")) %>%
        mutate(padj = p.adjust(co_p, method = "fdr")) %>%
        filter(padj <= 0.1)
    xx <- xx %>%
        select(LH, age, co) %>%
        unique()
    list(
        psds = psds, dseqres_isoline = dseqres_isoline,
        dseqres_age = dseqres_age, sigtab_age = sigtab_age,
        sigtab_iso = sigtab_iso, lh_res = xx
    )
})
objsave(da_all, "./data/processed/rarefaction_checks/da_all")


age_da_list <- lapply(da_all, function(da) {
    da$sigtab_age %>%
        group_by(Genus, Species) %>%
        summarise(n = n())
})


reshape2::melt(age_da_list, id.vars = c("Genus", "Species", "n")) %>%
    group_by(Genus, Species) %>%
    summarise(
        nx = length(unique(L1)),
        mean_n = mean(n),
        min_n = min(n),
        max_n = max(n)
    ) %>%
    arrange(-nx)

#   Genus             Species          nx mean_n min_n max_n
#   <chr>             <chr>         <int>  <dbl> <int> <int>
# 1 Acetobacter       aceti           100  34.4     28    41
# 2 Acetobacter       indonesiensis   100  16.2      8    26
# 3 Acetobacter       tropicalis      100  47       47    47
# 4 Levilactobacillus NA              100   1.96     1     4
# 5 Ralstonia         pickettii       100  41.0     39    42
# 6 Pseudomonas       yamanorum        30   1.03     1     2
# 7 Acetobacter       persici           9   1        1     1


reshape2::melt(lapply(da_all, function(da) {
    da$lh_res
}), id.vars = c("species", "LH", "age", "co")) %>%
    mutate(direction = sign(co)) %>%
    group_by(species, LH, age, direction) %>%
    summarise(n = n()) %>%
    arrange(-n)

#   species                   LH         age   direction     n
#   <chr>                     <chr>      <chr>     <dbl> <int>
# 1 Acetobacter ascendens     AvLF       late         -1   100
# 2 Acetobacter indonesiensis Rsen       early         1   100
# 3 Acetobacter pasteurianus  AvLF       late         -1   100
# 4 Acetobacter persici       ActuarialB late         -1   100
# 5 Acetobacter tropicalis    ActuarialB late          1   100
# 6 Acetobacter ascendens     CSSlope    early         1    57
# 7 Pseudomonas yamanorum     Rsen       early        -1     7
# 8 Pseudomonas yamanorum     Rsen       late         -1     1

cors_iso <- diversity_vals[[1]][[1]] %>%
    group_by(metric, age, Isoline, AvLF, CSearly, CSSlope, EarlyRS, Rsen) %>%
    summarise(div = mean(diversity)) %>%
    rename(
        `Average Lifespan` = AvLF,
        `Early Reproductive Success` = EarlyRS,
        `Early Climbing Speed` = CSearly,
        `Functional Aging` = CSSlope,
        `Reproductive Senescence` = Rsen
    ) %>%
    gather(key = "LH", value = "LHvalue", -metric, -age, -Isoline, -div)


cors_iso_p <- cors_iso %>%
    # group_by(metric, age, LH) %>%
    # summarise(co = cor(div,LHvalue,method='spearman'),
    #         co_p = cor.test(div,LHvalue,method='spearman')$p.val) %>%
    # filter(!(LH == "CSearly" & age == "late")) %>%
    # filter(!(LH == "EarlyRS" & age == "late")) %>%
    # mutate(padj = p.adjust(co_p, method = "fdr")) %>%
    # filter(padj<=0.1) %>%
    # left_join(cors_iso) %>%
    ggplot(aes(x = LHvalue, y = div, color = age)) +
    geom_point(size = 0.1, color = "gray") +
    geom_smooth(method = "lm") +
    facet_wrap(LH ~ metric, scales = "free", ncol = 5) +
    stat_cor(method = "spearman", size = 5 / pntnorm, cor.coef = "rho") +
    scale_color_manual(values = agecol) +
    theme_pubr(base_size = 5)

plotsave(cors_iso_p, "./results/rarefaction_checks/cors_iso_one", width = 16, height = 16)

cors_iso %>%
    group_by(LH, metric, age) %>%
    summarise(
        co = cor(div, LHvalue, method = "spearman"),
        co_p = cor.test(div, LHvalue, method = "spearman")$p.val
    ) %>%
    mutate(padj = p.adjust(co_p, method = "fdr")) %>%
    filter(padj <= 0.1)

lh_asvdiv <- lapply(diversity_vals, function(x) {
    x[[1]] %>%
        group_by(metric, age, Isoline, AvLF, F1quality, AvCS, CSearly, CSlate, CSSlope, CSexp, EarlyRS, Rsen) %>%
        summarise(div = mean(diversity)) %>%
        gather(key = "LH", value = "LHvalue", -metric, -age, -Isoline, -div) %>%
        group_by(metric, age, LH) %>%
        summarise(
            co = cor(div, LHvalue, method = "spearman"),
            co_p = cor.test(div, LHvalue, method = "spearman")$p.val
        ) %>%
        filter(!(LH == "CSearly" & age == "late")) %>%
        filter(!(LH == "EarlyRS" & age == "late")) %>%
        mutate(padj = p.adjust(co_p, method = "fdr")) %>%
        filter(padj <= 0.1)
})

reshape2::melt(lh_asvdiv, id.vars = c("metric", "age", "LH", "co", "co_p", "padj")) %>%
    mutate(direction = sign(co)) %>%
    group_by(metric, age, LH, direction) %>%
    summarise(n = n())

#   metric     age   LH        direction     n
#   <chr>      <chr> <chr>         <dbl> <int>
# 1 InvSimpson late  F1quality         1    34
# 2 Shannon    late  F1quality         1   100
# 3 Simpson    late  F1quality         1    56


lh_genusdiv <- lapply(diversity_vals, function(x) {
    x[[2]] %>%
        group_by(metric, age, Isoline, AvLF, F1quality, AvCS, CSearly, CSlate, CSSlope, CSexp, EarlyRS, Rsen) %>%
        summarise(div = mean(diversity)) %>%
        gather(key = "LH", value = "LHvalue", -metric, -age, -Isoline, -div) %>%
        group_by(metric, age, LH) %>%
        summarise(
            co = cor(div, LHvalue, method = "spearman"),
            co_p = cor.test(div, LHvalue, method = "spearman")$p.val
        ) %>%
        filter(!(LH == "CSearly" & age == "late")) %>%
        filter(!(LH == "EarlyRS" & age == "late")) %>%
        mutate(padj = p.adjust(co_p, method = "fdr")) %>%
        filter(padj <= 0.1)
})

reshape2::melt(lh_genusdiv, id.vars = c("metric", "age", "LH", "co", "co_p", "padj")) %>%
    mutate(direction = sign(co)) %>%
    group_by(metric, age, LH, direction) %>%
    summarise(n = n())

#   metric age   LH    direction     n
# 1 Chao1  early Rsen         -1     1

lh_betadiv <- lapply(diversity_vals, function(x) {
    betadat <- x[["beta"]]
    pairwise.wilcox.test(betadat$beta, betadat$agecomp, p.adjust.method = "fdr")
})

all(sapply(lh_betadiv, function(x) {
    sum(x$p.val <= 0.1, na.rm = T) == 3
}))
# all comparisons have significant betadiv differences between groups (wilcox test, corrected by FDR).

source("scripts/00-setup.r")
library(phyloseq)
library(DESeq2)
ps <- readRDS("processed/ps_pruned_rared10k.rds")
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
objsave(psds, "./data/processed/DESeq2/dseq_obj")

dseqres_isoline <- DESeq(psds, "LRT", reduced = ~ run + DNAbatch + age)
objsave(dseqres_isoline, "./data/processed/DESeq2/dseqres_isoline")
dseqres_age <- DESeq(psds, "LRT", reduced = ~ run + DNAbatch + Isoline)
objsave(dseqres_age, "./data/processed/DESeq2/dseqres_age")

dseqres_isoline <- readRDS("./data/processed/DESeq2/dseqres_isoline.rds")
dseqres_age <- readRDS("./data/processed/DESeq2/dseqres_age.rds")

res <- results(dseqres_age, cooksCutoff = FALSE)
alpha <- 0.1
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
sigtab$ASV <- rownames(sigtab)
rownames(sigtab) <- NULL
sigtab %>%
    mutate(direction = sign(log2FoldChange)) %>%
    group_by(direction, Genus, Species) %>%
    summarise(n = length(unique(ASV)))

annotcol <- data.frame(sample_data(ps)[, 4])
annotrow <- select(data.frame(sigtab), Genus, Species, ASV) %>%
    mutate(Species = ifelse(is.na(Species), "", Species)) %>%
    mutate(Species = paste(Genus, Species, sep = " ")) %>%
    select(-Genus)
rownames(annotrow) <- annotrow$ASV
annotrow$ASV <- NULL
coldat <- list(age = agecol)
pheatmap::pheatmap(log2(1 + t(otu_table(prune_taxa(sigtab$ASV, ps)))),
    show_colnames = F,
    show_rownames = F, annotation_col = annotcol, annotation_row = annotrow,
    cellwidth = 3, cellheight = 3, fontsize = 10, border_color = NA,
    color = colorRampPalette(RColorBrewer::brewer.pal(8, "Purples"))(100),
    annotation_colors = coldat, cutree_rows = 4, cutree_cols = 5
)

pheatmap::pheatmap(log2(1 + t(otu_table(prune_taxa(sigtab$ASV, ps)))),
    show_colnames = F,
    show_rownames = F, annotation_col = annotcol, annotation_row = annotrow,
    cellwidth = 3, cellheight = 3, fontsize = 10, border_color = NA,
    color = colorRampPalette(RColorBrewer::brewer.pal(8, "Purples"))(100),
    annotation_colors = coldat, cutree_rows = 4, cutree_cols = 5,
    filename = "./results/DA/DA_age_heatmap.pdf"
)


res <- results(dseqres_isoline, cooksCutoff = FALSE)
alpha <- 0.1
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
unique(paste(sigtab[, "Genus"], sigtab[, "Species"]))
sigtab %>%
    mutate(ASV = rownames(sigtab)) %>%
    mutate(direction = sign(log2FoldChange)) %>%
    group_by(direction, Genus, Species) %>%
    summarise(n = length(unique(ASV))) %>%
    arrange(Genus, Species)
IsolineSignif <- prune_taxa(rownames(sigtab), ps)
# all are significant for isoline effect!!!




IsolineSignif_species <- tax_glom(IsolineSignif, "Species")
# IsolineSignif_species <- IsolineSignif
asvmat <- otu_table(IsolineSignif_species)
allres <- reshape2::melt(asvmat) %>%
    set_names(c("sampleID", "ASV", "count")) %>%
    left_join(sample_data(IsolineSignif_species))
allres <- select(allres, -F1quality, -AvCS, -CSlate, -CSexp, -ActuarialB)

allres2 <- allres %>%
    dplyr::rename(
        `Average Lifespan` = AvLF,
        `Early Reproductive Success` = EarlyRS,
        `Early Climbing Speed` = CSearly,
        `Functional Aging` = CSSlope,
        `Reproductive Senescence` = Rsen
    ) %>%
    gather(
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



plab <- xx %>%
    mutate(plabel = paste("rho=", round(co, 2), ", adj.p=", round(padj, 3), sep = "")) %>%
    group_by(plabel, age, species, LH) %>%
    summarise(LHvalue = min(LHvalue, na.rm = T), count = min(count, na.rm = T))

allsignif <- xx %>%
    ggplot(aes(y = count + 1, x = LHvalue, color = age)) +
    geom_point(size = 0.1, color = "gray") +
    geom_smooth(method = "lm", alpha = 0.2, linewidth = 0.5) +
    stat_summary(geom='point',  aes(group = Isoline), fun = "mean", size = 0.2) +
    facet_wrap(LH ~ species, scales = "free", ncol = 6, nrow = 1) +
    scale_y_log10() +
    geom_text(
        data = plab,
        aes(label = plabel), size = 5 / pntnorm, hjust = 0, vjust = 0,
        color = "gray25", fontface = "bold"
    ) +
    xlab("Life History Trait Value") +
    ylab("Abundance") +
    scale_color_manual(values = agecol) +
    guides(color = guide_legend(title = NULL)) +
    theme_pubr(base_size = 5) +
    theme(
        legend.key.size = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_rect(linewidth = 0.1)
    )
allsignif

plotsave(allsignif, "./results/DA/allsignif", width = 11, height = 3)

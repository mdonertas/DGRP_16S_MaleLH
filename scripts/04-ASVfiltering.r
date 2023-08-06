source("scripts/00-setup.R")
library(phyloseq)

ps <- readRDS("./processed/ps.rds")

# Filter based on ASV length
asvlengths <- nchar(colnames(otu_table(ps)))

data.frame(asv_length = asvlengths) %>%
    ggplot(aes(x = asv_length)) +
    geom_histogram()

summary(asvlengths)

# exclude everything below 400 and above 465
# lp = length pruned
ps_lp <- prune_taxa(asvlengths >= 400 & asvlengths <= 465, ps)

# remove asvs without phyla assignment
ps_lp_phy <- prune_taxa(c(!is.na(tax_table(ps_lp)[, "Phylum"])), ps_lp)

library(ggtree)
library(treeio)
# root the tree
phy_tree(ps_lp_phy) <- root(phy_tree(ps_lp_phy),
    sample(taxa_names(ps_lp_phy), 1),
    resolve.root = TRUE
)

# create a phylo4d object for plotting
p4_tree_b_qc <- phylobase::phylo4d(phy_tree(ps_lp_phy),
    tip.data = data.frame(tax_table(ps_lp_phy))
)

objsave(p4_tree_b_qc, "./processed/tree_before_treepruning")

# plot the tree with the cutoff for taxa pruning
tree_before_prune <- ggplot(p4_tree_b_qc, aes(x, y)) +
    geom_tree(color = "gray70", alpha = 0.7, size = 0.2) +
    geom_tippoint(aes(color = Phylum), size = 0.15) +
    geom_vline(
        xintercept = 3, linetype = "dotted", color = "darkred",
        linewidth = 0.3
    ) +
    annotate("text",
        x = 3.2, y = 700, label = "cutoff for taxa pruning",
        size = 6 / pntnorm, hjust = 0, color = "darkred"
    ) +
    theme_void(base_size = 6) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("Before pruning")
tree_before_prune

# save the tree plot
plotsave(tree_before_prune, "./results/EDA/phylotree_before_treepruning",
    width = 8, height = 12
)

# prune the tree for divergent taxa
taxa2prune <- (ggtree::fortify(p4_tree_b_qc) %>%
    filter(x >= 3) %>% # nolint: indentation_linter.
    filter(isTip))$label

# tp = tree pruned
ps_lp_phy_tp <- prune_taxa(!taxa_names(ps_lp_phy) %in% taxa2prune, ps_lp_phy)

# root the tree after pruning
phy_tree(ps_lp_phy_tp) <- root(phy_tree(ps_lp_phy_tp),
    sample(taxa_names(ps_lp_phy_tp), 1),
    resolve.root = TRUE
)

# create a phylo4d object for plotting
p4_tree_a_qc <- phylobase::phylo4d(phy_tree(ps_lp_phy_tp),
    tip.data = data.frame(tax_table(ps_lp_phy_tp))
)

# plot the tree after pruning
tree_after_prune <- ggplot(p4_tree_a_qc, aes(x, y)) +
    geom_tree(color = "gray70", alpha = 0.7, size = 0.2) +
    geom_tippoint(aes(color = Phylum), size = 0.15) +
    theme_void(base_size = 6) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("After pruning")
tree_after_prune

plotsave(tree_after_prune, "./results/EDA/phylotree_after_treepruning",
    width = 8, height = 12
)

# save pruned phyloseq objects after each step
objsave(ps_lp, "./processed/ps_lengthpruned")
objsave(ps_lp_phy, "./processed/ps_length_phylapruned")
objsave(ps_lp_phy_tp, "./processed/ps_length_phyla_treepruned")

objsave(p4_tree_a_qc, "./processed/tree_after_treepruning")

tree_before_after <- ggarrange(tree_before_prune + ggtitle(NULL),
    tree_after_prune + ggtitle(NULL),
    ncol = 2, nrow = 1, common.legend = TRUE, legend = "right",
    labels = "auto", font.label = list(size = 8)
)

plotsave(tree_before_after, "./results/EDA/tree_before_and_after",
    width = 16,
    height = 12
)

# apply prevelance filtering - filter out all taxa that are not present in at
# least 1% of samples
ps_lp_phy_tp_prev <- prune_taxa(names(which(colSums(otu_table(
    ps_lp_phy_tp
) > 0) >= 0.01)), ps_lp_phy_tp)

summary(colMeans(otu_table(ps_lp_phy_tp_prev) > 0))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01130 0.02825 0.11299 0.27309 0.45904 0.98305

# save the phyloseq object after prevalence filtering
objsave(ps_lp_phy_tp_prev, "./processed/ps_length_phyla_treepruned_prev")

# create exploratory plots for ASV length, taxa assignment rate, and tree before
# and after each pruning step

asvlength_before <- data.frame(asvlen = nchar(colnames(otu_table(ps)))) %>%
    ggplot(aes(x = asvlen)) +
    geom_histogram() +
    xlab("ASV length")

plotsave(asvlength_before, "./results/EDA/asvlength_before_qc", height = 4)

asvlength_after <- data.frame(asvlen = nchar(colnames(otu_table(ps_lp)))) %>%
    ggplot(aes(x = asvlen)) +
    geom_histogram() +
    xlab("ASV length")

plotsave(asvlength_after, "./results/EDA/asvlength_after_qc", height = 4)

asvlength_before_after <- ggarrange(asvlength_before,
    asvlength_after,
    ncol = 2, nrow = 1,
    labels = "auto", font.label = list(size = 8)
)

plotsave(asvlength_before_after, "./results/EDA/asvlength_before_and_after",
    width = 16,
    height = 4
)

assignrate_before <- tibble(
    assignrate = colMeans(!is.na(tax_table(ps_lp))),
    taxalevel = factor(colnames(tax_table(ps_lp)),
        levels = colnames(tax_table(ps_lp))
    )
) %>%
    ggplot(aes(x = taxalevel, y = assignrate)) +
    geom_bar(stat = "identity") +
    xlab("Taxanomy Level") +
    ylab("Assignment Rate")

plotsave(assignrate_before, "./results/EDA/taxaassign_before_qc", height = 4)

assignrate_after <- tibble(
    assignrate = colMeans(!is.na(tax_table(ps_lp_phy))),
    taxalevel = factor(colnames(tax_table(ps_lp_phy)),
        levels = colnames(tax_table(ps_lp_phy))
    )
) %>%
    ggplot(aes(x = taxalevel, y = assignrate)) +
    geom_bar(stat = "identity") +
    xlab("Taxanomy Level") +
    ylab("Assignment Rate")

plotsave(assignrate_after, "./results/EDA/taxaassign_after_qc", height = 4)

assignrate_before_after <- ggarrange(assignrate_before,
    assignrate_after,
    ncol = 2, nrow = 1,
    labels = "auto", font.label = list(size = 8)
)

asvlength_end <- data.frame(asvlen = nchar(
    colnames(otu_table(ps_lp_phy_tp_prev))
)) %>%
    ggplot(aes(x = asvlen)) +
    geom_histogram() +
    xlab("ASV length")

assignrate_end <- tibble(
    assignrate = colMeans(!is.na(tax_table(ps_lp_phy_tp_prev))),
    taxalevel = factor(colnames(tax_table(ps_lp_phy_tp_prev)),
        levels = colnames(tax_table(ps_lp_phy_tp_prev))
    )
) %>%
    ggplot(aes(x = taxalevel, y = assignrate)) +
    geom_bar(stat = "identity") +
    xlab("Taxanomy Level") +
    ylab("Assignment Rate")

p4_tree_end <- phylobase::phylo4d(phy_tree(ps_lp_phy_tp_prev),
    tip.data = data.frame(tax_table(ps_lp_phy_tp_prev))
)

# plot the tree after pruning
tree_end <- ggplot(p4_tree_end, aes(x, y)) +
    geom_tree(color = "gray70", alpha = 0.7, size = 0.2) +
    geom_tippoint(aes(color = Phylum), size = 0.15) +
    theme_void(base_size = 6) +
    xlab(NULL) +
    ylab(NULL)
tree_end

plotsave(tree_end, "./results/EDA/phylotree_end",
    width = 8, height = 12
)

qc_summary_plot <- ggarrange(
    ggarrange(asvlength_before, asvlength_after,
        ncol = 2, nrow = 1,
        labels = c("a.", "b."), font.label = list(size = 8)
    ),
    ggarrange(assignrate_before, assignrate_after,
        ncol = 2, nrow = 1,
        labels = c("c.", "d."), font.label = list(size = 8)
    ),
    ggarrange(tree_before_prune + ggtitle(NULL),
        tree_after_prune + ggtitle(NULL),
        ncol = 2, nrow = 1,
        common.legend = TRUE, legend = "bottom",
        labels = c("e.", "f."), font.label = list(size = 8)
    ),
    ggarrange(tree_end, ggarrange(asvlength_end, assignrate_end,
        ncol = 1, nrow = 2,
        labels = c("h.", "i."), font.label = list(size = 8)
    ), ncol = 2, nrow = 1, labels = c("g.", ""), font.label = list(size = 8)),
    ncol = 1, nrow = 4, heights = c(1, 1, 1.75, 1.75), align = "hv"
)

plotsave(qc_summary_plot, "./results/EDA/asvqc_summary_plot",
    height = 20,
    width = 18
)

pruned_taxa <- setdiff(
    colnames(otu_table(ps)),
    colnames(otu_table(ps_lp_phy_tp_prev))
)

filtered_ps <- prune_taxa(taxa_names(ps) %in% pruned_taxa, ps)

filtered_ps_info <- data.frame(
    count = sample_sums(filtered_ps),
    num_filtered = rowSums(otu_table(filtered_ps) > 0),
    sampleID = sample_names(filtered_ps)
) %>%
    left_join(sample_data(filtered_ps))

# plot the number of reads and ASVs filtered per batch
tot_read <- ggplot(filtered_ps_info, aes(
    x = DNAbatch, y = count,
    color = run
)) +
    geom_boxplot(outlier.size = 0.2) +
    ylab("Total reads") +
    xlab("DNA Extraction Batch") +
    ggtitle("Total # reads filtered")

tot_asvs <- ggplot(filtered_ps_info, aes(
    x = DNAbatch, y = num_filtered,
    color = run
)) +
    geom_boxplot(outlier.size = 0.2) +
    ylab("Number of ASVs") +
    xlab("DNA Extraction Batch") +
    ggtitle("# of filtered ASVs")

plotsave(tot_read, "./results/EDA/prunedreads_per_batch")

plotsave(tot_asvs, "./results/EDA/prunedasvcount_per_batch")

filtered_asv_per_batch <- ggarrange(tot_read, tot_asvs,
    ncol = 2, nrow = 1, common.legend = TRUE, legend = "right",
    labels = "auto", font.label = list(size = 8)
)
plotsave(filtered_asv_per_batch, "./results/EDA/filtered_asv_per_batch",
    width = 16
)

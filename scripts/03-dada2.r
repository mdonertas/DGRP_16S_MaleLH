source("scripts/00-setup.R")
library(dada2)

metadata <- readRDS("./processed/metadata.rds")

run1_path <- "./raw/RUNGEN52_2019/fastq"
run2_path <- "./raw/RUNGEN66/fastq"

run1_fn_fs <- sort(list.files(run1_path,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
))
run1_fn_rs <- sort(list.files(run1_path,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
))

run2_fn_fs <- sort(list.files(run2_path,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
))
run2_fn_rs <- sort(list.files(run2_path,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
))


run1_samples <- sapply(strsplit(basename(run1_fn_fs), "_"), `[`, 1)
run2_samples <- sapply(strsplit(basename(run2_fn_fs), "_"), `[`, 1)

run1ind <- run1_samples %in% metadata$sampleID
run2ind <- run2_samples %in% metadata$sampleID

run1_fn_fs <- run1_fn_fs[run1ind]
run1_fn_rs <- run1_fn_rs[run1ind]

run2_fn_fs <- run2_fn_fs[run2ind]
run2_fn_rs <- run2_fn_rs[run2ind]

run1_samples <- run1_samples[run1ind]
run2_samples <- run2_samples[run2ind]

plotQualityProfile(run1_fn_fs[1:2])
plotQualityProfile(run2_fn_fs[1:2])

plotQualityProfile(run1_fn_rs[1:2])
plotQualityProfile(run2_fn_rs[1:2])

sapply(seq(0, length(run1_fn_fs), by = 6), function(n) {
    print(n)
    if (n + 6 <= length(run1_fn_fs)) {
        qualp_run1_fs <- plotQualityProfile(run1_fn_fs[(n + 1):(n + 6)])
    } else {
        qualp_run1_fs <- plotQualityProfile(run1_fn_fs[(n + 1):(
            length(run1_fn_fs))]) # nolint: indentation_linter.
    }
    plotsave(qualp_run1_fs, paste("./results/dada2/qualityProfiles/run1_F", n,
        sep = "" # nolint: indentation_linter.
    ), width = 16, height = 12)
})

sapply(seq(0, length(run2_fn_fs), by = 6), function(n) {
    print(n)
    if (n + 6 <= length(run2_fn_fs)) {
        qualp_run2_fs <- plotQualityProfile(run2_fn_fs[(n + 1):(n + 6)])
    } else {
        qualp_run2_fs <- plotQualityProfile(run2_fn_fs[(n + 1):(
            length(run2_fn_fs))]) # nolint: indentation_linter.
    }
    plotsave(qualp_run2_fs, paste("./results/dada2/qualityProfiles/run2_F", n,
        sep = "" # nolint: indentation_linter.
    ), width = 16, height = 12)
})

sapply(seq(0, length(run1_fn_rs), by = 6), function(n) {
    print(n)
    if (n + 6 <= length(run1_fn_rs)) {
        qualp_run1_rs <- plotQualityProfile(run1_fn_rs[(n + 1):(n + 6)])
    } else {
        qualp_run1_rs <- plotQualityProfile(run1_fn_rs[(n + 1):(
            length(run1_fn_rs))]) # nolint: indentation_linter.
    }
    plotsave(qualp_run1_rs, paste("./results/dada2/qualityProfiles/run1_R", n,
        sep = "" # nolint: indentation_linter.
    ), width = 16, height = 12)
})

sapply(seq(0, length(run2_fn_rs), by = 6), function(n) {
    print(n)
    if (n + 6 <= length(run2_fn_rs)) {
        qualp_run2_rs <- plotQualityProfile(run2_fn_rs[(n + 1):(n + 6)])
    } else {
        qualp_run2_rs <- plotQualityProfile(run2_fn_rs[(n + 1):(
            length(run2_fn_rs))]) # nolint: indentation_linter.
    }
    plotsave(qualp_run2_rs, paste("./results/dada2/qualityProfiles/run2_R", n,
        sep = "" # nolint: indentation_linter.
    ), width = 16, height = 12)
})


run1_filt_fs <- file.path(
    "./processed/RUNGEN52_2019", "filtered",
    paste0(run1_samples, "_F_filt.fastq.gz")
)
run1_filt_rs <- file.path(
    "./processed/RUNGEN52_2019", "filtered",
    paste0(run1_samples, "_R_filt.fastq.gz")
)

run2_filt_fs <- file.path(
    "./processed/RUNGEN66", "filtered",
    paste0(run2_samples, "_F_filt.fastq.gz")
)
run2_filt_rs <- file.path(
    "./processed/RUNGEN66", "filtered",
    paste0(run2_samples, "_R_filt.fastq.gz")
)

names(run1_filt_fs) <- run1_samples
names(run1_filt_rs) <- run1_samples

names(run2_filt_fs) <- run2_samples
names(run2_filt_rs) <- run2_samples

run1_out <- filterAndTrim(run1_fn_fs, run1_filt_fs, run1_fn_rs, run1_filt_rs,
    truncLen = c(280, 250),
    maxN = 0, rm.phix = TRUE,
    compress = TRUE, multithread = TRUE
)
head(run1_out)

summary(run1_out[, 2] / run1_out[, 1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9401  0.9942  0.9983  0.9936  0.9987  0.9990

# there are 3base long N in the beginning of the reads in run 2, reverse reads,
# so we trim 3 bases from left (without trimming, dada2 doesn't like Ns and
# filters out almost half of the reads)

run2_out <- filterAndTrim(run2_fn_fs, run2_filt_fs, run2_fn_rs, run2_filt_rs,
    truncLen = c(280, 250), trimLeft = c(0, 3),
    maxN = 0, rm.phix = TRUE,
    compress = TRUE, multithread = TRUE
)
head(run2_out)

summary(run2_out[, 2] / run2_out[, 1])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9413  0.9979  0.9989  0.9958  0.9993  1.0000


# to have comparable ASVs and make merging easier we apply the same trimming to
# run 1
run1_out <- filterAndTrim(run1_fn_fs, run1_filt_fs, run1_fn_rs, run1_filt_rs,
    truncLen = c(280, 250), trimLeft = c(0, 3),
    maxN = 0, rm.phix = TRUE,
    compress = TRUE, multithread = TRUE
)
head(run1_out)

summary(run1_out[, 2] / run1_out[, 1])

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9401  0.9942  0.9983  0.9936  0.9987  0.9990

run1_err_f <- learnErrors(run1_filt_fs, multithread = TRUE)
run1_err_r <- learnErrors(run1_filt_rs, multithread = TRUE)

run2_err_f <- learnErrors(run2_filt_fs, multithread = TRUE)
run2_err_r <- learnErrors(run2_filt_rs, multithread = TRUE)

pdf("./results/errorProfiles.pdf")
plotErrors(run1_err_f, nominalQ = TRUE)
plotErrors(run1_err_r, nominalQ = TRUE)

plotErrors(run2_err_f, nominalQ = TRUE)
plotErrors(run2_err_r, nominalQ = TRUE)
dev.off()

run1_dada_fs <- dada(run1_filt_fs, err = run1_err_f, multithread = TRUE)
run1_dada_rs <- dada(run1_filt_rs, err = run1_err_r, multithread = TRUE)

run2_dada_fs <- dada(run2_filt_fs, err = run2_err_f, multithread = TRUE)
run2_dada_rs <- dada(run2_filt_rs, err = run2_err_r, multithread = TRUE)

run1_mergers <- mergePairs(run1_dada_fs, run1_filt_fs, run1_dada_rs,
    run1_filt_rs,
    verbose = TRUE
)

run2_mergers <- mergePairs(run2_dada_fs, run2_filt_fs, run2_dada_rs,
    run2_filt_rs,
    verbose = TRUE
)

run1_seqtab <- makeSequenceTable(run1_mergers)
run2_seqtab <- makeSequenceTable(run2_mergers)

table(nchar(getSequences(run1_seqtab)))

run1_asv_len_p <- data.frame(asv_length = nchar(getSequences(run1_seqtab))) %>%
    ggplot(aes(x = asv_length)) +
    geom_histogram()

plotsave(run1_asv_len_p, "./results/dada2/merged_seq/asv_len_run1")

table(nchar(getSequences(run2_seqtab)))

run2_asv_len_p <- data.frame(asv_length = nchar(getSequences(run2_seqtab))) %>%
    ggplot(aes(x = asv_length)) +
    geom_histogram()

plotsave(run2_asv_len_p, "./results/dada2/merged_seq/asv_len_run2")

# both runs have 2 peaks at length 432 and 462

seqtab <- mergeSequenceTables(run1_seqtab, run2_seqtab)
table(nchar(getSequences(seqtab)))

asv_len_p <- data.frame(asv_length = nchar(getSequences(seqtab))) %>%
    ggplot(aes(x = asv_length)) +
    geom_histogram()

plotsave(asv_len_p, "./results/dada2/merged_seq/asv_len_bothruns")

seqtab_nochim <- removeBimeraDenovo(seqtab,
    method = "consensus",
    multithread = TRUE, verbose = TRUE
)
dim(seqtab_nochim)
# 177 2020

ncol(seqtab_nochim) / ncol(seqtab)
# [1] 0.2956674 # only 30% of the ASV sequences are maintained after chimera rem

sum(seqtab_nochim) / sum(seqtab)
# [1] 0.7931039 # 79.3% of reads are not chimeric.

get_n <- function(x) sum(getUniques(x))
track_run1 <- cbind(
    run1_out, sapply(run1_dada_fs, get_n),
    sapply(run1_dada_rs, get_n),
    sapply(run1_mergers, get_n)
)

colnames(track_run1) <- c(
    "input", "filtered", "denoisedF", "denoisedR",
    "merged"
)
rownames(track_run1) <- run1_samples

track_run2 <- cbind(
    run2_out, sapply(run2_dada_fs, get_n),
    sapply(run2_dada_rs, get_n),
    sapply(run2_mergers, get_n)
)
colnames(track_run2) <- c(
    "input", "filtered", "denoisedF", "denoisedR",
    "merged"
)
rownames(track_run2) <- run2_samples

objsave(track_run1, "./processed/dada2_run1_trackreads")
objsave(track_run2, "./processed/dada2_run2_trackreads")

data.frame(track_run1) %>%
    mutate(
        filtperc = filtered / input,
        mergeperc = merged / filtered
    ) %>%
    summary()

#     input           filtered        denoisedF        denoisedR
# Min.   : 19038   Min.   : 18753   Min.   : 17679   Min.   : 18124
# 1st Qu.: 87368   1st Qu.: 87238   1st Qu.: 85582   1st Qu.: 85922
# Median : 97950   Median : 96073   Median : 93746   Median : 94124
# Mean   :100294   Mean   : 99750   Mean   : 97508   Mean   : 97920
# 3rd Qu.:119032   3rd Qu.:118909   3rd Qu.:116451   3rd Qu.:117460
# Max.   :157279   Max.   :156903   Max.   :154221   Max.   :154512

#     merged          filtperc        mergeperc
# Min.   : 17183   Min.   :0.9401   Min.   :0.8728
# 1st Qu.: 83224   1st Qu.:0.9942   1st Qu.:0.9442
# Median : 91716   Median :0.9983   Median :0.9527
# Mean   : 94805   Mean   :0.9936   Mean   :0.9481
# 3rd Qu.:113233   3rd Qu.:0.9987   3rd Qu.:0.9598
# Max.   :150855   Max.   :0.9990   Max.   :0.9687

# run1 stats are pretty good.

data.frame(track_run2) %>%
    mutate(
        filtperc = filtered / input,
        mergeperc = merged / filtered
    ) %>%
    summary()

#     input            filtered         denoisedF         denoisedR
# Min.   :    131   Min.   :    130   Min.   :    104   Min.   :     96
# 1st Qu.:  49836   1st Qu.:  49792   1st Qu.:  49145   1st Qu.:  48730
# Median :  63823   Median :  63756   Median :  62783   Median :  62295
# Mean   :  91555   Mean   :  91359   Mean   :  89998   Mean   :  89453
# 3rd Qu.: 100391   3rd Qu.: 100308   3rd Qu.:  99248   3rd Qu.:  98618
# Max.   :1048436   Max.   :1047510   Max.   :1038435   Max.   :1028840

#       merged          filtperc        mergeperc
#   Min.   :    88   Min.   :0.9413   Min.   :0.6769
#   1st Qu.: 47121   1st Qu.:0.9979   1st Qu.:0.9170
#   Median : 60261   Median :0.9989   Median :0.9398
#   Mean   : 85598   Mean   :0.9958   Mean   :0.9219
#   3rd Qu.: 92530   3rd Qu.:0.9993   3rd Qu.:0.9539
#   Max.   :999271   Max.   :1.0000   Max.   :0.9696

# in run 2 we lose some of the reads for some samples during merging but on
# average it is not bad

taxa <- assignTaxonomy(seqtab_nochim,
    "./raw/SILVA/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
    multithread = TRUE
)

taxa_print <- taxa # Removing sequence rownames for display only
rownames(taxa_print) <- NULL

colMeans(!is.na(taxa_print))
#  Kingdom    Phylum     Class     Order    Family     Genus   Species
# 0.9757426 0.8544554 0.8247525 0.8242574 0.8059406 0.7396040 0.4935644

# this is a pretty good result, suggesting 49.5% of our ASVs can be assigned to
# species and 85% of all ASVs have phyla level annotation. We will probably
# filter out the unassigned 13% at the phyla level after checking the tree.

library(DECIPHER)
library(phangorn)
seqs <- colnames(seqtab_nochim)
names(seqs) <- seqs
align <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
phang_align <- phyDat(as(align, "matrix"), type = "DNA")
dml <- dist.ml(phang_align)
nj_tree <- NJ(dml) # Note, tip order != sequence order
fit <- pml(nj_tree, data = phang_align)
fit_gtr <- update(fit, k = 4, inv = 0.2)
fit_gtr <- optim.pml(fit_gtr,
    model = "GTR", optInv = TRUE, optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
)
otumat <- seqtab_nochim

rownames(metadata) <- metadata$sampleID

library(phyloseq)
ps <- phyloseq(
    tax_table(taxa), sample_data(metadata),
    otu_table(otumat, taxa_are_rows = FALSE), phy_tree(fit_gtr$tree)
)

phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1),
    resolve.root = TRUE
)

objsave(ps, "./processed/ps")

save(list = ls(), file = "./processed/dada2_allobjects.RData")

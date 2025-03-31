source("scripts/00-setup.r")

library(caret)

ps <- readRDS("processed/ps_pruned_rared10k.rds")

ps <- tax_glom(ps, taxrank = "Species")

ps <- microbiome::transform(ps, "clr")
metadata <- data.frame(sample_data(ps))


df <- as.data.frame((otu_table(ps)))

df$age <- metadata[rownames(df), "age"]

age <- metadata %>%
    select(age, Isoline) %>%
    unique()

age <- setNames(age$age, age$Isoline)
isolines <- setNames(metadata$Isoline, metadata$sampleID)

train_isolines <- sample(unique(isolines), 15, replace = FALSE)
test_isolines <- setdiff(unique(isolines), train_isolines)

train_samples <- names(isolines[isolines %in% train_isolines])
test_samples <- names(isolines[isolines %in% test_isolines])

training <- df[train_samples, ]
testing <- df[test_samples, ]

fit_ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 5
)

rf <- train(age ~ .,
    data = training,
    method = "rf",
    trControl = fit_ctrl,
    verbose = FALSE,
    tuneLength = 25
)

test_pred <- predict(rf, testing)
conf <- confusionMatrix(test_pred, factor(testing$age))
conf$overall["Accuracy"]

objsave(rf, "processed/age_pred/rf")
objsave(test_pred, "processed/age_pred/test_predictions")
objsave(testing, "processed/age_pred/test_data")
objsave(training, "processed/age_pred/training_data")
objsave(conf, "processed/age_pred/confusion_matrix")

conf_matrix <- data.frame(actual_class = testing$age, predicted_class = test_pred, val = 1) %>%
    group_by(actual_class, predicted_class) %>%
    summarise(count = sum(val)) %>%
    ggplot(aes(x = actual_class, fill = predicted_class, y = count)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(aes(group = predicted_class, label = round(count, 2)), position = position_fill(vjust = 0.5), color = "white", size = 6 / pntnorm) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Actual Age", y = NULL, fill = "Predicted Age") +
    scale_fill_manual(values = agecol) +
    theme_pubr(base_size = 5) +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.1, "cm")
    )


plotsave(conf_matrix, "results/age_pred/confusion_matrix", width = 4, height = 4)




mostimportant_5feat <- rownames(varImp(rf)$importance %>%
    mutate(ASV = rownames(.)) %>%
    arrange(-Overall) %>%
    head(5))

mostimportant_5feat_names <- unname(apply(tax_table(ps)[mostimportant_5feat, ], 1, function(x) {
    paste(x[!is.na(x)], collapse = ";")
}))

mostimportant_5feat_names

# [1] "Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;tropicalis"
# [2] "Bacteria;Proteobacteria;Gammaproteobacteria;Burkholderiales;Burkholderiaceae;Ralstonia;pickettii"
# [3] "Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;indonesiensis"
# [4] "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Leuconostoc;pseudomesenteroides"
# [5] "Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;persici"

rf <- readRDS("processed/age_pred/rf.rds")
test_pred <- readRDS("processed/age_pred/test_predictions.rds")
testing <- readRDS("processed/age_pred/test_data.rds")
training <- readRDS("processed/age_pred/training_data.rds")
conf <- readRDS("processed/age_pred/confusion_matrix.rds")

mostimportant_5feat <- rownames(varImp(rf)$importance %>% arrange(-Overall) %>% head(5))
library(phyloseq)
varimp_plot <- varImp(rf)$importance %>%
    mutate(ASV = rownames(.)) %>%
    left_join((data.frame(tax_table(ps)) %>% mutate(ASV = rownames(.)))) %>%
    arrange(-Overall) %>%
    select(Genus, Species, Overall) %>%
    mutate(Species = paste(Genus, Species, sep = " ")) %>%
    select(-Genus) %>%
    head(5) %>%
    ggplot(aes(x = reorder(Species, -Overall), y = Overall)) +
    geom_bar(stat = "identity") +
    theme_pubr(base_size = 5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Species") +
    ylab("Feature Importance")

plotsave(varimp_plot, "results/age_pred/varimp_plot", width = 3, height = 8)


## what are those late point samples that are being predicted as early?

rownames(testing[which(as.character(test_pred) != testing$age), ])
# Isolines 11 (4 samples) and 15 (3 samples) were predicted as early time points
# although they are late time point. Is there something special with these
# Isolines? One sample from Isoline 22 was predicted as late point, although it
# is an early time point.

isolinedata <- data.frame(sample_data(ps)) %>%
    select(Isoline, AvLF, CSearly, CSSlope, EarlyRS, Rsen, F1quality) %>%
    unique() %>%
    mutate(relAge = 25 / AvLF)

isolinedata_all <- isolinedata %>%
    filter(!Isoline %in% c(11, 15, 22)) %>%
    select(-Isoline)

specialIsoRes <- t(apply(isolinedata %>%
    filter(Isoline %in% c(11, 15, 22)) %>%
    select(-Isoline), 1, function(x) {
    sapply(seq_len(ncol(isolinedata_all)), function(i) {
        min(c(
            mean(isolinedata_all[, i] >= x[[i]]),
            mean(isolinedata_all[, i] <= x[[i]])
        ))
    })
}))
colnames(specialIsoRes) <- colnames(isolinedata_all)
specialIsoRes
#           AvLF   CSearly   CSSlope   EarlyRS      Rsen F1quality    relAge
# E11A 0.3157895 0.2631579 0.3157895 0.3157895 0.3157895 0.4210526 0.3157895
# E15A 0.2105263 0.0000000 0.1578947 0.4736842 0.3157895 0.1578947 0.2105263
# E22A 0.3157895 0.4210526 0.4210526 0.4210526 0.2105263 0.2105263 0.3157895

# It doesn't seem like there is anything special about these Isolines.
# Let's check the distance between young and old samples of the same Isoline.

early_late_beta_p <- reshape2::melt(as.matrix(distance(ps, method = "euclidean"))) %>%
    filter(Var1 != Var2) %>%
    mutate(sampleID = Var1) %>%
    left_join(sample_data(ps)) %>%
    mutate(Isoline1 = Isoline, age1 = age) %>%
    select(Var1, Var2, Isoline1, age1, value) %>%
    mutate(sampleID = Var2) %>%
    left_join(sample_data(ps)) %>%
    mutate(Isoline2 = Isoline, age2 = age) %>%
    select(Var1, Var2, Isoline1, age1, Isoline2, age2, value) %>%
    filter(Isoline1 == Isoline2) %>%
    mutate(sameAge = age1 == age2) %>%
    ggplot(aes()) +
    ggridges::geom_density_ridges2(aes(x = value, fill = sameAge, y = Isoline1), alpha = 0.5) +
    Rvislib::theme_rvis(base_size = 8)



plotsave(early_late_beta_p, "results/age_pred/early_late_beta_p", width = 8, height = 12)

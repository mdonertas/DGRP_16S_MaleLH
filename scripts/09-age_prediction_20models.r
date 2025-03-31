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

manymodels <- list()

for (modeli in 1:20) {
    print(modeli / 20)
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
        tuneLength = 10
    )
    test_pred <- predict(rf, testing)
    performancemetrics <- confusionMatrix(test_pred, factor(testing$age))
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
    manymodels[[modeli]] <- list(
        rf = rf,
        test_pred = test_pred,
        performancemetrics = performancemetrics,
        conf_matrix = conf_matrix
    )
}
objsave(manymodels, "processed/age_pred/manymodels")

summary(sapply(
    readRDS("processed/age_pred/manymodels.rds"),
    function(x) x$performancemetrics$overall["Accuracy"]
))

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7292  0.8217  0.8627  0.8592  0.8974  0.9600



sort(table(c(unlist(sapply(
    readRDS("processed/age_pred/manymodels.rds"),
    function(x) {
        library(caret)
        rf <- x$rf
        mostimportant_5feat <- (varImp(rf)$importance %>%
            mutate(ASV = rownames(.)) %>%
            arrange(-Overall) %>%
            head(5))$ASV

        mostimportant_5feat_names <- unname(apply(tax_table(ps)[mostimportant_5feat, ], 1, function(x) {
            paste(x[!is.na(x)], collapse = ";")
        }))
    }
)))), dec = T)

#   Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;ascendens
#                                                                                                   20
# Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;pasteurianus
#                                                                                                   20
#     Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;persici
#                                                                                                   20
#  Bacteria;Proteobacteria;Alphaproteobacteria;Acetobacterales;Acetobacteraceae;Acetobacter;tropicalis
#                                                                                                   20
#     Bacteria;Proteobacteria;Gammaproteobacteria;Burkholderiales;Burkholderiaceae;Ralstonia;pickettii
#                                                                                                   20

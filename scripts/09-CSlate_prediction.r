source("scripts/00-setup.r")
library(caret)

ps <- readRDS("processed/ps_pruned_rared10k.rds")
ps <- tax_glom(ps, taxrank = "Species")
ps <- microbiome::transform(ps, "clr")
metadata <- data.frame(sample_data(ps))
df <- as.data.frame((otu_table(ps)))
df$CSlate <- metadata[rownames(df), "CSlate"]
CSlate <- metadata %>%
    select(CSlate, Isoline) %>%
    unique()
CSlate <- setNames(CSlate$CSlate, CSlate$Isoline)
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
    rf <- train(CSlate ~ .,
        data = training,
        method = "rf",
        trControl = fit_ctrl,
        verbose = FALSE,
        tuneLength = 10
    )
    test_pred <- predict(rf, testing)
    rsq <- summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared
    mae <- mean(abs(predict(rf, testing) - testing$CSlate))
    noinfmae <- mean(abs(testing$CSlate - median(training$CSlate)))
    conf_matrix <- data.frame(y = test_pred, x = testing$CSlate, type = "test") %>%
        rbind(data.frame(y = predict(rf, training), x = training$CSlate, type = "train")) %>%
        ggplot(aes(x = x, y = y, color = type)) +
        geom_abline() +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 2.5) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        ggtitle(paste(
            "MAE:", round(mean(abs(predict(rf, testing) - testing$CSlate)), 2), "R2:",
            round(summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared, 2)
        ))
    conf_matrix
    manymodels[[modeli]] <- list(
        rf = rf,
        test_pred = test_pred,
        rsq = rsq,
        mae = mae,
        noinfmae = noinfmae,
        conf_matrix = conf_matrix
    )
}
objsave(manymodels, "processed/CSlate_pred/manymodels")

rm(list = ls())

source("scripts/00-setup.r")
library(caret)

ps <- readRDS("processed/ps_pruned_rared10k.rds")
ps <- prune_samples(sample_data(ps)$age == "late", ps)
ps <- tax_glom(ps, taxrank = "Species")
ps <- microbiome::transform(ps, "clr")
metadata <- data.frame(sample_data(ps))
df <- as.data.frame((otu_table(ps)))
df$CSlate <- metadata[rownames(df), "CSlate"]
CSlate <- metadata %>%
    select(CSlate, Isoline) %>%
    unique()
CSlate <- setNames(CSlate$CSlate, CSlate$Isoline)
isolines <- setNames(metadata$Isoline, metadata$sampleID)

manymodels_late <- list()

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
    rf <- train(CSlate ~ .,
        data = training,
        method = "rf",
        trControl = fit_ctrl,
        verbose = FALSE,
        tuneLength = 10
    )
    test_pred <- predict(rf, testing)
    rsq <- summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared
    mae <- mean(abs(predict(rf, testing) - testing$CSlate))
    noinfmae <- mean(abs(testing$CSlate - median(training$CSlate)))
    conf_matrix <- data.frame(y = test_pred, x = testing$CSlate, type = "test") %>%
        rbind(data.frame(y = predict(rf, training), x = training$CSlate, type = "train")) %>%
        ggplot(aes(x = x, y = y, color = type)) +
        geom_abline() +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 2.5) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        ggtitle(paste(
            "MAE:", round(mean(abs(predict(rf, testing) - testing$CSlate)), 2), "R2:",
            round(summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared, 2)
        ))
    conf_matrix
    manymodels_late[[modeli]] <- list(
        rf = rf,
        test_pred = test_pred,
        rsq = rsq,
        mae = mae,
        noinfmae = noinfmae,
        conf_matrix = conf_matrix
    )
}
objsave(manymodels_late, "processed/CSlate_pred/manymodels_lateonly")

rm(list = ls())

source("scripts/00-setup.r")
library(caret)

ps <- readRDS("processed/ps_pruned_rared10k.rds")
ps <- prune_samples(sample_data(ps)$age == "early", ps)
ps <- tax_glom(ps, taxrank = "Species")
ps <- microbiome::transform(ps, "clr")
metadata <- data.frame(sample_data(ps))
df <- as.data.frame((otu_table(ps)))
df$CSlate <- metadata[rownames(df), "CSlate"]
CSlate <- metadata %>%
    select(CSlate, Isoline) %>%
    unique()
CSlate <- setNames(CSlate$CSlate, CSlate$Isoline)
isolines <- setNames(metadata$Isoline, metadata$sampleID)

manymodels_early <- list()

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
    rf <- train(CSlate ~ .,
        data = training,
        method = "rf",
        trControl = fit_ctrl,
        verbose = FALSE,
        tuneLength = 10
    )
    test_pred <- predict(rf, testing)
    rsq <- summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared
    mae <- mean(abs(predict(rf, testing) - testing$CSlate))
    noinfmae <- mean(abs(testing$CSlate - median(training$CSlate)))
    conf_matrix <- data.frame(y = test_pred, x = testing$CSlate, type = "test") %>%
        rbind(data.frame(y = predict(rf, training), x = training$CSlate, type = "train")) %>%
        ggplot(aes(x = x, y = y, color = type)) +
        geom_abline() +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 2.5) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        ggtitle(paste(
            "MAE:", round(mean(abs(predict(rf, testing) - testing$CSlate)), 2), "R2:",
            round(summary(lm(predict(rf, testing) ~ testing$CSlate))$r.squared, 2)
        ))
    conf_matrix
    manymodels_early[[modeli]] <- list(
        rf = rf,
        test_pred = test_pred,
        rsq = rsq,
        mae = mae,
        noinfmae = noinfmae,
        conf_matrix = conf_matrix
    )
}
objsave(manymodels_early, "processed/CSlate_pred/manymodels_earlyonly")

summary(sapply(readRDS("processed/CSlate_pred/manymodels.rds"), function(x) x$rsq))
summary(sapply(readRDS("processed/CSlate_pred/manymodels_earlyonly.rds"), function(x) x$rsq))
summary(sapply(readRDS("processed/CSlate_pred/manymodels_lateonly.rds"), function(x) x$rsq))
rm(list = ls())

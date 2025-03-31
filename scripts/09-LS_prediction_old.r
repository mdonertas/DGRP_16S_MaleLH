source("scripts/00-setup.r")

library(caret)

get_best_result <- function(caret_fit) {
    best <- which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
    best_result <- caret_fit$results[best, ]
    rownames(best_result) <- NULL
    best_result
}

ps <- readRDS("processed/ps_pruned_rared10k.rds")
ps <- prune_samples(sample_data(ps)$age == "late", ps)
ps <- tax_glom(ps, taxrank = "Species")

otu_table(ps) <- otu_table((otu_table(ps) > quantile(c(otu_table(ps)), 0.50)) + 1 - 1)

# ps <- microbiome::transform(ps, "clr")
metadata <- data.frame(sample_data(ps))


df <- as.data.frame((otu_table(ps)))

df$Lifespan <- metadata[rownames(df), "AvLF"]

lifespan <- metadata %>%
    select(AvLF, Isoline) %>%
    unique()

lifespan <- setNames(lifespan$AvLF, lifespan$Isoline)
isolines <- setNames(metadata$Isoline, metadata$sampleID)


manymodels <- list()

for (modeli in 1:20) {
    print(modeli / 20)

    train_isolines <- names(lifespan[createDataPartition(lifespan,
        p = 0.8,
        list = FALSE
    )[, 1]])
    test_isolines <- setdiff(unique(isolines), train_isolines)

    train_samples <- names(isolines[isolines %in% train_isolines])
    test_samples <- names(isolines[isolines %in% test_isolines])

    training <- df[train_samples, ]
    testing <- df[test_samples, ]
    train_lifespan <- lifespan[train_isolines]
    fitset <- sapply(1:20, function(i) {
        which(rownames(training) %in% names(isolines[isolines %in% names(train_lifespan[createDataPartition(train_lifespan, p = 15 / 22, list = FALSE)[, 1]])]))
    })

    fit_ctrl <- trainControl(
        method = "repeatedcv",
        # index = fitset,
        number = 5,
        repeats = 5
    )
    rf <- train(Lifespan ~ .,
        data = training,
        method = "rf",
        trControl = fit_ctrl,
        verbose = FALSE,
        tuneLength = 10,
        metric = "MAE"
    )

    get_best_result(rf)
    rsq <- summary(lm(predict(rf, testing) ~ testing$Lifespan))$r.squared
    mae <- mean(abs(predict(rf, testing) - testing$Lifespan))
    noinfmae <- mean(abs(testing$Lifespan - median(training$Lifespan)))
    data.frame(y = predict(rf, testing), x = testing$Lifespan, type = "test") %>%
        rbind(data.frame(y = predict(rf, training), x = training$Lifespan, type = "train")) %>%
        ggplot(aes(x = x, y = y, color = type)) +
        geom_abline() +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 2.5) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        ggtitle(paste(
            "MAE:", round(mean(abs(predict(rf, testing) - testing$Lifespan)), 2), "R2:",
            round(summary(lm(predict(rf, testing) ~ testing$Lifespan))$r.squared, 2)
        ))

    manymodels[[modeli]] <- list(
        mae, rsq, noinfmae,
        rf, training, testing
    )
}




hist(unlist(sapply(manymodels, function(x) x[[2]])))
hist(unlist(sapply(manymodels, function(x) x[[1]])))
hist(unlist(sapply(manymodels, function(x) x[[3]] - x[[1]])))

summary(sapply(manymodels, function(x) x[[2]]))
summary(sapply(manymodels, function(x) x[[1]]))

objsave(manymodels, "processed/ls_prediction/manymodels.rds")

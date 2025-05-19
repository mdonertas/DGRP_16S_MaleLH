source("scripts/00-setup.r")
library(caret)

rbind(
    as.data.frame(t(sapply(readRDS("processed/CSearly_pred/manymodels.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Climbing Speed") %>%
        mutate(type = "all samples"),
    as.data.frame(t(sapply(readRDS("processed/CSearly_pred/manymodels_earlyonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Climbing Speed") %>%
        mutate(type = "early samples only"),
    as.data.frame(t(sapply(readRDS("processed/CSearly_pred/manymodels_lateonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Climbing Speed") %>%
        mutate(type = "late samples only"),
    as.data.frame(t(sapply(readRDS("processed/CSSlope_pred/manymodels.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Functional Aging") %>%
        mutate(type = "all samples"),
    as.data.frame(t(sapply(readRDS("processed/CSSlope_pred/manymodels_earlyonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Functional Aging") %>%
        mutate(type = "early samples only"),
    as.data.frame(t(sapply(readRDS("processed/CSSlope_pred/manymodels_lateonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Functional Aging") %>%
        mutate(type = "late samples only"),
    as.data.frame(t(sapply(readRDS("processed/EarlyRS_pred/manymodels.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Reproductive Success") %>%
        mutate(type = "all samples"),
    as.data.frame(t(sapply(readRDS("processed/EarlyRS_pred/manymodels_earlyonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Reproductive Success") %>%
        mutate(type = "early samples only"),
    as.data.frame(t(sapply(readRDS("processed/EarlyRS_pred/manymodels_lateonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Early Reproductive Success") %>%
        mutate(type = "late samples only"),
    as.data.frame(t(sapply(readRDS("processed/AvLF_pred/manymodels.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Average Lifespan") %>%
        mutate(type = "all samples"),
    as.data.frame(t(sapply(readRDS("processed/AvLF_pred/manymodels_earlyonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Average Lifespan") %>%
        mutate(type = "early samples only"),
    as.data.frame(t(sapply(readRDS("processed/AvLF_pred/manymodels_lateonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Average Lifespan") %>%
        mutate(type = "late samples only"),
    as.data.frame(t(sapply(readRDS("processed/Rsen_pred/manymodels.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Reproductive Senescence") %>%
        mutate(type = "all samples"),
    as.data.frame(t(sapply(readRDS("processed/Rsen_pred/manymodels_earlyonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Reproductive Senescence") %>%
        mutate(type = "early samples only"),
    as.data.frame(t(sapply(readRDS("processed/Rsen_pred/manymodels_lateonly.rds"), function(x) c(x$rsq, x$mae, x$noinfmae)))) %>%
        set_names(c("R2", "MAE", "MAE_noinf")) %>%
        mutate(LH = "Reproductive Senescence") %>%
        mutate(type = "late samples only")
) %>%
    tablesave("results/LH_prediction_allresults")

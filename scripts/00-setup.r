# compilation of common libraries, functions, visualisation settings across
# scripts. each script starts by running this setup script.

#### Folder Organisation ####
# sapply(
#     c("docs", "processed", "raw", "scripts", "results", "notebooks"),
#     function(folder) {
#         system(paste("mkdir -p", folder))
#     }
# ) # run only at the beginning of the project

#### initialize renv ####
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}

# check if renv folder exists, if not initialize renv
if (!file.exists("renv")) {
    renv::init()
}

#### Libraries ####
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(Rvislib)

#### Visualisation Settings ####
theme_set(theme_pubr(base_size = 8))
pntnorm <- (1 / 0.352777778)
geom_point2 <- function(...) ggplot2::geom_point(shape = 21, ...)

## Colors ##
agecol <- setNames(c("#2D78AB", "#A83F3F"), c("early", "late"))

#### Functions ####

plotsave <- function(ggobj, prefix, width = 8, height = 8, ...) {
    path <- strsplit(prefix, "/")[[1]]
    path <- paste(path[-length(path)], collapse = "/")
    if (!file.exists(path)) {
        system(paste("mkdir -p", path))
    }
    saveRDS(object = ggobj, file = paste(prefix, ".rds", sep = ""))
    ggplot2::ggsave(
        file = paste(prefix, ".pdf", sep = ""), plot = ggobj, units = "cm",
        width = width, height = height, useDingbats = FALSE, limitsize = FALSE
    )
    ggplot2::ggsave(
        file = paste(prefix, ".png", sep = ""), plot = ggobj, units = "cm",
        width = width, height = height, limitsize = FALSE
    )
}

tablesave <- function(tib, prefix, ...) {
    path <- strsplit(prefix, "/")[[1]]
    path <- paste(path[-length(path)], collapse = "/")
    if (!file.exists(path)) {
        system(paste("mkdir -p", path))
    }
    readr::write_csv(tib,
        file = paste(prefix, ".csv", sep = ""),
        append = FALSE
    )
    readr::write_tsv(tib,
        file = paste(prefix, ".tsv", sep = ""),
        append = FALSE
    )
    saveRDS(object = tib, file = paste(prefix, ".rds", sep = ""))
}

objsave <- function(obj, prefix) {
    path <- strsplit(prefix, "/")[[1]]
    path <- paste(path[-length(path)], collapse = "/")
    if (!file.exists(path)) {
        system(paste("mkdir -p", path))
    }
    saveRDS(object = obj, file = paste(prefix, ".rds", sep = ""))
}

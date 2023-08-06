source("scripts/00-setup.R")

sample_list <- strsplit(readLines("./raw/sampleList"), ":")
sample_list <- setNames(
    strsplit(sapply(sample_list, function(x) x[2]), " "),
    sapply(sample_list, function(x) x[1])
)

sample_list <- reshape2::melt(sample_list) %>%
    set_names(c("sampleID", "run"))

dnaext <- readxl::read_xlsx("./raw/DNA extraction batches.xlsx") %>%
    set_names(c("sampleID", "DNAbatch"))

sample_list <- sample_list %>%
    inner_join(dnaext)

sample_list$age <- substr(sample_list$sampleID, 1, 1)
sample_list$Isoline <- substr(
    sample_list$sampleID, 2,
    nchar(sample_list$sampleID) - 1
)

lifehist <- read_tsv("./raw/LifeHistoryVar.txt") %>%
    mutate(Isoline = as.character(Isoline))

sample_list %>%
    str()

# 'data.frame':   177 obs. of  5 variables:
#  $ sampleID: chr  "E10A" "E10B" "E10C" "E11A" ...
#  $ run     : chr  "RUNGEN52_2019" "RUNGEN52_2019" "RUNGEN52_2019"
#  $ DNAbatch: chr  "group9" "group6" "group7" "group9" ...
#  $ age     : chr  "E" "E" "E" "E" ...
#  $ Isoline : chr  "10" "10" "10" "11" ...

apply(sample_list, 2, function(x) unique(x))

# $sampleID
#   [1] "E10A"  "E10B"  "E10C"  "E11A"  "E11B"  "E11C"  "E12A"  "E12B"  "E12C"
# "E14A"  "E14B"  "E14D"  "E15A"  "E15C"  "E15D"  "E17A"
#  [17] "E17B"  "E17C"  "E17D"  "E19A"  "E19B"  "E19C"  "E19D"  "E20A"  "E20B"
# "E20C"  "E22A"  "E22B"  "E22D"  "E23B"  "E23C"  "E23D"
#  [33] "E24A"  "E24B"  "E24C"  "E6A"   "E6B"   "E6D"   "L10A"  "L10B"  "L10C"
# "L10D"  "L11A"  "L11B"  "L11C"  "L11D"  "L12A"  "L12B"
#  [49] "L12C"  "L12D"  "L14A"  "L14B"  "L14C"  "L14D"  "L15A"  "L15B"  "L15C"
# "L17A"  "L17B"  "L17C"  "L17D"  "L19A"  "L19B"  "L19C"
#  [65] "L19D"  "L20A"  "L20B"  "L20CD" "L22A"  "L22B"  "L22C"  "L22D"  "L23A"
# "L23D"  "L24A"  "L24B"  "L6A"   "L6B"   "L6C"   "L6D"
#  [81] "E22C"  "E23A"  "E25A"  "E25B"  "E25C"  "E25D"  "E26A"  "E26B"  "E26C"
# "E26D"  "E27A"  "E27B"  "E27C"  "E28A"  "E28B"  "E28C"
#  [97] "E28D"  "E29A"  "E29B"  "E29C"  "E29D"  "E30A"  "E30B"  "E30C"  "E30D"
# "E31A"  "E31B"  "E31C"  "E31D"  "E33A"  "E33B"  "E33C"
# [113] "E33D"  "E35A"  "E35B"  "E35C"  "E35D"  "E36A"  "E36B"  "E36C"  "E36D"
# "E38A"  "E38B"  "E38C"  "E38D"  "E39A"  "E39B"  "E39C"
# [129] "E39D"  "L23B"  "L24C"  "L25A"  "L25B"  "L25C"  "L25D"  "L26A"  "L26B"
# "L26C"  "L26D"  "L27A"  "L27B"  "L27C"  "L28A"  "L28B"
# [145] "L28C"  "L28D"  "L29A"  "L29B"  "L29C"  "L29D"  "L30A"  "L30B"  "L30C"
# "L30D"  "L31A"  "L31B"  "L31C"  "L31D"  "L33A"  "L33B"
# [161] "L33C"  "L33D"  "L35B"  "L35C"  "L35D"  "L36A"  "L36B"  "L36C"  "L36D"
# "L38A"  "L38B"  "L38C"  "L38D"  "L39A"  "L39B"  "L39C"
# [177] "L39D"

# $run
# [1] "RUNGEN52_2019" "RUNGEN66"

# $DNAbatch
#  [1] "group9"  "group6"  "group7"  "group5"  "group10" "group2"  "group1"
# "group8"  "group4"  "group11" "group3"  "group13" "group12"
# [14] "NA"      "group15" "group18" "group14" "group16" "group17"

# $age
# [1] "E" "L"

# $Isoline
#  [1] "10"  "11"  "12"  "14"  "15"  "17"  "19"  "20"  "22"  "23"  "24"  "6"
# "20C" "25"  "26"  "27"  "28"  "29"  "30"  "31"  "33"
# [22] "35"  "36"  "38"  "39"

# there is an NA in DNAbatch, it will throw a warning when we run the following
# code as it isn't a number after removing group pre-fix
sample_list <- sample_list %>%
    left_join(lifehist) %>%
    mutate(DNAbatch = as.factor(as.numeric(gsub("group", "", DNAbatch)))) %>%
    mutate(run = as.factor(as.numeric(as.factor(run)))) %>%
    mutate(age = c("early", "late")[1 + (age == "L")])

objsave(sample_list, "./processed/metadata")

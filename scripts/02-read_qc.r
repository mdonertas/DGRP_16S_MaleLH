source("scripts/00-setup.R")
library(fastqcr)

run1_path <- "./raw/RUNGEN52_2019/fastq"
run2_path <- "./raw/RUNGEN66/fastq"

# run fastqc
fastqc(fq.dir = run1_path, qc.dir = "./processed/fastqc/run1")
fastqc(fq.dir = run2_path, qc.dir = "./processed/fastqc/run2")

# aggregate fastqc reports
qc_run1 <- qc_aggregate("./processed/fastqc/run1", progressbar = FALSE)
qc_run2 <- qc_aggregate("./processed/fastqc/run2",
    progressbar = FALSE,
    show_col_types = FALSE
)

# select only the samples in our study and save their QC report
metadata <- readRDS("./processed/metadata.rds")
qc_run1 <- qc_run1 %>%
    separate(sample, into = c("SampleID", "FR"), remove = FALSE) %>%
    filter(SampleID %in% metadata$sampleID)
qc_run2 <- qc_run2 %>%
    separate(sample, into = c("SampleID", "FR"), remove = FALSE) %>%
    filter(SampleID %in% metadata$sampleID)

objsave(qc_run1, "./processed/fastqc/run1_summary")
objsave(qc_run2, "./processed/fastqc/run2_summary")

qc_run1 %>%
    filter(status == "FAIL") %>%
    group_by(module, FR) %>%
    summarise(n = length(unique(SampleID)))

qc_run2 %>%
    filter(status == "FAIL") %>%
    group_by(module, FR) %>%
    summarise(n = length(unique(SampleID)))

qc_run1 %>%
    filter(module == "Adapter Content") %>%
    group_by(FR) %>%
    summarise(percent_pass = mean(status == "PASS") * 100)

# No adaptors detected. Reads are trimmed for primers/adaptors. There is a
# decline in quality towards the end of reads as expected, we will truncate the
# reads where significant drop starts while running the DADA2 pipeline.

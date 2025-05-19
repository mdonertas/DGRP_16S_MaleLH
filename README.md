# README

This readme file was generated on [2023-06-08].

## GENERAL INFORMATION

Title of Project: Age-dependent gut microbiota dynamics and their association with male fitness traits in Drosophila melanogaster

Project Abstract:
Growing evidence suggests that the gut microbiota plays a key role in shaping life history in a wide range of species, including well-studied model organisms like Drosophila melanogaster. Although recent studies have explored the relationship between gut microbiota and female life history, the link between gut microbiota and male life history remains understudied. In this study, we explored the role of gut microbiota in shaping male life history traits by correlating variation in life history traits across genetically homogeneous isolines with their naturally occurring gut microbiota. Using 22 isolines from the Drosophila melanogaster Genetic Reference Panel (DGRP), we measured lifespan, early/late-life reproduction, and early/late-life physiological performance. We characterized the gut microbiota composition in young (5 days old) and old (26 days old) flies using 16S rDNA sequencing. We observed significant variation in male life history traits across isolines, as well as age-related changes in gut microbiota composition. Using machine learning, we showed that gut microbiota composition could predict the age of the organisms with high accuracy. Associations between gut microbiota and life history traits were notable, particularly involving the Acetobacter genus. In early life, the abundance of Acetobacter ascendens was associated with functional aging, while Acetobacter indonesiensis was linked to reproductive senescence. In late life, higher abundances of A. ascendens and Acetobacter pasteurianus were negatively associated with lifespan. These findings highlight the potential role of gut microbiota, especially the Acetobacter genus, in male fitness and aging. 

[Github link](https://github.com/mdonertas/DGRP_16S_MaleLH)

## DATA

### RAW DATA

Location of the raw data used for this study:

- Local: ./raw/
- Publicly available data repository: ENA PRJEB88786 

### PROCESSED DATA

Location of the processed data:

- Local: ./processed/
- Publicly available data repository: BioStudies S-BSST2029

## OTHER DOCUMENTATION

## SOFTWARE

In this project, from beginning to end, `renv` R package was used. All required packages, including the versions are listed in the `renv.lock` file. To install all the packages, run the following command in R:

```r
renv::restore()
```

All analysis was done within R programming environment. The only additional tool used was FastQC (v 0.12.1) for quality control of raw reads. All scripts are available in the `scripts/` folder. The scripts are numbered in the order they were run.

## CONTACTS

- Name: Melike Dönertaş
  - Contribution[^1]: Data preprocessing, analysis, visualization, conceptualization of analysis, code writing, documenting
  - ORCID: [0000-0002-9788-6535](https://orcid.org/0000-0002-9788-6535)
  - Institution: Leibniz Institute on Aging - Fritz Lipmann Institute (FLI), Jena, Germany
  - Email: <melike.donertas@leibniz-fli.de>
- Name: Zahida Sultanova
  - Contribution[^1]: Conceptualization of analysis
  - ORCID:
  - Institution: School of Biological Sciences, University of East Anglia, Norwich, UK
  - Email: <Zahida.Sultanova@uea.ac.uk>

[^1]: Please note that this is not necessarily the same as the manuscript as the manuscript involves wetlab work and other people. Contributions listed here are only related to the analysis presented in this repo. For more information about the project, please refer to the paper.

## FUNDING

This study was supported by the Spanish Ministry of Economy and Competitivity grants: PGC2018-099344-B-I00 (AL), CGL2014-58722-P (PC and JILL), RYC-2013-12998 (PC), and RYC-2012-11872 (JILL). HMD is funded by Carl-Zeiss-Stiftung (P2021-00-007). PA was supported by FCT through a PhD grant (number 2021.05611). ZS was supported by the Atracció de Talent Fellowship and the Leverhulme Trust Early Career Fellowship (ECF-2022-214).
##### R Code for the updated systematic search of "Meta-analysis on reporting practices as a source of heterogeneity in in vitro cancer research" ######
##### by Timo Sander, Joly Ghanawi, Emma Wilson, Sajjad Muhammad, Malcolm MacLeod and Ulf Dietrich Kahlert #####
##### April 2022 #####


### Loading/preparing data exported from SyRF

# Setting working directory:
setwd("/users/tsand/desktop/SERRANO FULL TEXT Annotations/")

# Loading libraries used for further calculations:
library(dplyr)
library(plyr)
library(readr)
library(tidyverse)
library(scales)
library(metafor)
library(clubSandwich)
library(orchaRd)
library(meta)
library(stringr)

# Loading study annotations:
SA <- read.csv("SERRANO Update Annotations Study.csv")

# Loading model annotations:
MA <- read.csv("SERRANO Update Annotations Model.csv")



### Analysis of Full-Text Screening:

# How many studies in total:
SA %>% select(StudyIdStr) %>% unique() %>% nrow()


## Are there two annotations for each article:

# Adding a column with joined StudyID and Relevance:
SA$joined_StudyID_Annotator <- paste(SA$StudyIdStr, SA$AnnotatorIdStr)

# Checking if there are two annotations for each article:
SA %>% select(joined_StudyID_Annotator) %>% unique() %>% nrow()
# as there are 198 every article has 2 screenings (99 x 2)


# Unique screening decisions per article:
SA %>% distinct(StudyIdStr, X7b03292b.4886.4d92.aeaa.2a37b54c33f1_Is.this.article.relevant.for.our.review...Full.Text.Screening.) %>% select(StudyIdStr) %>% count()

# Articles screeners agreed on:
SA_Agree <- SA %>% distinct(StudyIdStr, X7b03292b.4886.4d92.aeaa.2a37b54c33f1_Is.this.article.relevant.for.our.review...Full.Text.Screening.) %>% select(StudyIdStr) %>% count() %>% filter(freq < 2) %>% select(StudyIdStr)

# Articles screeners disagreed on:
SA_Disagree <- SA %>% distinct(StudyIdStr, X7b03292b.4886.4d92.aeaa.2a37b54c33f1_Is.this.article.relevant.for.our.review...Full.Text.Screening.) %>% select(StudyIdStr) %>% count() %>% filter(freq > 1) %>% select(StudyIdStr)


# Filtering the screenings with disagreements:
SA %>% filter(is.element(StudyIdStr, SA_Disagree$StudyIdStr)) %>% write.csv("SERRANO Update Screening disagreements1.csv", row.names = F)

# Filtering the articles that are already included:
SA_already_in <-SA %>% filter(is.element(StudyIdStr, SA_Agree$StudyIdStr)) %>% filter(X7b03292b.4886.4d92.aeaa.2a37b54c33f1_Is.this.article.relevant.for.our.review...Full.Text.Screening. == TRUE) %>% select(StudyIdStr) %>% unique()

# Filtering the articles that are already excluded:
SA_already_ex <-SA %>% filter(is.element(StudyIdStr, SA_Agree$StudyIdStr)) %>% filter(X7b03292b.4886.4d92.aeaa.2a37b54c33f1_Is.this.article.relevant.for.our.review...Full.Text.Screening. == FALSE) %>% select(StudyIdStr) %>% unique()


## Screening reconciliation by E. W.:

# Emma excluded every article except the one with the ID "4a5de13a-a78e-4f61-ad28-59d4e30680c0".

# Adding a column to SA_Disagree with Emmas screening decision:
SA_Disagree$emmas_screening <- "ex"
SA_Disagree[SA_Disagree$StudyIdStr == "4a5de13a-a78e-4f61-ad28-59d4e30680c0","emmas_screening"] <- "in"



### Analysis of reasons for exclusion

# Final list of excluded articles:
SA_ex <- as.data.frame(c(SA_already_ex$StudyIdStr, SA_Disagree$StudyIdStr[SA_Disagree$emmas_screening == "ex"]))

# Getting reasons for exclusion:
SA %>% filter(is.element(StudyIdStr, SA_ex$`c(SA_already_ex$StudyIdStr, SA_Disagree$StudyIdStr[SA_Disagree$emmas_screening == "ex"])`)) %>% write.csv("Serrano Update Reasons for exclusion1.csv", row.names = F)




### Analysing extracted data

## Preparations:

# Filtering the articles that are already included:
MA %>% filter(is.element(StudyIdStr, SA_already_in$StudyIdStr))

# Exporting:
MA %>% filter(is.element(StudyIdStr, SA_already_in$StudyIdStr)) %>% write.csv("MA_already_in1.csv", row.names = F)


## Importing the reconciled annotations:
MA <- read.csv("SERRANO Update Annotations FINAL incl reconciliation.csv")


## Getting one row per article:

# Filling in Emmas rows where she did not have to annotate:
for (a in seq(3,123,3)) {
  MA[a, seq(4,56,2)] <- ifelse(is.na(MA[a, seq(4,56,2)]) | MA[a, seq(4,56,2)] == "", MA[a-1, seq(4,56,2)], MA[a, seq(4,56,2)])
}

# Replacing "1" with TRUE
for (a in seq(3,123,3)) {
  MA[a, seq(4,56,2)] <- ifelse(MA[a, seq(4,56,2)] == 1, TRUE, MA[a, seq(4,56,2)])
}

# Getting StudyIDs in the reconciled rows:
for (a in seq(3,123,3)) {
  MA[a,2] <- MA[a-1,2]
}

# Selecting only the reconciled rows:
MA <- MA %>% filter(AnnotatorIdStr == "Emma")


## Reporting frequencies:

# Glucose
table(MA$X084199b6.59b8.4552.b67a.e9916fa36798_Was.the.glucose.level.of.the.DMEM.culture.medium.reported..high.low.no.glucose..)

# Sample Size Calculation
table(MA$X10a27061.dba1.4890.bd76.366c4e7c892f_Was.a.sample.size.calculation.for.the.number.of.cells.exp..replicates.reported.)

# Conflicts of interests statement
table(MA$X1647bbf1.9e29.4688.9cba.b39d9c5659f8_Is.a.statement.about.conflicts.of.interests.available.)

# FBS use (yes)
table(MA$X1aae7fe5.c90a.46a3.82f5.e9a95e92c589_Was.fetal.bovine..calf..serum..FBS.FCS..used.in.cell.culture.)

# Volume added TMZ
table(MA$X1edb75dd.60a8.48e4.a01c.bcb2e45e1538_Was.the.volume.of.added.Temozolomide.treatment.reported.)

# Antibiotics
table(MA$X335e2583.e959.4e86.80c2.3348619318ef_Was.reported.whether.antibiotics.were.supplemented.to.the.cell.culture.medium.)

# Volume added control
table(MA$X36396819.de6d.4959.a803.57d1fd0e9359_Was.the.volume.of.added.untreated.control.suspension.reported.)

# Protocol
table(MA$X3f678a88.9dc8.47ed.a2f6.f8c15ac56a54_Is.an.open.access.study.protocol.available.)

# Clear number ind. Exp. And repl.
table(MA$X460000cd.a8a2.428c.a2fe.01ef9ade6a60_Is.the.number.of.independent.experiments.and.replications.per.experiment.clear.)

# Volume per well
table(MA$X4f7202d0.27f7.48d6.af51.4f177ca678b5_Was.the.volume.per.well.reported.)

# Random allocation of cells to treatment and control group
table(MA$X535efa6b.de3a.46b0.9916.2d179b39098d_Was.reported.if.they.randomly.assigned.the.cells.to.control.and.intervention.)

# Blinded outcome assessment
table(MA$e3a76d9e.2e64.4e87.9a2c.29bcb8879970_Was.the.cell.viability.proliferation.measured.blinded.)

# U87 source
table(MA$X83f65c80.aeba.406e.b115.b410b3bae878_Was.the.source.of.U.87.MG.cells.reported.)

# U87 age
table(MA$X979d983b.ad64.4bde.b417.d03536ff3b02_Was.the.age.of.U.87.MG.cells.used.for.experiments.reported.)

# Errortype
table(MA$X98360455.8dfe.4cd6.899c.58ba8755c21f_Was.the.errortype..SD.SEM..in.results.of.U.87.MG.and.Temozolomide.reported.)

# Type of untreated control
table(MA$ac053e6e.9e19.4c37.8ab4.6ce95ee6aa1b_Was.type.of.untreated.control.in.Temozolomide.U.87.MG.viability.exp..reported.)

# Mycoplasma exclusion
table(MA$b257b57a.dcb9.4c26.8fe4.92e876a6386d_Was.reported.whether.a.mycoplasma.contamination.check.of.cell.culture.was.done.)

# U87 authentication
table(MA$b412c32b.ba2d.42b4.814f.3ca7f9e01f0e_Was.reported.whether.an.authentication.of.U.87.MG.cells.was.performed.)

# TMZ source
table(MA$b5b66773.1f3d.4810.829b.5406a710bc67_Was.the.source.of.Temozolomide..TMZ..reported..manufacturer..)

# FBS source
table(MA$c0e325d2.ce11.4356.9afc.52b4c2dbbbbd_Was.the.source.of.FBS.FCS.reported.)

# Cells per well
table(MA$c8ea8e2c.9fc9.4707.af00.455f7d4875f7_Was.the.number.of.cells.per.well.reported.)

# TMZ conc.
table(MA$ccbcebea.2846.4acc.95a3.ed440ad02714_Was.the.concentration.of.Temozolomide.in.U.87.MG.viability.experiments.reported.)

# Treatment duration
table(MA$da52bc21.d2b4.41be.aeeb.8bc79fb755ee_Was.the.treatment.duration.of.U.87.MG.cells.with.Temozolomide.reported.)

# Number of experiments
table(MA$dcbe3947.894c.4f8f.9d2d.5d5e7d691ba1_Was.the.number.of.experiments.replications.of.U.87.MG.and.Temozolomide.reported.)

# Clear way of calculation of mean and error values
table(MA$e302b000.3756.4206.ba0c.dc95be33031d_Is.the.way.of.calculation.for.cell.viability.average.and.error.values.clear.)

# Criteria for cell passaging
table(MA$e7fd642a.891c.4af0.b4cf.60ca4705b430_Where.criteria.for.cell.passaging.reported...e.g..time..number.of.passage.)

# Cell concentration
table(MA$f43c07d5.ff37.4e37.aeba.509b172a7acd_Was.the.concentration.of.U.87.MG.cells.reported.)

# CAVE: calculable cell concentration (if cells per well and volume per well is reported):
MA %>% filter(f43c07d5.ff37.4e37.aeba.509b172a7acd_Was.the.concentration.of.U.87.MG.cells.reported. == FALSE 
              & c8ea8e2c.9fc9.4707.af00.455f7d4875f7_Was.the.number.of.cells.per.well.reported. == TRUE 
              & X4f7202d0.27f7.48d6.af51.4f177ca678b5_Was.the.volume.per.well.reported. == TRUE) 
# For this six articles, it is possible to calculate cell concentration. Therefore, cell concentration is considered as reported for those six additional articles.


## Calculating articles reporting score for each article:

# Converting some columns into logical values:
for (b in 1:nrow(MA)) {
  MA[b,12] <- ifelse(MA[b,12] == "TRUE",1,0)
}
for (b in 1:nrow(MA)) {
  MA[b,16] <- ifelse(MA[b,16] == "TRUE",1,0)
}
for (b in 1:nrow(MA)) {
  MA[b,32] <- ifelse(MA[b,32] == "TRUE",1,0)
}
for (b in 1:nrow(MA)) {
  MA[b,50] <- ifelse(MA[b,50] == "TRUE",1,0)
}
for (b in 1:nrow(MA)) {
  MA[b,56] <- ifelse(MA[b,56] == "TRUE",1,0)
}

# Converting the columns into numerical format:
MA[,12] <- as.numeric(MA[,12])
MA[,16] <- as.numeric(MA[,16])
MA[,32] <- as.numeric(MA[,32])
MA[,50] <- as.numeric(MA[,50])
MA[,56] <- as.numeric(MA[,56])

# Converting the remaining columns into numerical values:
MA[,c(4,8,14,26,28,30,34,36,38,40,44,46,48,54)] <- ifelse(MA[,c(4,8,14,26,28,30,34,36,38,40,44,46,48,54)] == TRUE, 1, 0) 
# 1 means "reported", 0 means "not reported reported"

# Adding a column with number of reported study parameters:
MA$reported_study_parameters <- rowSums(MA[,c(4,8,12,16,14,26,28,30,32,34,36,38,40,44,46,48,54,56)])

# Adding a column with proportion of reported study parameters:
MA$proportion_reported_study_parameters <- MA$reported_study_parameters / 18


## Calculating prevalence of risks of bias for each article:

# Adding the columns for missing data of experiments and within experiments:
MA_already_in <- read.csv("MA_already_in_updated.csv")
all(MA$StudyIdStr == MA_already_in$StudyIdStr)
MA$data_for_every_TMZ_conc_duration <- MA_already_in$X59133e30.17f5.4610.9025.71c888faff1c_Where.data.for.every.used.drug.conc..and.time.within.an.experiment.reported.
MA$data_for_every_experiment <- MA_already_in$daa4a231.a325.4a4c.aeb8.2c9e85e5f05c_Were.data.presented.for.every.U.87.MG.and.Temozolomide.viability.experiment.

# Converting the remaining columns into numerical values:
MA[,c(6,18,20,24,52,60,61)] <- ifelse(MA[,c(6,18,20,24,52,60,61)] == "TRUE",0,1)
# 1 means "risk of bias prevalent", 0 means "risk of bias not prevalent"

# Adding a column with number of prevalent risks of bias:
MA$prevalent_risks_of_bias <- rowSums(MA[,c(6,18,20,24,50,52,60,61)])

# Adding a column with proportion of prevalent risks of bias:
MA$proportion_prevalent_risks_of_bias <- MA$prevalent_risks_of_bias / 8

# Exporting:
write.csv(MA, file = "MA_with_reporting_prevalences.csv", row.names = F)



### Reconciled clear numbers of independent experiments:

# Loading reconciled data:
Rec_IE <- read.csv("SERRANO Update Independent experiments clear reconciled.csv")

# Same StudyIDs?
all(Rec_IE$StudyID == MA$StudyIdStr)

# Updating reconciled data in MA:
MA$X460000cd.a8a2.428c.a2fe.01ef9ade6a60_Is.the.number.of.independent.experiments.and.replications.per.experiment.clear. <- Rec_IE$Reconciled

# Converting into numerical data:
MA$X460000cd.a8a2.428c.a2fe.01ef9ade6a60_Is.the.number.of.independent.experiments.and.replications.per.experiment.clear. <- ifelse(MA$X460000cd.a8a2.428c.a2fe.01ef9ade6a60_Is.the.number.of.independent.experiments.and.replications.per.experiment.clear. == "TRUE", 1, 0)




### Comparison of overall reporting:

# Loading articles reporting and risks of bias prevalence:
old <- read.csv("SERRANO reporting and rob until 2020.csv")
new <- read.csv("SERRANO reporting and rob after 2020.csv")

# Loading combined csv:
comb <- read.csv("SERRANO update comparison reporting and rob.csv")

# Changing reconciled clear numbers of independent experiments:
new[1,3] <- 0.5
new[2,3] <- 0.5
new[6,3] <- 0.625
new[28,3] <- 0.5
new[31,3] <- 0.625
comb[1,3] <- 0.5
comb[2,3] <- 0.5
comb[6,3] <- 0.625
comb[28,3] <- 0.5
comb[31,3] <- 0.625


## Mann-Whitney-U-Test (nichtparametrisch):

# Reporting quality
wilcox.test(comb$Reporting.Study.Parameters.after.2020 ~ comb$Search, conf.int = TRUE)

# Prevalence of risks of bias:
wilcox.test(comb$Prevalence.Risks.of.Bias.after.2020 ~ comb$Search)


## Mean and median Reporting:
mean(new$Reporting.Study.Parameters.after.2020)
median(new$Reporting.Study.Parameters.after.2020)
mean(old$Reporting.Study.Parameters.until.2020)
median(old$Reporting.Study.Parameters.until.2020)


## Mean and median Risks of bias prevalence:
mean(new$Prevalence.Risks.of.Bias.after.2020)
median(new$Prevalence.Risks.of.Bias.after.2020)
mean(old$Prevalence.Risks.of.Bias.until.2020)
median(old$Prevalence.Risks.of.Bias.until.2020)


## Plotting:

# Reporting quality:
pdf(file = "Figure Change in reporting quality of articles.pdf", width = 7, height = 7)
boxplot(comb$Reporting.Study.Parameters.after.2020*100 ~ comb$Search,
        names = c("Articles published since August 2020", "Systematic search in August 2020"),
        main = "Change in reporting quality of articles",
        xlab = "",
        ylab = "Ratio or reported parameters per article (%)",
        ylim = c(0,100),
        cex.lab = 1.1,
        cex.axis = 1.1)
dev.off()

# Prevalence of risks of bias:
pdf(file = "Figure Change in risks of bias prevalence.pdf", width = 7, height = 7)
boxplot(comb$Prevalence.Risks.of.Bias.after.2020*100 ~ comb$Search,
        names = c("Articles published since August 2020", "Systematic search in August 2020"),
        main = "Change in risks of bias prevalence",
        xlab = "",
        ylab = "Ratio or risks of bias per article (%)",
        ylim = c(0,100),
        cex.lab = 1.1,
        cex.axis = 1.1)
dev.off()
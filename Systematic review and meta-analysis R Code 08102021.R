##### R Code for "Meta-analysis on reporting practices as a source of heterogeneity in in vitro cancer research" ######
##### by Timo Sander, Joly Ghanawi, Emma Wilson, Sajjad Muhammad, Malcolm MacLeod and Ulf Dietrich Kahlert #####
##### October 2021 #####


# Setting working directory:
setwd("/users/tsand/desktop/SERRANO - Data for R Code/")

# Loading libraries used for further calculations:
library(dplyr)
library(plyr)
library(readr)
library(tidyverse)
library(scales)
library(metafor)
library(clubSandwich)
library(orchaRd)



#### Loading the raw data downloaded from the SyRF webapp after screening articles and extracting data

Araw <- read.csv("SERRANO_outcomes_final_with_additional_articles.csv")

# Switching the first and second column:
A <- Araw[,c(2,1)]

# Inserting a column that says which annotator has assessed this data:
A <- add_column(A, "Annotator" = ifelse(
  A$AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b", "Emma", ifelse(
    A$AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6", "Timo", "Joly")
), .after = 2)

# Adding more columns to th new dataframe based on the raw data in Araw:
A <- data.frame(A, Araw[,c(43:44,42,40,46,13:15,21:25,26:31,33:38,41,45,47,50:55,57:58,60:70,72:78)], check.names = FALSE)

# Filtering the reconciled annotations (the annotator Emma reconciled the extracted data from the first two reviewers Joly and Timo):
AR <- A %>% filter(A$Annotator == "Emma")

# Loading the list of included articles based on the screening results:
included_articles <- read.csv("Final_list_of_included_articles.csv")

# Removing excluded articles:
AR <- AR %>% filter(is.element(StudyIdStr, included_articles$StudyID))

# Are they reconciled data for every included article?
included_articles %>% filter(is.element(StudyID, AR$StudyIdStr)) %>% nrow()
# yes, 137 rows for 137 included articles!



#### Preparing the data format for further analyses


### Number formatting:

# Commas to points in decimal numbers:
AR$Average <- gsub(",",".",AR$Average)
AR$Time <- gsub(",",".",AR$Time)
AR$Error <- gsub(",",".",AR$Error)

# Formatting to numerical values:
AR$Time <- as.numeric(AR$Time)
AR$Average <- as.numeric(AR$Average)
AR$Error <- as.numeric(AR$Error)

# Correcting a typo:
AR$Time[AR$StudyIdStr == "672fe86b-5b54-48fe-b85e-20e6e9af8a37" & AR$Time == 4] <- 3


### Temozolomide (TMZ) concentrations:

# Creating a new dataframe with the focus on the outcome data:
ARD <- AR[1:nrow(AR),c(1,seq(4,11,1),15,28)]

# Adding a row to the ARD dataframe with the TMZ concentration presented in microM:
ARD$tmz_concentration_in_microM <- ifelse(ARD$CohortLabel == "control" |
                                            ARD$CohortLabel == "control 1" |
                                            ARD$CohortLabel == "control bz" |
                                            ARD$CohortLabel == "control cc" |
                                            ARD$CohortLabel == "control dmso" |
                                            ARD$CohortLabel == "control mts" |
                                            ARD$CohortLabel == "control mtt" |
                                            ARD$CohortLabel == "control tb" |
                                            ARD$CohortLabel == "control thymidine", 0, 
                                          ifelse(ARD$CohortLabel == "tmz 100 uM", 100, 
                                                 ifelse(ARD$CohortLabel == "tmz 200 uM", 200, 
                                                        ifelse(ARD$CohortLabel == "tmz 500 uM", 500, 
                                                               ifelse(ARD$CohortLabel == "tmz 50 uM", 50,
                                                                      ifelse(ARD$CohortLabel == "tmz 20 uM", 20,
                                                                             ifelse(ARD$CohortLabel == "tmz 10 uM", 10,
                                                                                    ifelse(ARD$CohortLabel == "tmz 0.5 uM", 0.5,
                                                                                           ifelse(ARD$CohortLabel == "tmz 5 uM", 5,
                                                                                                  ifelse(ARD$CohortLabel == "tmz 15 uM", 15,
                                                                                                         ifelse(ARD$CohortLabel == "tmz 2.5 uM", 2.5,
                                                                                                                ifelse(ARD$CohortLabel == "tmz 250 uM", 250,
                                                                                                                       ifelse(ARD$CohortLabel == "tmz 1000 uM", 1000,
                                                                                                                              ifelse(ARD$CohortLabel == "tmz 25 uM", 25,
                                                                                                                                     ifelse(ARD$CohortLabel == "tmz 4000 uM", 4000,
                                                                                                                                            ifelse(ARD$CohortLabel == "tmz 1 uM", 1,
                                                                                                                                                   ifelse(ARD$CohortLabel == "tmz 12.5 uM", 12.5,
                                                                                                                                                          ifelse(ARD$CohortLabel == "tmz 125 uM", 125,
                                                                                                                                                                 ifelse(ARD$CohortLabel == "tmz 150 uM", 150, 
                                                                                                                                                                        ifelse(ARD$CohortLabel == "tmz 2000 uM", 2000, 
                                                                                                                                                                               ifelse(ARD$CohortLabel == "tmz 300 uM", 300, 
                                                                                                                                                                                      ifelse(ARD$CohortLabel == "tmz 400 uM", 400, 
                                                                                                                                                                                             ifelse(ARD$CohortLabel == "tmz 40 uM", 40, 
                                                                                                                                                                                                    ifelse(ARD$CohortLabel == "tmz 60 uM", 60, 
                                                                                                                                                                                                           ifelse(ARD$CohortLabel == "tmz 70 uM", 70,
                                                                                                                                                                                                                  ifelse(ARD$CohortLabel == "tmz 75 uM", 75,
                                                                                                                                                                                                                         ifelse(ARD$CohortLabel == "tmz 7.5 uM", 7.5, 
                                                                                                                                                                                                                                ifelse(ARD$CohortLabel == "tmz 750 uM", 750,
                                                                                                                                                                                                                                       ifelse(ARD$CohortLabel == "tmz 800 uM", 800, ARD$CohortLabel))))))))))))))))))))))))))))
)

# Adding a column to the ARD dataframe where the TMZ concentrations get transformed to numerical data:
ARD$tmz_concentration_in_microM_numeric <- as.numeric(ARD$tmz_concentration_in_microM)

# Loading the finished manual tmz concentration values:
ARD_tmz_conc_manually_finished <- read.csv("ARD_tmz_conc_manually after reconciliation finished.csv")

# Ordering all by StudyID and CohortLabel:
ARD_tmz_conc_manually_finished <- ARD_tmz_conc_manually_finished[order(ARD_tmz_conc_manually_finished$StudyIdStr, ARD_tmz_conc_manually_finished$CohortLabel),]
ARD <- ARD[order(ARD$StudyIdStr, ARD$CohortLabel),]

# Checking if the StudyIDs and CohortLabels are identical:
ARD_tmz_conc_manually_finished$StudyIdStr == ARD$StudyIdStr[is.na(ARD$tmz_concentration_in_microM_numeric) == TRUE] & ARD_tmz_conc_manually_finished$CohortLabel == ARD$CohortLabel[is.na(ARD$tmz_concentration_in_microM_numeric) == TRUE]
# the three non-identical values are not relevant!

# Adding the manually assessed TMZ concentration values to the ARD dataframe:
ARD$tmz_concentration_in_microM[is.na(ARD$tmz_concentration_in_microM_numeric) == TRUE] <- ARD_tmz_conc_manually_finished$tmz_concentration_in_microM_numeric_manually

# Updating the numeric TMZ concentrations column:
ARD$tmz_concentration_in_microM_numeric <- as.numeric(ARD$tmz_concentration_in_microM)

# Checking if it was successful:
table(is.na(ARD$tmz_concentration_in_microM_numeric))
# Now there are only 15 columns left with the tmz concentration as NA (because in this studies the tmz concentration is unknown)

# Ordering AR by StudyID and CohortLabel:
AR <- AR[order(AR$StudyIdStr, AR$CohortLabel),]

# Checking if the StudyIDs and CohortLabels are identical:
all(AR$StudyIdStr == ARD$StudyIdStr & AR$CohortLabel == AR$CohortLabel)

# Adding the completed/corrected TMZ concentration values to the AR dataframe:
AR$tmz_concentration_in_microM_numeric <- ARD$tmz_concentration_in_microM_numeric

# Adding a new control column:
AR$Control <- ifelse(AR$tmz_concentration_in_microM_numeric == 0, TRUE, FALSE)

# Checking if the Intervention Control and the Control column of AR are identical:
table(AR$InterventionControl == AR$Control)
# There are several rows were they are not identical (probably because of mistakes during Syrf annotating)
# --> therefore we continue using the correct "Control" column in the following code




#### Updating the extracted data with the information obtained through contacting the authors of the included articles

# Creating a new dataframe with the updated parameters:
AR_updated <- AR

# Updating the parameters:
AR_updated$What.is.the.type.of.control.treatment.[AR_updated$StudyIdStr == "64a785ad-4714-4ec6-b731-4d66128602f7" | AR_updated$StudyIdStr == "4b98515f-057a-4d17-af4a-527e9b58ec36"] <- "other"
AR_updated$What.is.the.type.of.control.treatment.[AR_updated$StudyIdStr == "86d1eac3-4f54-4d81-bd3d-673f115de939" |
                                                    AR_updated$StudyIdStr == "e8c88cd2-830a-4b62-95e8-18fdeaa25347" |
                                                    AR_updated$StudyIdStr == "6f3a6613-897c-4bb6-9f5a-a510bf1ec6ad" |
                                                    AR_updated$StudyIdStr == "bb056e51-f833-4bbc-a8f1-eecfe62bf9c9"] <- "DMSO"
AR_updated$What.other.control.did.they.use.[AR_updated$StudyIdStr == "64a785ad-4714-4ec6-b731-4d66128602f7" | AR_updated$StudyIdStr == "4b98515f-057a-4d17-af4a-527e9b58ec36"] <- "medium only"
AR_updated$Was.a.cell.line.authentification.of.U87.cells.performed.[AR_updated$StudyIdStr == "64a785ad-4714-4ec6-b731-4d66128602f7" |
                                                                      AR_updated$StudyIdStr == "0bcb0299-c385-4142-93f6-6ee7d5767ee6" |
                                                                      AR_updated$StudyIdStr == "e8c88cd2-830a-4b62-95e8-18fdeaa25347"] <- "TRUE"
AR_updated$Was.a.cell.line.authentification.of.U87.cells.performed.[AR_updated$StudyIdStr == "86d1eac3-4f54-4d81-bd3d-673f115de939" |
                                                                      AR_updated$StudyIdStr == "9cb3f84c-86a8-42e7-9237-77684a4af6ec" |
                                                                      AR_updated$StudyIdStr == "399b4c31-88e5-4430-8324-090dc87864db" |
                                                                      AR_updated$StudyIdStr == "4b98515f-057a-4d17-af4a-527e9b58ec36" |
                                                                      AR_updated$StudyIdStr == "c2320d09-2c4f-413e-b3e3-a2f5a896a357" |
                                                                      AR_updated$StudyIdStr == "6f3a6613-897c-4bb6-9f5a-a510bf1ec6ad" |
                                                                      AR_updated$StudyIdStr == "bb056e51-f833-4bbc-a8f1-eecfe62bf9c9"] <- "FALSE"
AR_updated$If.the.age.of.the.cell.line.is.reported..hwo.many.cell.passages.where.driven.[AR_updated$StudyIdStr == "64a785ad-4714-4ec6-b731-4d66128602f7" |
                                                                                           AR_updated$StudyIdStr == "a47caeef-d61b-4b62-9781-d0a59f1bdf0c"] <- 15
AR_updated$If.the.age.of.the.cell.line.is.reported..hwo.many.cell.passages.where.driven.[AR_updated$StudyIdStr == "9cb3f84c-86a8-42e7-9237-77684a4af6ec" |
                                                                                           AR_updated$StudyIdStr == "0bcb0299-c385-4142-93f6-6ee7d5767ee6"] <- 10
AR_updated$If.the.age.of.the.cell.line.is.reported..hwo.many.cell.passages.where.driven.[AR_updated$StudyIdStr == "4b98515f-057a-4d17-af4a-527e9b58ec36"] <- 3
AR_updated$How.is.the.glucose.level.of.the.DMEM.medium.for.U87.[AR_updated$StudyIdStr == "9cb3f84c-86a8-42e7-9237-77684a4af6ec" |
                                                                  AR_updated$StudyIdStr == "c2320d09-2c4f-413e-b3e3-a2f5a896a357" |
                                                                  AR_updated$StudyIdStr == "e8c88cd2-830a-4b62-95e8-18fdeaa25347" |
                                                                  AR_updated$StudyIdStr == "6f3a6613-897c-4bb6-9f5a-a510bf1ec6ad" |
                                                                  AR_updated$StudyIdStr == "bb056e51-f833-4bbc-a8f1-eecfe62bf9c9" |
                                                                  AR_updated$StudyIdStr == "a47caeef-d61b-4b62-9781-d0a59f1bdf0c"] <- "high glucose & pyruvate"
AR_updated$How.is.the.glucose.level.of.the.DMEM.medium.for.U87.[AR_updated$StudyIdStr == "0bcb0299-c385-4142-93f6-6ee7d5767ee6"] <- "no glucose"
AR_updated$Number.of.animals.in.cohort[AR_updated$StudyIdStr == 1] <- 1

# For the analysis of the parameter phenotypes and meta-analysis the dataframe "AR_updated" is used.
# For the analysis of the reporting quality the dataframe "AR" is used.




#### Preparations for meta-analysis

# Adding a row number column to AR dataframe:
AR_updated$rn <- row_number(AR_updated$StudyIdStr)

# Selecting outcome data columns (AR_OD is abbr. for "Annotations reconciled outcome data"):
AR_OD <- AR_updated %>% select_(1,4,5,6,7,8,9,10,11,58,29,60,28,15)

# Filtering the Control data:
AR_OD_C <- AR_OD %>% filter(tmz_concentration_in_microM_numeric == 0)

# Filtering the tmz data:
AR_OD_T <- AR_OD %>% filter(!is.element(AR_OD$rn, AR_OD_C$rn))

# Correcting data that were not plausible and checked again in the original article:
AR_OD_T$Average[AR_OD_T$StudyIdStr == "de2d96d2-bbff-462a-b344-fa179031f126" & AR_OD_T$tmz_concentration_in_microM_numeric == 400 & AR_OD_T$Time == 2] <- 71.8
AR_OD_T$Greater.is.worse.[AR_OD_T$StudyIdStr == "b1c376e5-7358-4b31-b62c-d878e7c9ce68"] <- TRUE
AR_OD_T$Greater.is.worse.[AR_OD_T$StudyIdStr == "65125fda-7b0f-44e5-8f71-2505146d4f2e" & AR_OD_T$tmz_concentration_in_microM_numeric == 100 & AR_OD_T$Time == 3] <- FALSE
AR_OD_T$Average[AR_OD_T$StudyIdStr == "d033365b-a968-4559-ab12-2ea17e5908e1"] <- 0.247
AR_OD_C$Average[AR_OD_C$StudyIdStr == "d033365b-a968-4559-ab12-2ea17e5908e1"] <- 0.383



### Transformung standard error of the mean (SEM) into standard deviation (SD)

# Filtering the AR_OD for experiments with errortype = SEM:
AR_OD %>% filter(Outcome.measure.error.type == "SEM") %>% nrow()

# Adding a column with the standard deviation for each outcome data row to AR_OD:
AR_OD$SD <- ifelse(AR_OD$Outcome.measure.error.type == "SD", AR_OD$Error, AR_OD$Error * sqrt(AR_OD$Number.of.animals.in.cohort))

# If SD = 0, print NA because then the Number of experiments is zero (which means it is not reported):
AR_OD$SD <- ifelse(AR_OD$SD == 0, NA, AR_OD$SD) 

# Checking if there are any data rows left with SD as "NA" and the number of experimente != 0:
AR_OD %>% filter(is.na(SD) == TRUE) %>% filter(Number.of.animals.in.cohort != 0)

# as they are all control or irrelavant data we continue
# the two articles with NA values in the outcome data are not relevant for the meta-analysis!

# Filtering the control data again:
AR_OD_C <- AR_OD %>% filter(tmz_concentration_in_microM_numeric == 0)

# Filtering the tmz data again:
AR_OD_T <- AR_OD %>% filter(!is.element(AR_OD$rn, AR_OD_C$rn))

# Adding the control row number to AR_OD_C:
AR_OD_C$rn_C <- row_number(AR_OD_C$StudyIdStr)



### Dealing with missing error data for control cell viability values:

# If the Error of the control values is 0, the mean error of the TMZ values of the same experiment gets taken as the missing error value of control as the best possible assumption:

# For how many values does this matter?:
AR_OD_C %>% filter(Error == 0) %>% nrow()

# Creating a new dataframe for these values:
AR_OD_C_E0 <- AR_OD_C %>% filter(Error == 0)

# Adding a column with the joined StudyID, ExpLabel and Time to the control data:
AR_OD_C_E0$joined_StudyID_EL_Time <- paste(AR_OD_C_E0$StudyIdStr, AR_OD_C_E0$ExperimentLabel, AR_OD_C_E0$Time)

# Adding a column with the joined StudyID, ExpLabel and Time to the treatment data:
AR_OD_T$joined_StudyID_EL_Time <- paste(AR_OD_T$StudyIdStr, AR_OD_T$ExperimentLabel, AR_OD_T$Time)

# getting the errors of the TMZ outcome data in the experiments where the error value of control is 0:
Listce1 <- list()
for (c in 1:nrow(AR_OD_C_E0)) {
  dummyce <- AR_OD_T %>% filter(is.element(AR_OD_T$joined_StudyID_EL_Time, AR_OD_C_E0[c,17])) %>% select_(15) %>% print()
  Listce1[[length(Listce1)+1]] = dummyce
}
# Here is one study with NA values because the number of exp. is 0!

# Getting the means of errors of the TMZ outcome data in the experiments where the error of control is 0:
Listce2 <- list()
for (d in 1:length(Listce1)) {
  dummyce2 <- Listce1[d] %>% unlist() %>% mean() %>% print()
  Listce2[[length(Listce2)+1]] = dummyce2
}

# Checking if the number of the calculated errors is the number of errors = 0:
length(Listce2) == nrow(AR_OD_C_E0)
# Should be TRUE!

# converting Listce2 into a dataframe with one column:
ace1 <- data.frame(Listce2)
bce <- data.frame(t(ace1))

# Adding the joined StudyID and ExpLabel and the rn column:
bce$joined_StudyID_EL <- AR_OD_C_E0$joined_StudyID_EL
bce$rn_C <- AR_OD_C_E0$rn_C
bce$rn <- AR_OD_C_E0$rn
# Manually checked some random values -> they are correct!

# Replacing the SDs = 0 with the calculated SDs in the AR_OD_C dataframe:
for (cei in 1:nrow(bce)) {
  AR_OD_C[bce[cei, 3], 15] <- bce[cei, 1]
}

# Replacing the SDs = 0 with the calculated SDs in the overall AR_OD dataframe:
for (ceio in 1:nrow(bce)) {
  AR_OD[bce[ceio, 4], 15] <- bce[ceio, 1]
}




#### Exclusions from Meta-analysis

# Preliminary exclusions (as they do not match the inclusion criteria):
AR_OD_MA <- AR_OD %>% filter(StudyIdStr != "8348eaa2-160b-49e9-9a8a-cab250a0a730" & StudyIdStr != "1ef05e18-9607-4927-9137-69fd633b2773")
# reason for exclusion for study 8348eaa2-160b-49e9-9a8a-cab250a0a730: no outcome data of relevant experiments are presented
# reason for exclusion for study 1ef05e18-9607-4927-9137-69fd633b2773 no usable outcome data of relevant experiments are presented



### Number of experiments not reported

# For the first condition we filter the experiments that have number of experiments = 0 as these studies that do not report the number of experiments:

# Creating a joined StudyID and ExperimentLabel column:
AR_OD$joined_StudyID_EL <- paste(AR_OD$StudyIdStr, AR_OD$ExperimentLabel)

# Creating the new dataframe "AR_OD_ENNNR" (abbr. for annotations reconciled outcome data experiment number not reported):
AR_OD_ENNR <- AR_OD %>% filter(Number.of.animals.in.cohort == 0) %>% select_(16) %>% unique()

# Adding a column with the joined StudyID and ExperimentLabel to AR_OD_MA:
AR_OD_MA$joined_StudyID_EL <- paste(AR_OD_MA$StudyIdStr, AR_OD_MA$ExperimentLabel)

# Excluding them from the AR_OD_MA:
AR_OD_MA <- AR_OD_MA %>% filter(!is.element(AR_OD_MA$joined_StudyID_EL, AR_OD_ENNR$joined_StudyID_EL))

# Checking if there are data left in AR_OD_MA that have the number of experiments = 0:
AR_OD_MA %>% filter(Number.of.animals.in.cohort == 0) %>% nrow()
# as this is zero we can continue!



###  Errortype not reported

# The studies with unclear reporting of error values as either SD or SEM are not excluded, but SEM was assumed as the more conservative approach.

# Loading the articles that have the errortype not reported (afer reviewing by at least two independent reviewers and reconcilitation by a third independent reviewer):
AR_ET_nr_S_final <- read.csv("AR_ET_nr_S_final.csv")

# Filtering AR_OD studies if there are an element of the studies that have the errortype not reported:
AR_OD_ETNR <- AR_OD_MA %>% filter(is.element(AR_OD_MA$StudyIdStr, AR_ET_nr_S_final$x))

# Assuming SEM as the more conservative approach leading to a higher variance and therefore a lower weight in the meta-analysis:
AR_OD_MA$SD[is.element(AR_OD_MA$StudyIdStr, AR_OD_ETNR$StudyIdStr)] <- AR_OD_MA$Error[is.element(AR_OD_MA$StudyIdStr, AR_OD_ETNR$StudyIdStr)] * sqrt(AR_OD_MA$Number.of.animals.in.cohort[is.element(AR_OD_MA$StudyIdStr, AR_OD_ETNR$StudyIdStr)])



### Time (= treatment duration) not reported

# The basis for this are the annotation notes during data extraction with SyRF:

# Loading the dataframe with the annotations notes:
annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments <- read.csv("annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments.csv")

# List of experiments where the time of outcome measurement is not reported:
annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments$X = NULL
annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments$joined_StudyID_Exp <- paste(annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments$StudyIdStr, annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments$ExperimentLabel)
annotations_studies_experiments_timenotreported <- annotations_experiments_cohortnotes_studieswithnotes_manually_adjustments %>% filter(Time.not.reported == "x")

# Deduplicating:
annotations_studies_experiments_timenotreported_dedup <- annotations_studies_experiments_timenotreported[!duplicated(annotations_studies_experiments_timenotreported$joined_StudyID_Exp),]

# Filtering the AR_OD experiments if there are an element of the studies that have the tmz exposure time not reported:
AR_OD_tmzETnr <- AR_OD %>% filter(is.element(joined_StudyID_EL, annotations_studies_experiments_timenotreported_dedup$joined_StudyID_Exp))

# How many data rows have the time not reported?
nrow(AR_OD_tmzETnr)

# Excluding them from the AR_OD_MA:
AR_OD_MA <- AR_OD_MA %>% filter(!is.element(AR_OD_MA$joined_StudyID_EL, AR_OD_tmzETnr$joined_StudyID_EL))



### TMZ concentration not reported

# Creating a new dataframe:
AR_OD_MA_TMZCnr <- AR_OD_MA %>% filter(is.na(tmz_concentration_in_microM_numeric))

# Excluding them from the AR_OD_MA:
AR_OD_MA <- AR_OD_MA %>% filter(!is.element(AR_OD_MA$rn, AR_OD_MA_TMZCnr$rn))




#### Combination of belonging together control and intervention (TMZ) cell viability data

# Adding a column to AR_OD_MA with the joined Study_ID, ExpLabel & Time:
AR_OD_MA$joined_StudyID_EL_Time <- paste(AR_OD_MA$joined_StudyID_EL, AR_OD_MA$Time)

# Filtering the control data of AR_OD_MA:
AR_OD_MA_C <- AR_OD_MA %>% filter(tmz_concentration_in_microM_numeric == 0)



### Checking if there is only one control value for every study:

# number of control data per experiment and time:
AR_OD_MA_C %>% select_(17) %>% unique() %>% nrow()
table(AR_OD_MA_C$joined_StudyID_EL_Time) %>% as.data.frame() %>% arrange(desc(Freq))
# There is one study "35e10def-f7d7-4c8c-a02b-922b44510337" that has to control values for the same experiment. 
# This is because they performed two controls (untreated & DMSO). 
# For the meta-analysis we take the DMSO control value as this is the more rare control type

# deleting the other control data of this experiment to get one control value for every experiment:
AR_OD_MA_C <- AR_OD_MA_C %>% filter(rn != 87)

# checking if this was successful:
AR_OD_MA_C %>% select_(17) %>% unique() %>% nrow() == nrow(AR_OD_MA_C)
# Should be TRUE!

#  Filtering the TMZ data of AR_OD_MA:
AR_OD_MA_T <- AR_OD_MA %>% filter(tmz_concentration_in_microM_numeric != 0)

# Checking if the row numbers are correct:
nrow(AR_OD_MA_T) + nrow(AR_OD_MA_C)
nrow(AR_OD_MA)
# CAVE: there is one more row in AR_OD_MA because for one experiment we have two different types of control (see above)!



### Getting the belonging together control values if available:

# Filtering TMZ data that have corresponding control data:
AR_OD_MA_T_CC <- AR_OD_MA_T %>% filter(is.element(AR_OD_MA_T$joined_StudyID_EL_Time, AR_OD_MA_C$joined_StudyID_EL_Time))

# How many tmz data have corresponding control data:
nrow(AR_OD_MA_T_CC)

# Checking if the unique joined StudyID, ExpLabel and Times are equal in control and tmz data:
AR_OD_MA_T_CC %>% select_(17) %>% unique() %>% nrow()
AR_OD_MA_C %>% select_(17) %>% unique() %>% nrow()
# there have to be control values without corresponding relevant treatment values

# Filtering control data that have no corresponding relevant tmz experiment data in AR_OD_MA_T (becauese the experiments are not relevant):
AR_OD_MA_C_CC <- AR_OD_MA_C %>% filter(is.element(joined_StudyID_EL_Time, AR_OD_MA_T_CC$joined_StudyID_EL_Time))

# Again checking if the unique joined StudyID, ExpLabel and Times are equal in control and tmz data:
AR_OD_MA_T_CC %>% select_(17) %>% unique() %>% nrow() == AR_OD_MA_C_CC %>% select_(17) %>% unique() %>% nrow()
# both have the same number of rows so we can continue!

# Adding a column to AR_OD_MA_T_CC with the corresponding control mean data:
for (ccv in 1:nrow(AR_OD_MA_T_CC)) {
  AR_OD_MA_T_CC[ccv,18] <- AR_OD_MA_C_CC %>% filter(joined_StudyID_EL_Time == AR_OD_MA_T_CC[ccv,17]) %>% select_(8)
}
# Adding a column to AR_OD_MA_T_CC with the corresponding control error data:
for (ccv in 1:nrow(AR_OD_MA_T_CC)) {
  AR_OD_MA_T_CC[ccv,19] <- AR_OD_MA_C_CC %>% filter(joined_StudyID_EL_Time == AR_OD_MA_T_CC[ccv,17]) %>% select_(15)
}
# Adding a column to AR_OD_MA_T_CC with the corresponding control number of experiments data:
for (ccv in 1:nrow(AR_OD_MA_T_CC)) {
  AR_OD_MA_T_CC[ccv,20] <- AR_OD_MA_C_CC %>% filter(joined_StudyID_EL_Time == AR_OD_MA_T_CC[ccv,17]) %>% select_(5)
}



### Filtering the experiment data that have no corresponding control values

# Creating the new dataframe AR_OD_MA_T_nCC (abbr. for "annotations reconciled outcome data meta-analysis treatment no corresponding control):
AR_OD_MA_T_nCC <- AR_OD_MA_T %>% filter(!is.element(rn, AR_OD_MA_T_CC$rn))

# Checking if all data are considered:
nrow(AR_OD_MA_T_CC) + nrow(AR_OD_MA_T_nCC) == nrow(AR_OD_MA_T)
# if it is TRUE we can continue!



### Filtering the treatment data that have no corresponding control values but are presented as percentage of control:

# Getting the different outcome units of AR_OD_MA_T_nCC:
table(AR_OD_MA_T_nCC$Outcome.measure.units)



### Adding a column to AR_OD_MA_T_nCC that states if the outcome unit is in relation to control:

# Looked manually in every article again to see if the data were presented in relation to an untreated control:
AR_OD_MA_T_nCC$outcome_unit_in_relation_to_control <- ifelse(AR_OD_MA_T_nCC$Outcome.measure.units == "Cell viability (% of control)" |
                                                               AR_OD_MA_T_nCC$Outcome.measure.units == "Proliferation (% of control)" |
                                                               AR_OD_MA_T_nCC$Outcome.measure.units == "Cell survival (%) / MTT response" |
                                                               AR_OD_MA_T_nCC$Outcome.measure.units == "Cell survival rate (%)" |
                                                               AR_OD_MA_T_nCC$Outcome.measure.units == "Cell viability (%)",
                                                             TRUE, FALSE)
AR_OD_MA_T_nCC$outcome_unit_in_relation_to_control[AR_OD_MA_T_nCC$StudyIdStr == "1e165179-cc1f-4fef-bbdb-97756651583d"] <- FALSE

# Creating a dataframe AR_OD_MA_T_nCC_RtC (abbr. for AR_OD_MA_T_nCC relation to control) with the data that have no corresponding control values but are presented in relation to control:
AR_OD_MA_T_nCC_RtC <- AR_OD_MA_T_nCC %>% filter(outcome_unit_in_relation_to_control == TRUE)

# Adding a column with the current row number:
AR_OD_MA_T_nCC_RtC$rn_rtc <- row_number(AR_OD_MA_T_nCC_RtC$StudyIdStr)



### Getting the mean SD of the intervention (TMZ) values of the same experiment and time as an estimate for the control SD:

# Getting the SDs of the TMZ outcome data in the experiments where the error of control is 0:
Listrtc1 <- list()
for (rtc in 1:nrow(AR_OD_MA_T_nCC_RtC)) {
  dummyrtc <- AR_OD_T %>% filter(is.element(AR_OD_T$joined_StudyID_EL_Time, AR_OD_MA_T_nCC_RtC[rtc,17])) %>% select_(15) %>% print()
  Listrtc1[[length(Listrtc1)+1]] = dummyrtc
}

# Getting the means of SDs of the TMZ outcome data in the experiments where the error of control is 0:
Listrtc2 <- list()
for (rtcm in 1:length(Listrtc1)) {
  dummyrtcm <- Listrtc1[rtcm] %>% unlist() %>% mean() %>% print()
  Listrtc2[[length(Listrtc2)+1]] = dummyrtcm
}

# Checking if the number of the calculated errors is the number of errors = 0:
length(Listrtc2) == nrow(AR_OD_MA_T_nCC_RtC)
# Should be TRUE!

# converting Listrtc2 into a dataframe with one column:
artc1 <- data.frame(Listrtc2)
brtc <- data.frame(t(artc1))

# Adding the joined StudyID and ExpLabel and the rn column:
brtc$joined_StudyID_EL_Time <- AR_OD_MA_T_nCC_RtC$joined_StudyID_EL_Time
brtc$rn_rtc <- AR_OD_MA_T_nCC_RtC$rn_rtc
brtc$rn <- AR_OD_MA_T_nCC_RtC$rn
# Manually checked some random values -> they are correct!

# Checking if the experiments are in an identical order:
all(AR_OD_MA_T_nCC_RtC$joined_StudyID_EL_Time == brtc$joined_StudyID_EL_Time)
# all values should be TRUE!

# Adding the calculated average SD values to the AR_OD_MA_T_nCC_RtC dataframe:
AR_OD_MA_T_nCC_RtC$SD_control <- brtc$t.artc1.



### Getting the control cell viability means

# They can be assumed as either 100% or 1 depending on the scale of the intervention outcome values (manually checked):
AR_OD_MA_T_nCC_RtC$Average_control <- 100
AR_OD_MA_T_nCC_RtC$Average_control[AR_OD_MA_T_nCC_RtC$joined_StudyID_EL_Time == "e3d07e37-7606-446e-bac2-39fbda6d045d e 5"] <- 1

# The number of experiments of control gets assumed as the number of experiments of intervention (TMZ):
AR_OD_MA_T_nCC_RtC$N_control <- AR_OD_MA_T_nCC_RtC$Number.of.animals.in.cohort



### Combining the dataframe AR_OD_MA_T_CC and AR_OD_MA_T_nCC_RTC

colnames(AR_OD_MA_T_CC)
colnames(AR_OD_MA_T_nCC_RtC)

# creating new dataframes with the relevant columns of AR_OD_MA_T_nCC_RTC:
AR_OD_MA_T_nCC_RtC_f <- AR_OD_MA_T_nCC_RtC[,c(1:17,21,20,22)]

# Changing column names of AR_OD_MA_T_nCC_RTC_f:
colnames(AR_OD_MA_T_nCC_RtC_f)[18] <- paste("Average.1")
colnames(AR_OD_MA_T_nCC_RtC_f)[19] <- paste("SD.1")
colnames(AR_OD_MA_T_nCC_RtC_f)[20] <- paste("Number.of.animals.in.cohort.1")

# Changing a column name of AR_OD_MA_T_CC:
colnames(AR_OD_MA_T_CC)[19] <- paste("SD.1")

# Checking if the colnames are identical:
all(colnames(AR_OD_MA_T_CC) == colnames(AR_OD_MA_T_nCC_RtC_f))
# all values should be TRUE

# combining the two dataframes:
AR_OD_MA_T_f <- rbind(AR_OD_MA_T_CC, AR_OD_MA_T_nCC_RtC_f)
# the final dataframe with the relevant outcome data: "AR_OD_MA_T_f" (abbr. for "annotations reconciled outcome data treatment final")




#### Manual corrections of typing mistakes during data extraction (extremely outlying data were all checked again based on the original full-text articles)

# Control average should be 100 and not 1:
AR_OD_MA_T_f$Average.1[AR_OD_MA_T_f$joined_StudyID_EL == "b1c376e5-7358-4b31-b62c-d878e7c9ce68 e"] <- 100

# Typo:
AR_OD_MA_T_f$Average[AR_OD_MA_T_f$StudyIdStr == "de2d96d2-bbff-462a-b344-fa179031f126" & AR_OD_MA_T_f$tmz_concentration_in_microM_numeric == 400 & AR_OD_MA_T_f$Time == 2] <- 71.8

# Wrong outcome unit:
AR_OD_MA_T_f$Greater.is.worse.[AR_OD_MA_T_f$joined_StudyID_EL == "65125fda-7b0f-44e5-8f71-2505146d4f2e e2"] <- FALSE
AR_OD_MA_T_f$Outcome.measure.units[AR_OD_MA_T_f$joined_StudyID_EL == "65125fda-7b0f-44e5-8f71-2505146d4f2e e2"] <- paste("Inhibition Rate (100%)")

# Wrong position of decimal:
AR_OD_MA_T_f$SD.1[AR_OD_MA_T_f$StudyIdStr == "cea51952-de58-4120-8cc1-dc6a2a441661" & AR_OD_MA_T_f$Time == 0] <- 0.065

# Wrong position of decimal:
AR_OD_MA_T_f$SD[AR_OD_MA_T_f$StudyIdStr == "d033365b-a968-4559-ab12-2ea17e5908e1"] <- 4
AR_OD_MA_T_f$SD.1[AR_OD_MA_T_f$StudyIdStr == "d033365b-a968-4559-ab12-2ea17e5908e1"] <- 4

# Wrong position of decimal:
AR_OD_MA_T_f$SD.1[AR_OD_MA_T_f$joined_StudyID_EL_Time == "803a5cdf-ab6b-4f76-b286-35a14c443612 e 3"] <- 0.8660254

# Wrong position of decimal:
AR_OD_MA_T_f$tmz_concentration_in_microM_numeric[AR_OD_MA_T_f$StudyIdStr == "4f7e8b69-fe49-476a-bf2c-b8e84d4229ba" & AR_OD_MA_T_f$CohortLabel == "tmz 1.875 uM"] <- 1.875




##### Parameter phenotypes analysis 

### Checking if the parameters that should be the same for each row of one study are indeed the same for each study:

# Counting how many rows each study has in a new dataframe ARSC:
ARSC <- count_(AR, "StudyIdStr")

# Getting the number of unique values of experiment parameters for each study with a for-loop:
# Every result for each parameter gets stored in a list "List(Number)" with the helping data "dummy(Number)":
# Note that number 6 is missing because one parameter was present twice and removed afterwards.
List1 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy1 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(13) %>% unique() %>% nrow() 
  List1[[length(List1)+1]] = dummy1
}
List2 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy2 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(17) %>% unique() %>% nrow() 
  List2[[length(List2)+1]] = dummy2
}
List3 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy3 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(18) %>% unique() %>% nrow() 
  List3[[length(List3)+1]] = dummy3
}
List4 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy4 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(19) %>% unique() %>% nrow() 
  List4[[length(List4)+1]] = dummy4
}
List5 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy5 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(20) %>% unique() %>% nrow() 
  List5[[length(List5)+1]] = dummy5
}
List7 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy7 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(21) %>% unique() %>% nrow() 
  List7[[length(List7)+1]] = dummy7
}
List8 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy8 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(22) %>% unique() %>% nrow() 
  List8[[length(List8)+1]] = dummy8
}
List9 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy9 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(25) %>% unique() %>% nrow() 
  List9[[length(List9)+1]] = dummy9
}
List10 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy10 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(26) %>% unique() %>% nrow() 
  List10[[length(List10)+1]] = dummy10
}
List11 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy11 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(33) %>% unique() %>% nrow() 
  List11[[length(List11)+1]] = dummy11
}
List12 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy12 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(34) %>% unique() %>% nrow() 
  List12[[length(List12)+1]] = dummy12
}
List13 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy13 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(36) %>% unique() %>% nrow() 
  List13[[length(List13)+1]] = dummy13
}
List14 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy14 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(37) %>% unique() %>% nrow() 
  List14[[length(List14)+1]] = dummy14
}
List15 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy15 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(39) %>% unique() %>% nrow() 
  List15[[length(List15)+1]] = dummy15
}
List16 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy16 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(27) %>% unique() %>% nrow() 
  List16[[length(List16)+1]] = dummy16
}
List17 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy17 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(38) %>% unique() %>% nrow() 
  List17[[length(List17)+1]] = dummy17
}
List18 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy18 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(54) %>% unique() %>% nrow() 
  List18[[length(List18)+1]] = dummy18
}

# Every List gets converted into a dataframe "a(Number)":
a1 <- data.frame(List1)
a2 <- data.frame(List2)
a3 <- data.frame(List3)
a4 <- data.frame(List4)
a5 <- data.frame(List5)
a7 <- data.frame(List7)
a8 <- data.frame(List8)
a9 <- data.frame(List9)
a10 <- data.frame(List10)
a11 <- data.frame(List11)
a12 <- data.frame(List12)
a13 <- data.frame(List13)
a14 <- data.frame(List14)
a15 <- data.frame(List15)
a16 <- data.frame(List16)
a17 <- data.frame(List17)
a18 <- data.frame(List18)

# Creating a helping dataframe "b" with the column/row reversed dataframes of "a"
b <- data.frame(t(a1), t(a2), t(a3), t(a4), t(a5), t(a7), t(a8), t(a9), t(a10), t(a11), t(a12), t(a13), t(a14), t(a15), t(a16), t(a17), t(a18))

# Creating the dataframe ARPM (Annotations reconciled parameter match) with the column specifications for each parameter based on the "b" dataframe
ARPM <- data.frame(ARSC$StudyIdStr)
ARPM$CoI <- b$t.a1.
ARPM$U87otherSource <- b$t.a2.
ARPM$Passages <- b$t.a3.
ARPM$Authentification <- b$t.a4.
ARPM$Mycoplasma <- b$t.a5.
ARPM$U87age <- b$t.a7.
ARPM$Confluency <- b$t.a8.
ARPM$TimeCellPassaging <- b$t.a9.
ARPM$U87source <- b$t.a10.
ARPM$FBSuse <- b$t.a11.
ARPM$otherAntibiotics <- b$t.a12.
ARPM$Antibiotics <- b$t.a13.
ARPM$FBSsource <- b$t.a14.
ARPM$FBSotherSource <- b$t.a15.
ARPM$CellPassagingCriteria <- b$t.a16.
ARPM$Glucose <- b$t.a17.
ARPM$Cell_number_vol_conc <- b$t.a18.

# Creating a dataframe "ARPML" (Annotations reconciled parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARPML <- data.frame(ifelse(ARPM[,2:18] == 1, TRUE, FALSE))
# Adding a Column with the StudyIDs:
ARPML <- add_column(ARPML, "StudyIdStr" <- ARSC$StudyIdStr, .before = "CoI")
# Changing the column name of the studyIDs:
colnames(ARPML)[colnames(ARPML) == "\"StudyIdStr\" <- ARSC$StudyIdStr"] <- "StudyIdStr"

# Checking if the parameters are identical in all studies:
all(ARPML$CoI)
all(ARPML$U87otherSource)
all(ARPML$Passages)
all(ARPML$Authentification)
all(ARPML$Mycoplasma)
all(ARPML$U87age)
all(ARPML$Confluency)
all(ARPML$TimeCellPassaging)
all(ARPML$U87source)
all(ARPML$FBSuse)
all(ARPML$otherAntibiotics)
all(ARPML$Antibiotics)
all(ARPML$FBSsource)
all(ARPML$FBSotherSource)
all(ARPML$CellPassagingCriteria)
all(ARPML$Glucose)
all(ARPML$Cell_number_vol_conc)

## There are four parameters (mycoplasma exclusion, U87 source and other U87 source) with different phenotypes within the same study. 
# This parameters were checked again based on the original articles.

## Corrections after looking into the original articles:
# Parameter: U87otherSource
ARPML %>% filter(U87otherSource != 1) %>% select(StudyIdStr)
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",17]
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",17] <- "Invitrogen Life Technologies, Inc., Carlsbad, CA"
AR[AR$StudyIdStr == "4f7e8b69-fe49-476a-bf2c-b8e84d4229ba",17]
AR[AR$StudyIdStr == "4f7e8b69-fe49-476a-bf2c-b8e84d4229ba",17] <- "China Academia Sinica cell repository (Shanghai, China)"
# Parameter: Mycoplasma
ARPML %>% filter(Mycoplasma != 1) %>% select(StudyIdStr)
AR[AR$StudyIdStr == "65e0c64b-fd79-4f96-b170-ef785d6aa038",20]
AR[AR$StudyIdStr == "65e0c64b-fd79-4f96-b170-ef785d6aa038",20] <- TRUE
# Parameter: U87Source:
ARPML %>% filter(U87source != 1) %>% select(StudyIdStr)
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",26]
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",26] <- "other source"
# Parameter: Glucose
ARPML %>% filter(Glucose != 1) %>% select(StudyIdStr)
AR[AR$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606",38]
# -> this study contains experiments with both high and low glucose! -> no correction needed
# Parameter: number of animals (mistake seen randomly, therefore it gets corrected here):
AR$Number.of.animals.in.cohort[AR$StudyIdStr == "f0c825a5-5083-4867-8af2-38cf451117c6" &
                                 AR$ExperimentLabel == "e2 (tb)"] <- 4

## After corrections: Checking again if the parameters are identical in all studies:


# Getting the number of unique values of experiment parameters for each study with a for-loop:
# every result for each parameter gets storaged in a list "List(Number)" with the helping data "dummy(Number)":
# note that number 6 is missing because I had one parameter duplicate
List1 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy1 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(13) %>% unique() %>% nrow() 
  List1[[length(List1)+1]] = dummy1
}
List2 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy2 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(17) %>% unique() %>% nrow() 
  List2[[length(List2)+1]] = dummy2
}
List3 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy3 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(18) %>% unique() %>% nrow() 
  List3[[length(List3)+1]] = dummy3
}
List4 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy4 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(19) %>% unique() %>% nrow() 
  List4[[length(List4)+1]] = dummy4
}
List5 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy5 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(20) %>% unique() %>% nrow() 
  List5[[length(List5)+1]] = dummy5
}
List7 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy7 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(21) %>% unique() %>% nrow() 
  List7[[length(List7)+1]] = dummy7
}
List8 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy8 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(22) %>% unique() %>% nrow() 
  List8[[length(List8)+1]] = dummy8
}
List9 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy9 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(25) %>% unique() %>% nrow() 
  List9[[length(List9)+1]] = dummy9
}
List10 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy10 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(26) %>% unique() %>% nrow() 
  List10[[length(List10)+1]] = dummy10
}
List11 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy11 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(33) %>% unique() %>% nrow() 
  List11[[length(List11)+1]] = dummy11
}
List12 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy12 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(34) %>% unique() %>% nrow() 
  List12[[length(List12)+1]] = dummy12
}
List13 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy13 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(36) %>% unique() %>% nrow() 
  List13[[length(List13)+1]] = dummy13
}
List14 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy14 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(37) %>% unique() %>% nrow() 
  List14[[length(List14)+1]] = dummy14
}
List15 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy15 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(39) %>% unique() %>% nrow() 
  List15[[length(List15)+1]] = dummy15
}
List16 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy16 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(27) %>% unique() %>% nrow() 
  List16[[length(List16)+1]] = dummy16
}


# Every List gets converted into a dataframe "a(Number)":
a1 <- data.frame(List1)
a2 <- data.frame(List2)
a3 <- data.frame(List3)
a4 <- data.frame(List4)
a5 <- data.frame(List5)
a7 <- data.frame(List7)
a8 <- data.frame(List8)
a9 <- data.frame(List9)
a10 <- data.frame(List10)
a11 <- data.frame(List11)
a12 <- data.frame(List12)
a13 <- data.frame(List13)
a14 <- data.frame(List14)
a15 <- data.frame(List15)
a16 <- data.frame(List16)

# Creating a helping dataframe "b" with the column/row reversed dataframes of "a"
b <- data.frame(t(a1), t(a2), t(a3), t(a4), t(a5), t(a7), t(a8), t(a9), t(a10), t(a11), t(a12), t(a13), t(a14), t(a15), t(a16))

# Creating the dataframe ARPM (Annotations reconciled parameter match) with the column specifications for each parameter based on the "b" dataframe
ARPM <- data.frame(ARSC$StudyIdStr)
ARPM$CoI <- b$t.a1.
ARPM$U87otherSource <- b$t.a2.
ARPM$Passages <- b$t.a3.
ARPM$Authentification <- b$t.a4.
ARPM$Mycoplasma <- b$t.a5.
ARPM$U87age <- b$t.a7.
ARPM$Confluency <- b$t.a8.
ARPM$TimeCellPassaging <- b$t.a9.
ARPM$U87source <- b$t.a10.
ARPM$FBSuse <- b$t.a11.
ARPM$otherAntibiotics <- b$t.a12.
ARPM$Antibiotics <- b$t.a13.
ARPM$FBSsource <- b$t.a14.
ARPM$FBSotherSource <- b$t.a15.
ARPM$CellPassagingCriteria <- b$t.a16.

# Creating a dataframe "ARPML" (Annotations reconciled parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARPML <- data.frame(ifelse(ARPM[,2:16] == 1, TRUE, FALSE))
# Adding a Column with the StudyIDs:
ARPML <- add_column(ARPML, "StudyIdStr" <- ARSC$StudyIdStr, .before = "CoI")
# Changing the column name of the studyIDs:
colnames(ARPML)[colnames(ARPML) == "\"StudyIdStr\" <- ARSC$StudyIdStr"] <- "StudyIdStr"

# Checking if the parameters are identical in all studies:
all(ARPML$CoI)
all(ARPML$U87otherSource)
all(ARPML$Passages)
all(ARPML$Authentification)
all(ARPML$Mycoplasma)
all(ARPML$U87age)
all(ARPML$Confluency)
all(ARPML$TimeCellPassaging)
all(ARPML$U87source)
all(ARPML$FBSuse)
all(ARPML$otherAntibiotics)
all(ARPML$Antibiotics)
all(ARPML$FBSsource)
all(ARPML$FBSotherSource)
all(ARPML$CellPassagingCriteria)
# as all values despite glucose concentration (where actually two different phenotypes are used) are TRUE, we can continue!



### Parameter phenotype frequencies

# Creating a dataframe with only one row per article and the parameters that were checked for equality before:
# ARIP1R is an abbreviation for "annotations reconciled identical parameters one row"
ARIP1R <- AR_updated[, c(1,13,17:22,25:27,33:34,36:37,39,54)] %>% distinct(StudyIdStr, .keep_all = T)

## Frequencies of selected study parameter phenotypes:
table(ARIP1R$How.is.the.reporting.of.conflicts.of.interest.)
table(ARIP1R$Was.a.cell.line.authentification.of.U87.cells.performed.)
table(ARIP1R$Is.a.successfull.mycoplasma.contamination.exclusion.of.the.cell.line.reported.)
table(ARIP1R$Is.the.age.of.the.cell.line.reported.)
table(ARIP1R$If.the.age.of.the.cell.line.is.reported..hwo.many.cell.passages.where.driven.)
table(ARIP1R$At.which.confluency.level.did.they.perform.a.cell.passaging.)
table(ARIP1R$How.long.is.is.the.time.period.until.they.performed.a.cell.passaging.)
table(ARIP1R$How.ist.the.reporting.on.cell.passaging.criteria.)
table(ARIP1R$From.whom.did.they.get.the.U87.cells.)
table(ARIP1R$What.other.source.for.U87.cells.did.they.have.)
table(ARIP1R$Did.they.use.FCS.FBS.medium.for.U87.)
table(ARIP1R$Which.antibiotic.did.they.used.in.the.U87.cell.line.medium.)
table(ARIP1R$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it.)
table(ARIP1R$Is.the.number...concentration...volume.the.U87.cells.are.plated.reported.)


## Frequency of combined U87 sources (including other sources):

# Getting the current U87 source phenotypes:
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(26) %>% table()
# ATCC and "not reported" are already complete

# Adding another column to AR_updated with the combined U87 source:
AR_updated$combined_U87_source <- ifelse(AR_updated$From.whom.did.they.get.the.U87.cells. == "ATCC" | 
                                           AR_updated$From.whom.did.they.get.the.U87.cells. == "not reported" |
                                           AR_updated$From.whom.did.they.get.the.U87.cells. == "colleagues",
                                         AR_updated$From.whom.did.they.get.the.U87.cells., NA)

# Getting the current "other U87 sources":
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(17) %>% table()

# Dealing with the "colleagues" as source for U-87 MG cells:
AR_updated$combined_U87_source <- ifelse(AR_updated$What.other.source.for.U87.cells.did.they.have. == "Dr. Adriana da Silva Santos Duarte from Hemocenter, State University of Campinas, Campinas, SÃƒÂƒÃ‚Â£o Paulo, Brazil" |
                                           AR_updated$What.other.source.for.U87.cells.did.they.have. == "Dr. J. Ponten, University of Uppsala, Sweden" |
                                           AR_updated$What.other.source.for.U87.cells.did.they.have. == "Dr. R. Pieper, University of California at San Francisco",
                                         "colleagues",
                                         AR_updated$combined_U87_source)

# current combined U87 sources:
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(61) %>% table()

# dealing with the "Cell Bank of the Chinese Academy of Sciences (Shanghai, China)" as "other U87 source":
toMatch_U87_source_Chinese <- c("Shanghai","Chinese")
AR_updated$combined_U87_source <- ifelse(grepl(paste(toMatch_U87_source_Chinese,collapse="|"), AR_updated$What.other.source.for.U87.cells.did.they.have.),
                                         "Cell Bank of the Chinese Academy of Sciences (Shanghai, China)",
                                         AR_updated$combined_U87_source)

AR_updated$combined_U87_source <- ifelse(is.na(AR_updated$combined_U87_source), "other", AR_updated$combined_U87_source)

# final combined U87 sources parameter phenotypes count:
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(61) %>% table()


## Frequency of combined fetal bovine serum (FBS) sources (including other sources):

# Overview about "normal" FBS sources:
table(AR_updated$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it.)

# Overview about "other" FBS sources:
table(AR_updated$Which.other.source.did.they.have.for.FCS.FBS.)

# Adding a new column to AR_updated with the final combined FBS source:
AR_updated$combined_FBS_source <- ifelse(AR_updated$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it. == "Sigma",
                                         "Sigma-Aldrich, US", ifelse(
                                           AR_updated$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it. == "not reported",
                                           "not reported", ifelse(
                                             AR_updated$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it. == "Gibco & Therma",
                                             "Thermo Fischer Scientific, US", ifelse(
                                               AR_updated$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it. == "",
                                               "No use of FBS reported", NA
                                             )
                                           )
                                         ))

# Adding the "other" FBS sources:
AR_updated$combined_FBS_source <- ifelse(!is.na(AR_updated$combined_FBS_source), AR_updated$combined_FBS_source,
                                         ifelse(grepl("lone", AR_updated$Which.other.source.did.they.have.for.FCS.FBS.), "Hyclone Laboratories Inc, US",
                                                ifelse(grepl("Invitrogen", AR_updated$Which.other.source.did.they.have.for.FCS.FBS.), "Thermo Fischer Scientific, US",
                                                       ifelse(grepl("Techno", AR_updated$Which.other.source.did.they.have.for.FCS.FBS.), "Thermo Fischer Scientific, US",
                                                              ifelse(grepl("Thermo Fisher", AR_updated$Which.other.source.did.they.have.for.FCS.FBS.), "Thermo Fischer Scientific, US", AR_updated$combined_FBS_source)))))

# Adding the "other" other sources:
AR_updated$combined_FBS_source <- ifelse(is.na(AR_updated$combined_FBS_source), "Other source", AR_updated$combined_FBS_source)

# Overview of the distinct combined FBS sources:
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(62) %>% table()


## Frequency of combined antibiotics (including other antibiotics):

# changing the antibiotics from "other" to "not reported" because in these articles they did not name the antibiotics used (after manually checking again based on the original articles):
AR_updated$Which.antibiotic.did.they.used.in.the.U87.cell.line.medium.[AR_updated$StudyIdStr == "3c41250d-a6ef-47cd-bdce-1323921aa009" |
                                                                         AR_updated$StudyIdStr == "828c1097-b791-43fe-9bf8-3cc8b020bdd7" |
                                                                         AR_updated$StudyIdStr == "cc7b32c6-5485-4ed8-90cf-562706fa6955"] <- "not reported"

# Frequency of combined antibiotics:
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(36) %>% table()


## Frequency of glucose levels:

# frequency of the glucose level parameters
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(38) %>% table()

# which parameter is attributed to the study with both high and low glucose? (5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606)
AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% filter(StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606") %>% select_(38)
# so we have to subtract one "high glucose & pyruvate" and add one "both high and low glucose"

# creating a df with one row per study and glucose level:
AR1Rglucose <- AR_updated %>% distinct(StudyIdStr, .keep_all = T) %>% select_(1,38) %>% as.data.frame()
AR1Rglucose <- AR1Rglucose[order(AR1Rglucose$StudyIdStr),]


## Frequency of outcome assessment methods:

# Adding a new column to the AR_updated dataframe with the combined assays (first based on the "Which cell growth/viability assay is used?" column):
AR_updated$assay_combined <- ifelse(
  AR_updated$Which.cell.growth.viability.assay.is.used. == "Alamar-blue assay (colorimetric assay)", "Alamar Blue, colorimetric",
  ifelse(
    AR_updated$Which.cell.growth.viability.assay.is.used. == "CCK8 (Cell counting)", "CCK8, colorimetric",
    ifelse(
      AR_updated$Which.cell.growth.viability.assay.is.used. == "MTS (colorimetric assay)", "MTS, colorimetric",
      ifelse(
        AR_updated$Which.cell.growth.viability.assay.is.used. == "MTT (3-(4,5-dimethylthiazol-2-yl)-2,5-diphenyltetrazolium bromide; colorimetric assay)", "MTT, colorimetric",
        ifelse(
          AR_updated$Which.cell.growth.viability.assay.is.used. == "Trypan Blue exclusion test (Cell counting)", "Trypan Blue exclusion, cell counting",
          ifelse(
            AR_updated$Which.cell.growth.viability.assay.is.used. == "XTT (colorimetric assay)", "other",
            NA)
        )
      )
    )
  )
)

# Adding the values from the "other assays":
AR_updated$assay_combined <- ifelse(
  !is.na(AR_updated$assay_combined), AR_updated$assay_combined,
  ifelse(
    AR_updated$Which.other.outcome.measure.method.is.used. == " sulforhodamine B (SRB) colorimetric assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == " Sulforhodamine B (SRB) colorimetric method" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "SRB assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "SRB colorimetric assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "sulforhodamine B (SRB) assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "Sulforhodamine B (SRB) assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "sulforhodamine B (SRB; Sigma-Aldrich) colorimetric assay" |
      AR_updated$Which.other.outcome.measure.method.is.used. == "Sulforhodamine B assay (SRB)", "Sulforhodamine B (SRB), colorimetric",
    ifelse(
      AR_updated$Which.other.outcome.measure.method.is.used. == "WST-1 (4-[3-(4- iodophenyl)-2-(4-nitrophenyl)-2H-5-tetrazolio]-1,3-benzene disulphonate) colorimetric assay" |
        AR_updated$Which.other.outcome.measure.method.is.used. == "WST-1 assay" |
        AR_updated$Which.other.outcome.measure.method.is.used. == "WST-1 assay (Dojin Kagaku Corp., Kumamoto, Japan)" |
        AR_updated$Which.other.outcome.measure.method.is.used. == "WST-1 cell viability assay" |
        AR_updated$Which.other.outcome.measure.method.is.used. == "WST-1 cell viability assay kit V" |
        AR_updated$Which.other.outcome.measure.method.is.used. == "water soluble tetrazolium salts 1 (WSTÃ¢â‚¬â€˜1) assay", "WST-1 assay, colorimetric",
      "other")
  )
)

# Similar procedure as above (checking how many different outcome assessment methods are present in each study):
List18 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy18 <- AR_updated %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(63) %>% unique() %>% nrow() 
  List18[[length(List18)+1]] = dummy18
}

a18 <- data.frame(List18)

b18 <- data.frame(t(a18)) 

# Adding a column to ARPM that tells how many different outcome assessment methods are present in each study:
ARPM$assay <- b18$t.a18.

# How many studies have used multiple assay types?
ARPM %>% filter(assay > 1) %>% nrow()

# Studies with only one outcome assay:
Studies_with_one_assay <- ARPM %>% filter(assay == 1) %>% select_(1)

# Getting the assays of the studies with only one outcome assay:
AR_updated %>% filter(is.element(StudyIdStr, Studies_with_one_assay$ARSC.StudyIdStr)) %>% distinct(StudyIdStr, .keep_all = T) %>% select_(63) %>% table()
# CAVE: Some articles use multiple outcome assessment methods, these are checked manually in the csv file below.


## Frequency of Temozolomide parameter phenotypes:

# Checking if the TMZ source parameter is identical over all rows belonging to a particular study with Control == FALSE:

# Getting the unique combinations of StudyID and Control in a dataframe "ARCC" (abbreviation for annotations reconciled control count):
ARCC <- ddply(AR_updated, .(AR_updated$StudyIdStr, AR_updated$Control), nrow)

# Creating a new dataframe ARCFC (abbreviation for annotations reconciled control false count) with ARCC filtered by Control == FALSE:
ARCFC <- ARCC %>% filter(`AR_updated$Control` == FALSE)

# Adding a column to ARCFC with the joined StudyIDs and Control values:
ARCFC$joined_StudyID_Control <- paste(ARCFC$`AR_updated$StudyIdStr`, ARCFC$`AR_updated$Control`)

# Adding a column to AR with the joined StudyIDs, ExperimentLabels and InterventionLabels:
AR$joined_StudyID_Control <- paste(AR$StudyIdStr, AR$Control)
AR_updated$joined_StudyID_Control <- paste(AR_updated$StudyIdStr, AR_updated$Control)

# Getting the number of unique values of intervention parameters for each study with a for-loop:

# Every result for each parameter gets storaged in a list "Listc(Number)" with the helping data "dummyc(Number)":
Listc1 <- list()
for (o in seq(1, nrow(ARCFC), 1)) {
  dummyc1 <- AR_updated %>% filter(joined_StudyID_Control == ARCFC[o,4]) %>% select_(47) %>% unique() %>% nrow() 
  Listc1[[length(Listc1)+1]] = dummyc1
}
Listc2 <- list()
for (o in seq(1, nrow(ARCFC), 1)) {
  dummyc2 <- AR_updated %>% filter(joined_StudyID_Control == ARCFC[o,4]) %>% select_(45) %>% unique() %>% nrow() 
  Listc2[[length(Listc2)+1]] = dummyc2
}
Listc3 <- list()
for (o in seq(1, nrow(ARCFC), 1)) {
  dummyc3 <- AR_updated %>% filter(joined_StudyID_Control == ARCFC[o,4]) %>% select_(41) %>% unique() %>% nrow() 
  Listc3[[length(Listc3)+1]] = dummyc3
}

# Every List gets converted into a dataframe "ac(Number)":
ac1 <- data.frame(Listc1)
ac2 <- data.frame(Listc2)
ac3 <- data.frame(Listc3)

# Creating a helping dataframe "bc" with the column/row reversed dataframes of "ac"
bc <- data.frame(t(ac1), t(ac2), t(ac3))

# Creating the dataframe ARCPM (Annotations reconciled control parameter match) with the column specifications for each parameter based on the "bc" dataframe
ARCPM <- data.frame(ARCFC$joined_StudyID_Control)
ARCPM$TMZ_source <- bc$t.ac1.
ARCPM$TMZ_vol <- bc$t.ac2.
ARCPM$TMZ_vol_number <- bc$t.ac3.

# Creating a dataframe "ARCPML" (Annotations reconciled control parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARCPML <- data.frame(ifelse(ARCPM[,2:4] == 1, TRUE, FALSE))

# Checking if the parameter is identical in each treatment group row of each study:
all(ARCPML$TMZ_source) 
all(ARCPML$TMZ_vol)
all(ARCPML$TMZ_vol_number)

# Correcting some TMZ_source and TMZ_Volume parameters after checking again with looking into the original articles:
AR$Is.the.source.of.temozolomide.reported.[AR$StudyIdStr == "1e165179-cc1f-4fef-bbdb-97756651583d"] <- FALSE
AR$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR$StudyIdStr == "1e165179-cc1f-4fef-bbdb-97756651583d"] <- FALSE
AR$Is.the.source.of.temozolomide.reported.[AR$StudyIdStr == "3eedd5d4-fa1a-4ae8-ad65-8cab83302689" |
                                             AR$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606" |
                                             AR$StudyIdStr == "9deed9c3-3f5b-4d43-89d2-22f128a03258" |
                                             AR$StudyIdStr == "b3599eec-a611-4133-8030-42371542a2ef" |
                                             AR$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad" |
                                             AR$StudyIdStr == "9a9433d1-7245-4510-9f56-8032ad295081" |
                                             AR$StudyIdStr == "db764b96-86d1-4aee-961e-b85e0cda55cf"] <- TRUE
AR$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR$StudyIdStr == "3eedd5d4-fa1a-4ae8-ad65-8cab83302689" |
                                                                  AR$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606" |
                                                                  AR$StudyIdStr == "9deed9c3-3f5b-4d43-89d2-22f128a03258" |
                                                                  AR$StudyIdStr == "b3599eec-a611-4133-8030-42371542a2ef" |
                                                                  AR$StudyIdStr == "9a9433d1-7245-4510-9f56-8032ad295081" |
                                                                  AR$StudyIdStr == "db764b96-86d1-4aee-961e-b85e0cda55cf"] <- FALSE
AR$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad"] <- TRUE
AR$What.is.the.volume.of.the.added.Temozolomide.Suspension.[AR$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad"] <- 200

# Also correcting it in AR_updated:
AR_updated$Is.the.source.of.temozolomide.reported.[AR_updated$StudyIdStr == "1e165179-cc1f-4fef-bbdb-97756651583d"] <- FALSE
AR_updated$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR_updated$StudyIdStr == "1e165179-cc1f-4fef-bbdb-97756651583d"] <- FALSE
AR_updated$Is.the.source.of.temozolomide.reported.[AR_updated$StudyIdStr == "3eedd5d4-fa1a-4ae8-ad65-8cab83302689" |
                                                     AR_updated$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606" |
                                                     AR_updated$StudyIdStr == "9deed9c3-3f5b-4d43-89d2-22f128a03258" |
                                                     AR_updated$StudyIdStr == "b3599eec-a611-4133-8030-42371542a2ef" |
                                                     AR_updated$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad" |
                                                     AR_updated$StudyIdStr == "9a9433d1-7245-4510-9f56-8032ad295081" |
                                                     AR_updated$StudyIdStr == "db764b96-86d1-4aee-961e-b85e0cda55cf"] <- TRUE
AR_updated$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR_updated$StudyIdStr == "3eedd5d4-fa1a-4ae8-ad65-8cab83302689" |
                                                                          AR_updated$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606" |
                                                                          AR_updated$StudyIdStr == "9deed9c3-3f5b-4d43-89d2-22f128a03258" |
                                                                          AR_updated$StudyIdStr == "b3599eec-a611-4133-8030-42371542a2ef" |
                                                                          AR_updated$StudyIdStr == "9a9433d1-7245-4510-9f56-8032ad295081" |
                                                                          AR_updated$StudyIdStr == "db764b96-86d1-4aee-961e-b85e0cda55cf"] <- FALSE
AR_updated$Is.the.volume.of.the.added.Temozolomide.suspension.reported.[AR_updated$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad"] <- TRUE
AR_updated$What.is.the.volume.of.the.added.Temozolomide.Suspension.[AR_updated$StudyIdStr == "cd9cc8b6-ff29-4acd-a1d6-4b1c4764daad"] <- 200

# Checking if there are NA values in the TMZ values:
AR_updated %>% filter(Control == FALSE) %>% filter(is.na(Is.the.source.of.temozolomide.reported.) == TRUE) %>% select_(1)
# if there are no data remaining we can continue!

# Creating a dataframe with only one row per StudyID Control == FALSE:
# ARSCF1R is an abbreviation for "annotations reconciled study control false one row"
ARSCF1R <- AR_updated[, c(60,59,47,45,41,1,64)] %>% filter(Control == FALSE) %>% distinct(joined_StudyID_Control,.keep_all = T)

# Frequencies of temozolomide parameter phenotypes:
table(ARSCF1R$Is.the.source.of.temozolomide.reported.)
table(ARSCF1R$Is.the.volume.of.the.added.Temozolomide.suspension.reported.)
table(ARSCF1R$What.is.the.volume.of.the.added.Temozolomide.Suspension.)

# creating a df with one row per study and tmz source and tmz volume:
AR1Rtmz_source <- AR_updated %>% distinct(joined_StudyID_Control, .keep_all = TRUE) %>% select_(1,47,59,45) %>% filter(Control == FALSE)
AR1Rtmz_source <- AR1Rtmz_source[order(AR1Rtmz_source$StudyIdStr),]

# How many articles have information in the AR1Rtmz_source dataframe?
nrow(AR1Rtmz_source)
# Five articles informations are missing

# Manually adding the missing information based on the extracted data (plus checking again with looking in the original research articles:
AR1Rtmz_source[133,1:4] <- c("62452ec0-afbf-4afd-b275-20070bbf2e65", TRUE, FALSE, FALSE)
AR1Rtmz_source[134,1:4] <- c("9a9433d1-7245-4510-9f56-8032ad295081", TRUE, FALSE,	FALSE)
AR1Rtmz_source[135,1:4] <- c("a5eea77e-d942-42de-9d14-318ea8c5a802", TRUE, FALSE, FALSE)
AR1Rtmz_source[136,1:4] <- c("afeaea9e-1bb7-4314-881a-ea9811b22186", TRUE, FALSE, FALSE)
AR1Rtmz_source[137,1:4] <- c("db764b96-86d1-4aee-961e-b85e0cda55cf", TRUE, FALSE, FALSE)

# Ordering again by StudyID:
AR1Rtmz_source <- AR1Rtmz_source[order(AR1Rtmz_source$StudyIdStr),]


## Frequency of control parameter phenotypes

# Earlier corrections of control parameters of one Study (necessary because this study had no reconciled control parameters before, checked based on the original research articles):
AR$What.is.the.type.of.control.treatment.[AR$StudyIdStr == "f72cceb0-db8a-40e0-b989-887271dc61ec"] <- "not reported"
AR$Is.the.volume.of.the.added.control.suspension.reported.[AR$StudyIdStr == "f72cceb0-db8a-40e0-b989-887271dc61ec"] <- FALSE
AR_updated$What.is.the.type.of.control.treatment.[AR_updated$StudyIdStr == "f72cceb0-db8a-40e0-b989-887271dc61ec"] <- "not reported"
AR_updated$Is.the.volume.of.the.added.control.suspension.reported.[AR_updated$StudyIdStr == "f72cceb0-db8a-40e0-b989-887271dc61ec"] <- FALSE

# Creating a new dataframe ARCTC (abbreviation for annotations reconciled control true count) with ARCC filtered by Control == TRUE:
ARCTC <- ARCC %>% filter(`AR_updated$Control` == TRUE)

# Adding a column to ARCTC with the joined StudyIDs and Control values:
ARCTC$joined_StudyID_Control <- paste(ARCTC$`AR_updated$StudyIdStr`, ARCTC$`AR_updated$Control`)

# Getting the number of unique values of intervention parameters for each study with a for-loop:

# Every result for each parameter gets storaged in a list "Listc(Number)" with the helping data "dummyc(Number)":
Listc11 <- list()
for (p in seq(1, nrow(ARCTC), 1)) {
  dummyc11 <- AR_updated %>% filter(joined_StudyID_Control == ARCTC[p,4]) %>% select_(49) %>% unique() %>% nrow() 
  Listc11[[length(Listc11)+1]] = dummyc11
}
Listc12 <- list()
for (p in seq(1, nrow(ARCTC), 1)) {
  dummyc12 <- AR_updated %>% filter(joined_StudyID_Control == ARCTC[p,4]) %>% select_(44) %>% unique() %>% nrow() 
  Listc12[[length(Listc12)+1]] = dummyc12
}
Listc13 <- list()
for (p in seq(1, nrow(ARCTC), 1)) {
  dummyc13 <- AR_updated %>% filter(joined_StudyID_Control == ARCTC[p,4]) %>% select_(57) %>% unique() %>% nrow() 
  Listc13[[length(Listc13)+1]] = dummyc13
}
Listc14 <- list()
for (p in seq(1, nrow(ARCTC), 1)) {
  dummyc14 <- AR_updated %>% filter(joined_StudyID_Control == ARCTC[p,4]) %>% select_(52) %>% unique() %>% nrow() 
  Listc14[[length(Listc14)+1]] = dummyc14
}

# Every List gets converted into a dataframe "ac(Number)":
ac11 <- data.frame(Listc11)
ac12 <- data.frame(Listc12)
ac13 <- data.frame(Listc13)
ac14 <- data.frame(Listc14)

# Creating a helping dataframe "bct" with the column/row reversed dataframes of "ac"
bct <- data.frame(t(ac11), t(ac12), t(ac13), t(ac14))

# Creating the dataframe ARCTPM (Annotations reconciled control true parameter match) with the column specifications for each parameter based on the "bct" dataframe
ARCTPM <- data.frame(ARCTC$joined_StudyID_Control)
ARCTPM$ControlType <- bct$t.ac11.
ARCTPM$Control_other <- bct$t.ac12.
ARCTPM$Controlvol <- bct$t.ac13.
ARCTPM$Controlmatchedvol <- bct$t.ac14.

# Creating a dataframe "ARCTPML" (Annotations reconciled control true parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARCTPML <- data.frame(ifelse(ARCTPM[,2:5] == 1, TRUE, FALSE))

# Checking if the parameter is identical in each treatment group row of each study:
all(ARCTPML$ControlType)
all(ARCTPML$Control_other)
all(ARCTPML$Controlvol)
all(ARCTPML$Controlmatchedvol)

# Which articles have multiple control types?:
ARCTPM %>% filter(ControlType > 1) %>% select_(1)

# Adding a column to AR with the combined info of the type of control including others:
AR_updated$Control_type_combined <- ifelse(AR_updated$What.is.the.type.of.control.treatment. != "other" & AR_updated$What.is.the.type.of.control.treatment. != "media & PBS",
                                           AR_updated$What.is.the.type.of.control.treatment.,
                                           ifelse(AR_updated$What.is.the.type.of.control.treatment. == "media & PBS",
                                                  "medium without DMSO",
                                                  ifelse(!is.na(AR_updated$What.other.control.did.they.use.), "medium without DMSO", NA)))

# Adding one missing control type after looking in the original article:
AR[AR$StudyIdStr == "31268e7c-5b50-47cb-81cd-65a249bedf43", 49] <- "DMSO"
AR_updated$Control_type_combined[AR_updated$StudyIdStr == "31268e7c-5b50-47cb-81cd-65a249bedf43"] <- "DMSO"
AR[AR$StudyIdStr == "31268e7c-5b50-47cb-81cd-65a249bedf43", 57] <- FALSE
AR_updated[AR_updated$StudyIdStr == "31268e7c-5b50-47cb-81cd-65a249bedf43", 57] <- FALSE

# One row per study and control type and control volume:
ARC1R <- AR_updated %>% distinct(joined_StudyID_Control,.keep_all = TRUE) %>% select_(1,65,59,57) %>% filter(Control == TRUE)

# Table of the control types:
ARC1R %>% filter(Control == TRUE) %>% select_(2) %>% table()

# There are 12 studies with not reported control types missing:
Missing_control_type <- ARIP1R %>% filter(!is.element(StudyIdStr, ARC1R$StudyIdStr)) %>% select_(1) %>% as.data.frame()
Missing_control_type$Control_type_combined <- paste("not reported")
Missing_control_type$Control <- TRUE
Missing_control_type$Is.the.volume.of.the.added.control.suspension.reported. <- FALSE

# Merging ARC1R and Missing_control_types and using this as updated ARC1R:
ARC1R <- rbind(ARC1R, Missing_control_type)

# Ordering ARC1R by StudyID:
ARC1R <- ARC1R[order(ARC1R$StudyIdStr),]
# the ones that have no control type info are considered as "not reported". This affects the following number of studies:

# table of the control types:
ARC1R %>% filter(Control == TRUE) %>% select_(2) %>% table()

# table of the volume of added control suspension:
ARC1R %>% filter(Control == TRUE) %>% select_(4) %>% table()


## Frequency of concentration of cells plated parameter phenotypes:

# Adding a joined StudyID and ExperimentLabel column to AR:
AR_updated$joined_StudyID_EL <- paste(AR$StudyIdStr, AR$ExperimentLabel)

# Via creating a new column with the combined cell concentrations to AR

# For the data rows that have a defined concentration of cells reported:
AR_updated$combined_cell_concentration <- ifelse(AR_updated$Is.a.defined.concentration.or.a.range.of.conc..of.cells.plated.reported. == "defined concentration", AR_updated$What.is.the.defined.concentration.of.cells.plated., NA)

# Transforming the volume the cells are plated in into numerical values:
AR_updated$What.is.the.volume.the.cells.are.plated.in. <- as.numeric(AR_updated$What.is.the.volume.the.cells.are.plated.in.)

# For the data rows that have the number and the volume of the plated cells reported we can calculate the concentration:
AR_updated$combined_cell_concentration <-  ifelse(
  AR_updated$Is.the.number...concentration...volume.the.U87.cells.are.plated.reported. == "c(\"number of cells plated is reported\", \"volume the cells are plated in is reported\")" &
    AR_updated$Is.a.defined.number.or.a.range.of.cells.plated.reported. == "defined number",
  AR_updated$What.is.the.number.of.cells.plated.per.well. / AR_updated$What.is.the.volume.the.cells.are.plated.in., 
  AR_updated$combined_cell_concentration)
# data rows with ranges of numbers of cells were not considered!

# Checking if it was successful (how many studies have combined cell concentrations?):
AR_updated %>% distinct(StudyIdStr, .keep_all = TRUE) %>% select_(66) %>% table()
# seems to be successful

# Creating a df with one row for every study and cell concentration:
AR1Rcellconc <- AR_updated %>% distinct(StudyIdStr, .keep_all = TRUE) %>% select_(1,67)
AR1Rcellconc <- AR1Rcellconc[order(AR1Rcellconc$StudyIdStr),]

# Checking if the combined concentration is identical in every row of each study:

# Every result for each parameter gets storaged in a list "Listc(Number)" with the helping data "dummyc(Number)":
Listceco1 <- list()
for (ceco in seq(1, nrow(ARSC), 1)) {
  dummyceco1 <- AR_updated %>% filter(StudyIdStr == ARSC[ceco,1]) %>% select_(62) %>% unique() %>% nrow() 
  Listceco1[[length(Listceco1)+1]] = dummyceco1
}

# Every List gets converted into a dataframe "ac(Number)":
aceco1 <- data.frame(Listceco1)

# Creating a helping dataframe "bct" with the column/row reversed dataframes of "ac"
bceco <- data.frame(t(aceco1))

# Creating the dataframe ARCTPM (Annotations reconciled control true parameter match) with the column specifications for each parameter based on the "bct" dataframe
ARCeco <- data.frame(ARSC$StudyIdStr)
ARCeco$comb_cell_concentration <- bceco$t.aceco1.

# Creating a dataframe "ARCTPML" (Annotations reconciled control true parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARCeco <- data.frame(ifelse(ARCeco[,2] == 1, TRUE, FALSE))

# Checking if the parameter is identical in each treatment group row of each study:
all(ARCeco$ifelse.ARCeco...2.....1..TRUE..FALSE.)
# the combined concentration is identical in every row of each study!

# Frequency of combined concentration of cells plated parameter phenotypes:
AR_updated %>% distinct(StudyIdStr, .keep_all = TRUE) %>% select_(67) %>% table()

# Is there a study with a range of cells and volume reported?
AR_updated %>% filter(Is.a.defined.number.or.a.range.of.cells.plated.reported. == "range") %>% filter(!is.na(What.is.the.volume.the.cells.are.plated.in.)) %>% select_(1) %>% unique()
# there is one study (9f69cc7c-6213-447a-8b4f-d0b3d080596a)

# Calculating the range of concentrations for this study
AR_updated %>% filter(StudyIdStr == "9f69cc7c-6213-447a-8b4f-d0b3d080596a") %>% select_(46)
AR_updated %>% filter(StudyIdStr == "9f69cc7c-6213-447a-8b4f-d0b3d080596a") %>% select_(50)
# so the cell concentration ranges from 1000/80 cells/microl to 5000/80 cells/microl (= 12.5 - 62.5 cells/microl)
# The range is only relevant for the frequencies of phenotypes but not as a Meta-analysis moderator!




#### Reporting quality analysis

# General information: TRUE means "not reported" and FALSE means "reported"

# Creating a dataframe with one row per article using AR and not AR_updated as the basis (before contacting authors)
ARIP1R <- AR[, c(1,13,17:22,25:27,33:34,36:37,39,54)] %>% distinct(StudyIdStr, .keep_all = T)

# Orderung ARIP1R by StudyID:
ARIP1R <- ARIP1R[order(ARIP1R$StudyIdStr),]

# checking if the studyIDs are identical:
all(ARIP1R$StudyIdStr == AR1Rglucose$StudyIdStr)
all(ARIP1R$StudyIdStr == ARC1R$StudyIdStr)
all(ARIP1R$StudyIdStr == AR1Rcellconc$StudyIdStr)

# Adding some columns to ARIP1R:
ARIP1R$assay <- FALSE
ARIP1R$glucose <- AR1Rglucose$How.is.the.glucose.level.of.the.DMEM.medium.for.U87.
ARIP1R$TMZ_source <- AR1Rtmz_source$Is.the.source.of.temozolomide.reported.
ARIP1R$control_type <- ARC1R$Control_type_combined
ARIP1R$TMZ_vol <- AR1Rtmz_source$Is.the.volume.of.the.added.Temozolomide.suspension.reported.
ARIP1R$control_vol <- ARC1R$Is.the.volume.of.the.added.control.suspension.reported.
ARIP1R$cell_conc <- AR1Rcellconc$combined_cell_concentration

# Creating a helping dataframe "b" with the column/row reversed dataframes of "a"
b <- data.frame(t(a1), t(a2), t(a3), t(a4), t(a5), t(a7), t(a8), t(a9), t(a10), t(a11), t(a12), t(a13), t(a14), t(a15), t(a16), t(a17), t(a18))

# Creating the dataframe ARPM (Annotations reconciled parameter match) with the column specifications for each parameter based on the "b" dataframe
ARPM <- data.frame(ARSC$StudyIdStr)
ARPM$CoI <- b$t.a1.
ARPM$U87otherSource <- b$t.a2.
ARPM$Passages <- b$t.a3.
ARPM$Authentification <- b$t.a4.
ARPM$Mycoplasma <- b$t.a5.
ARPM$U87age <- b$t.a7.
ARPM$Confluency <- b$t.a8.
ARPM$TimeCellPassaging <- b$t.a9.
ARPM$U87source <- b$t.a10.
ARPM$FBSuse <- b$t.a11.
ARPM$otherAntibiotics <- b$t.a12.
ARPM$Antibiotics <- b$t.a13.
ARPM$FBSsource <- b$t.a14.
ARPM$FBSotherSource <- b$t.a15.
ARPM$CellPassagingCriteria <- b$t.a16.
ARPM$Glucose <- b$t.a17.
ARPM$Cell_number_vol_conc <- b$t.a18.

# Creating a dataframe "ARPML" (Annotations reconciled parameter match logical) where TRUE means that all the parameter are identical in all rows of a study
ARPML <- data.frame(ifelse(ARPM[,2:18] == 1, TRUE, FALSE))
# Adding a Column with the StudyIDs:
ARPML <- add_column(ARPML, "StudyIdStr" <- ARSC$StudyIdStr, .before = "CoI")
# Changing the column name of the studyIDs:
colnames(ARPML)[colnames(ARPML) == "\"StudyIdStr\" <- ARSC$StudyIdStr"] <- "StudyIdStr"

# Checking if the parameters are identical in all studies:
all(ARPML$CoI)
all(ARPML$U87otherSource)
all(ARPML$Passages)
all(ARPML$Authentification)
all(ARPML$Mycoplasma)
all(ARPML$U87age)
all(ARPML$Confluency)
all(ARPML$TimeCellPassaging)
all(ARPML$U87source)
all(ARPML$FBSuse)
all(ARPML$otherAntibiotics)
all(ARPML$Antibiotics)
all(ARPML$FBSsource)
all(ARPML$FBSotherSource)
all(ARPML$CellPassagingCriteria)
all(ARPML$Glucose)
all(ARPML$Cell_number_vol_conc)
# There are four parameters (mycoplasma exclusion, glucose, U87 source and other U87 source) with different phenotypes within the same study. 
# These parameters were checked again based on the original articles.


## Corrections after looking into the original articles:

# Parameter: U87otherSource
ARPML %>% filter(U87otherSource != 1) %>% select_(1)
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",17]
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",17] <- "Invitrogen Life Technologies, Inc., Carlsbad, CA"
AR[AR$StudyIdStr == "4f7e8b69-fe49-476a-bf2c-b8e84d4229ba",17]
AR[AR$StudyIdStr == "4f7e8b69-fe49-476a-bf2c-b8e84d4229ba",17] <- "China Academia Sinica cell repository (Shanghai, China)"

# Parameter: Mycoplasma
ARPML %>% filter(Mycoplasma != 1) %>% select_(1)
AR[AR$StudyIdStr == "65e0c64b-fd79-4f96-b170-ef785d6aa038",20]
AR[AR$StudyIdStr == "65e0c64b-fd79-4f96-b170-ef785d6aa038",20] <- TRUE

# Parameter: U87Source:
ARPML %>% filter(U87source != 1) %>% select_(1)
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",26]
AR[AR$StudyIdStr == "4460c6e8-45ac-4cf5-9d06-ebbdaafb7912",26] <- "other source"

# Parameter: Glucose
ARPML %>% filter(Glucose != 1) %>% select_(1)
AR[AR$StudyIdStr == "5d6b5ba3-32aa-4c8e-bcbc-0bc0b4efc606",38]
# -> this study contains experiments with both high and low glucose! -> no correction needed

# Parameter: number of animals 
# additional correction of an observed mistake in data extraction after looking again into the original article:
AR$Number.of.animals.in.cohort[AR$StudyIdStr == "f0c825a5-5083-4867-8af2-38cf451117c6" &
                                 AR$ExperimentLabel == "e2 (tb)"] <- 4


## Checking again if the parameters are identical in all studies (same code as before using loops):

# Getting the number of unique values of experiment parameters for each study with a for-loop:
# Every result for each parameter gets stored in a list "List(Number)" with the helping data "dummy(Number)":
# Note that number 6 is missing because one parameter was present twice and removed afterwards.
List1 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy1 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(13) %>% unique() %>% nrow() 
  List1[[length(List1)+1]] = dummy1
}
List2 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy2 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(17) %>% unique() %>% nrow() 
  List2[[length(List2)+1]] = dummy2
}
List3 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy3 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(18) %>% unique() %>% nrow() 
  List3[[length(List3)+1]] = dummy3
}
List4 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy4 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(19) %>% unique() %>% nrow() 
  List4[[length(List4)+1]] = dummy4
}
List5 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy5 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(20) %>% unique() %>% nrow() 
  List5[[length(List5)+1]] = dummy5
}
List7 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy7 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(21) %>% unique() %>% nrow() 
  List7[[length(List7)+1]] = dummy7
}
List8 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy8 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(22) %>% unique() %>% nrow() 
  List8[[length(List8)+1]] = dummy8
}
List9 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy9 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(25) %>% unique() %>% nrow() 
  List9[[length(List9)+1]] = dummy9
}
List10 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy10 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(26) %>% unique() %>% nrow() 
  List10[[length(List10)+1]] = dummy10
}
List11 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy11 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(33) %>% unique() %>% nrow() 
  List11[[length(List11)+1]] = dummy11
}
List12 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy12 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(34) %>% unique() %>% nrow() 
  List12[[length(List12)+1]] = dummy12
}
List13 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy13 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(36) %>% unique() %>% nrow() 
  List13[[length(List13)+1]] = dummy13
}
List14 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy14 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(37) %>% unique() %>% nrow() 
  List14[[length(List14)+1]] = dummy14
}
List15 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy15 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(39) %>% unique() %>% nrow() 
  List15[[length(List15)+1]] = dummy15
}
List16 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy16 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(27) %>% unique() %>% nrow() 
  List16[[length(List16)+1]] = dummy16
}
List17 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy17 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(38) %>% unique() %>% nrow() 
  List17[[length(List17)+1]] = dummy17
}
List18 <- list()
for (i in seq(1, nrow(ARSC), 1)) {
  dummy18 <- AR %>% filter(StudyIdStr == ARSC[i,1]) %>% select_(54) %>% unique() %>% nrow() 
  List18[[length(List18)+1]] = dummy18
}

# Every List gets converted into a dataframe "a(Number)":
a1 <- data.frame(List1)
a2 <- data.frame(List2)
a3 <- data.frame(List3)
a4 <- data.frame(List4)
a5 <- data.frame(List5)
a7 <- data.frame(List7)
a8 <- data.frame(List8)
a9 <- data.frame(List9)
a10 <- data.frame(List10)
a11 <- data.frame(List11)
a12 <- data.frame(List12)
a13 <- data.frame(List13)
a14 <- data.frame(List14)
a15 <- data.frame(List15)
a16 <- data.frame(List16)
a17 <- data.frame(List17)
a18 <- data.frame(List18)




###### Reporting quality analysis

# General information: TRUE means "not reported" and FALSE means "reported"

# Creating a dataframe with one row per article using AR and not AR_updated as the basis (before contacting authors)
ARIP1R <- AR[, c(1,13,17:22,25:27,33:34,36:37,39,54)] %>% distinct(StudyIdStr, .keep_all = T)

# Orderung ARIP1R by StudyID:
ARIP1R <- ARIP1R[order(ARIP1R$StudyIdStr),]

# checking if the studyIDs are identical:
all(ARIP1R$StudyIdStr == AR1Rglucose$StudyIdStr)
all(ARIP1R$StudyIdStr == AR1Rtmz_source$StudyIdStr)
all(ARIP1R$StudyIdStr == ARC1R$StudyIdStr)
all(ARIP1R$StudyIdStr == AR1Rcellconc$StudyIdStr)

# Adding some columns to ARIP1R:
ARIP1R$assay <- FALSE
ARIP1R$glucose <- AR1Rglucose$How.is.the.glucose.level.of.the.DMEM.medium.for.U87.
ARIP1R$TMZ_source <- AR1Rtmz_source$Is.the.source.of.temozolomide.reported.
ARIP1R$control_type <- ARC1R$Control_type_combined
ARIP1R$TMZ_vol <- AR1Rtmz_source$Is.the.volume.of.the.added.Temozolomide.suspension.reported.
ARIP1R$control_vol <- ARC1R$Is.the.volume.of.the.added.control.suspension.reported.
ARIP1R$cell_conc <- AR1Rcellconc$combined_cell_concentration


## Antibiotics:

# Some antibiotics parameters have to be changed to not reported because in the section "other antibiotics" they did not specify the antibiotic supplementation:
ARIP1R[ARIP1R$StudyIdStr == "3c41250d-a6ef-47cd-bdce-1323921aa009",14] <- "not reported"
ARIP1R[ARIP1R$StudyIdStr == "828c1097-b791-43fe-9bf8-3cc8b020bdd7",14] <- "not reported"
ARIP1R[ARIP1R$StudyIdStr == "cc7b32c6-5485-4ed8-90cf-562706fa6955",14] <- "not reported"


## Duration of treatment (time):

# It is already manually assessed in the crucial parameters reporting dataframe (all parameters assessed by two independent reviewers, reconciliation by a third independent reviewer)
Overview_about_critical_experiment_parameter_reporting <- read.csv("Overview_about_critical_experiment_parameter_reporting.csv")

# Ordering the overview of critical parameters reportng dataframe by StudyId:
OaCEPRordered <- Overview_about_critical_experiment_parameter_reporting[order(Overview_about_critical_experiment_parameter_reporting$StudyID),]

# Checking if the filtered StudyIDs are identical:
ARIP1R$StudyIdStr == OaCEPRordered$StudyID[is.element(OaCEPRordered$StudyID, ARIP1R$StudyIdStr)]

# Which studies are not an element of ARIP1R?:
OaCEPRordered %>% filter(!is.element(StudyID,ARIP1R$StudyIdStr)) %>% select_(1)

# These four studies have to be removed from OaCEPRordered:
OaCEPRordered <- OaCEPRordered %>% filter(is.element(StudyID,ARIP1R$StudyIdStr))
# Now the two later added articles are missing

# Manually adding them (they both have the treatment duration and TMZ concentration reported):
OaCEPRordered[136,1] <- "fa8c41e6-ab4d-4e97-a846-75208625a8cc"
OaCEPRordered[137,1] <- "31268e7c-5b50-47cb-81cd-65a249bedf43"

# ordering again:
OaCEPRordered <- OaCEPRordered[order(OaCEPRordered$StudyID),]
ARIP1R <- ARIP1R[order(ARIP1R$StudyIdStr),]

# checking again if the filtered StudyIDs are identical:
all(ARIP1R$StudyIdStr == OaCEPRordered$StudyID)
# all values should be TRUE!

# Adding it to ARIP1R:
ARIP1R$Treatment_duration <- OaCEPRordered$Time_of_exposure


## Temozolomide concentration reporting:
ARIP1R$TMZ_concentration <- OaCEPRordered$tmz_concentration


## Errorype reporting:
ARIP1R$errortype <- OaCEPRordered$Errortype

# Changing the parameter from an additional study to not reported (after checking again in the original research article):
ARIP1R[ARIP1R$StudyIdStr == "fa8c41e6-ab4d-4e97-a846-75208625a8cc", 26] <- "not reported"


## Experiment number:
ARIP1R$exp_number <- OaCEPRordered$Experiment_number

# are there data with n = 0 that are not in the list of OaCEPRordered$Experiment_number?
Studies_without_exp_num <- ARIP1R %>% filter(exp_number == "not reported") %>% select_(1)
Studies_without_exp_num_AR <- AR %>% filter(Number.of.animals.in.cohort == 0) %>% select_(1) %>% unique()
Studies_without_exp_num %>% filter(!is.element(StudyIdStr, Studies_without_exp_num_AR$StudyIdStr))
Studies_without_exp_num_AR %>% filter(!is.element(StudyIdStr, Studies_without_exp_num$StudyIdStr))

# Changing this two experiment number phenotypes to "not reported" (after checking again in the original research article):
ARIP1R[ARIP1R$StudyIdStr == "9a9433d1-7245-4510-9f56-8032ad295081", 28] <- "not reported"
ARIP1R[ARIP1R$StudyIdStr == "fa8c41e6-ab4d-4e97-a846-75208625a8cc", 28] <- "not reported"



### Creating the dataframe "ARRQSP" (abbr. for annotations reconciled reporting quality study parameters):
ARRQSP <- data.frame(ARIP1R$StudyIdStr,
                     ARIP1R$Is.a.successfull.mycoplasma.contamination.exclusion.of.the.cell.line.reported. == FALSE,
                     ARIP1R$Is.the.age.of.the.cell.line.reported. == FALSE,
                     ARIP1R$From.whom.did.they.get.the.U87.cells. == "not reported",
                     ARIP1R$How.ist.the.reporting.on.cell.passaging.criteria. == "not reported",
                     ARIP1R$Which.antibiotic.did.they.used.in.the.U87.cell.line.medium. == "not reported",
                     ARIP1R$If.they.used.FCS.FBS..from.which.manufactur.did.they.purchased.it. == "not reported",
                     ARIP1R$Was.a.cell.line.authentification.of.U87.cells.performed. == FALSE,
                     ARIP1R$assay,
                     ARIP1R$glucose == "not reported",
                     ARIP1R$TMZ_source == FALSE,
                     ARIP1R$TMZ_vol == FALSE,
                     ARIP1R$control_type == "not reported",
                     ARIP1R$control_vol == FALSE,
                     is.na(ARIP1R$cell_conc),
                     ARIP1R$Treatment_duration == "not reported",
                     ARIP1R$TMZ_concentration == "not reported",
                     ARIP1R$errortype == "not reported",
                     ARIP1R$exp_number == "not reported",
                     ARIP1R$How.is.the.reporting.of.conflicts.of.interest. == "not reported at all")

# Converting columns to logical values:
# this affects column 16 to 19
ARRQSP[, 16:19] <- ifelse(is.na(ARRQSP[, 16:19]), FALSE, TRUE)

# converting to numerical data (TRUE = 1, FALSE = 0):
ARRQSP[, 2:20] <- ifelse(ARRQSP[, 2:20] == TRUE, 1, 0) 
# 1 means "not reported", 0 means "reported"


## Adding the summary score and the reporting percentage:

# Viability assay as well as TMZ and Control volume get excluded from this calculations because the reported concentration can be assumed as working concentration of the whole solution (medium + drug) and reporting of cell viability assay was already an inclusion criterion and therefore reported 100%.
ARRQSP$Summary <- rowSums(ARRQSP[,c(2:8,10,11,13,15:20)])
ARRQSP$Reporting_percentage <- percent((16-ARRQSP$Summary) / 16)
ARRQSP$Reporting_ratio <- (16-ARRQSP$Summary) / 16



### Years of publication linked with reporting quality

# Adding the years of publication:
Publication_Years <- read.csv("Included Studies - Year and Correspondence.csv")[,c(1,3)]
Publication_Years_relevant <- Publication_Years %>% filter(is.element(Publication_Years$StudyID, ARRQSP$ARIP1R.StudyIdStr))

# Which articles have missing publication years?
ARRQSP %>% filter(!is.element(ARRQSP$ARIP1R.StudyIdStr, Publication_Years_relevant$StudyID)) %>% select_(1)

# Adding two missing publication years (after manually checking the publication year):
Publication_Years_relevant[136,1] <- "31268e7c-5b50-47cb-81cd-65a249bedf43"
Publication_Years_relevant[136,2] <- 2018
Publication_Years_relevant[137,1] <- "fa8c41e6-ab4d-4e97-a846-75208625a8cc"
Publication_Years_relevant[137,2] <- 2018

# Ordering both ARRQSP and Publication_Years_relevant by StudyID:
ARRQSP <- ARRQSP[order(ARRQSP$ARIP1R.StudyIdStr),]
Publication_Years_relevant <- Publication_Years_relevant[order(Publication_Years_relevant$StudyID),]

# Checking if the StudyIDs are identical:
all(ARRQSP$ARIP1R.StudyIdStr == Publication_Years_relevant$StudyID)
# if all values are TRUE we can continue

# Adding a year of publication column to the ARRQSP dataframe:
ARRQSP$Year <- Publication_Years_relevant$Year.of.publication

# table of years:
table(ARRQSP$Year)


### Calculating the mean reporting ratio for the different years (the years before 2012 are summarized because otherwise there would be to small groups)

# Storing the results in a new dataframe RQ_Years:
Years <- c("2003 - 2011", seq(2012,2020,1))

# Creating a new dataframe:
RQ_Years <- data.frame(Years)
RQ_Years$Years <- RQ_Years$Years %>% as.numeric()
RQ_Years[1,1] <- paste("2003 - 2011")
RQ_Years[1,2] <- ARRQSP %>% filter(Year < 2012) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[2,2] <- ARRQSP %>% filter(Year == 2012) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[3,2] <- ARRQSP %>% filter(Year == 2013) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[4,2] <- ARRQSP %>% filter(Year == 2014) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[5,2] <- ARRQSP %>% filter(Year == 2015) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[6,2] <- ARRQSP %>% filter(Year == 2016) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[7,2] <- ARRQSP %>% filter(Year == 2017) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[8,2] <- ARRQSP %>% filter(Year == 2018) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[9,2] <- ARRQSP %>% filter(Year == 2019) %>% select_(23) %>% unlist() %>% mean()
RQ_Years[10,2] <- ARRQSP %>% filter(Year == 2020) %>% select_(23) %>% unlist() %>% mean()

# changing a column name:
colnames(RQ_Years)[2] <- paste("Mean reporting ratio")

# Calculating the standard deviation of the mean reporting ratio for the different years:
RQ_Years[1,3] <- ARRQSP %>% filter(Year < 2012) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[2,3] <- ARRQSP %>% filter(Year == 2012) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[3,3] <- ARRQSP %>% filter(Year == 2013) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[4,3] <- ARRQSP %>% filter(Year == 2014) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[5,3] <- ARRQSP %>% filter(Year == 2015) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[6,3] <- ARRQSP %>% filter(Year == 2016) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[7,3] <- ARRQSP %>% filter(Year == 2017) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[8,3] <- ARRQSP %>% filter(Year == 2018) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[9,3] <- ARRQSP %>% filter(Year == 2019) %>% select_(23) %>% unlist() %>% sd()
RQ_Years[10,3] <- ARRQSP %>% filter(Year == 2020) %>% select_(23) %>% unlist() %>% sd()

# changing a column name:
colnames(RQ_Years)[3] <- paste("Standard deviation")

# Adding the number of studies per year:
RQ_Years[1,4] <- ARRQSP %>% filter(Year < 2012) %>% nrow()
RQ_Years[2,4] <- ARRQSP %>% filter(Year == 2012) %>% nrow()
RQ_Years[3,4] <- ARRQSP %>% filter(Year == 2013) %>% nrow()
RQ_Years[4,4] <- ARRQSP %>% filter(Year == 2014) %>% nrow()
RQ_Years[5,4] <- ARRQSP %>% filter(Year == 2015) %>% nrow()
RQ_Years[6,4] <- ARRQSP %>% filter(Year == 2016) %>% nrow()
RQ_Years[7,4] <- ARRQSP %>% filter(Year == 2017) %>% nrow()
RQ_Years[8,4] <- ARRQSP %>% filter(Year == 2018) %>% nrow()
RQ_Years[9,4] <- ARRQSP %>% filter(Year == 2019) %>% nrow()
RQ_Years[10,4] <- ARRQSP %>% filter(Year == 2020) %>% nrow()

# changing a column name:
colnames(RQ_Years)[4] <- paste("Number of studies") 



### Addition of Journal impact factors (JIF) of the journals the articles were published in in the respective year.

# The JIF factors were obtained through Clavivates.

# Load the journal impact factors:
JIF <- read.csv("SERRANO - journal impact factors.csv")
# CAVE: For the study "9cb3f84c-86a8-42e7-9237-77684a4af6ec" no Impact factor could be obtained!

# Filtering the journal impact factors that are an element of ARRQSP:
JIF <- JIF %>% filter(is.element(StudyID, ARRQSP$ARIP1R.StudyIdStr))

# Ordering by StudyID
JIF <- JIF[order(JIF$StudyID),]

# Ordering ARRQSP by StudyID:
ARRQSP <- ARRQSP[order(ARRQSP$ARIP1R.StudyIdStr),]

# Checking if StudyIDs of JIF and ARRQSP are equal:
all(JIF$StudyID == ARRQSP$ARIP1R.StudyIdStr)
# if all values are TRUE we can continue!

# Adding the journal impact factors to ARRQSP:
ARRQSP$JIF <- JIF$Impact.factor




#### Analysis of the correlation between articles reporting quality and the year they were published in as well as the JIF

### Correlation between the articles reporting quality and year of publication:

# Linear regression model with the year as independent variable:
model_RQ_years <- lm((ARRQSP$Reporting_ratio *100) ~ ARRQSP$Year)
summary(model_RQ_years)



### Correlation between the articles reporting quality and JIF:

# Linear regression model with the JIF as independent variable:
model_RQ_JIF <- lm(Reporting_ratio~JIF, data = ARRQSP, na.action = na.omit)
summary(model_RQ_JIF)





######## Meta-analysis


# Correcting some mistakes in data extraction after manually checking extreme values again in the original research articles:
AR_OD_MA_T_f$Error[AR_OD_MA_T_f$StudyIdStr == "e97857cd-aa4e-4af4-9448-6f6973fd6755" & AR_OD_MA_T_f$Time == 3 & AR_OD_MA_T_f$tmz_concentration_in_microM_numeric == 50] <- 23.5
AR_OD_MA_T_f$SD[AR_OD_MA_T_f$StudyIdStr == "e97857cd-aa4e-4af4-9448-6f6973fd6755" & AR_OD_MA_T_f$Time == 3 & AR_OD_MA_T_f$tmz_concentration_in_microM_numeric == 50] <- 23.5
AR_OD_MA_T_f$SD.1[AR_OD_MA_T_f$StudyIdStr == "65b98083-98fc-469f-bd94-4a2074894498" & AR_OD_MA_T_f$Time == 3] <- 2.191
AR_OD_MA_T_f$SD.1[AR_OD_MA_T_f$StudyIdStr == "35c0a928-2ad6-4912-919e-46eab215fafb"] <- 2



#### Calculation of the effect size and its variance


### Effect size D

# Raw mean difference is defined as the cell viability reduction caused by Temozolomide in comparison to the untreated control (in percent of the untreated control viability)

# Creating a new dataframe "AR_OD_MA_T_f raw mean difference":
AR_OD_MA_T_f_RMD <- AR_OD_MA_T_f

# Correcting some  mistakes in data extraction after manually checking extreme value outliers based on the original research articles data:
AR_OD_MA_T_f_RMD$Error[AR_OD_MA_T_f_RMD$StudyIdStr == "e97857cd-aa4e-4af4-9448-6f6973fd6755" & AR_OD_MA_T_f_RMD$Time == 3 & AR_OD_MA_T_f_RMD$tmz_concentration_in_microM_numeric == 50] <- 23.5
AR_OD_MA_T_f_RMD$SD[AR_OD_MA_T_f_RMD$StudyIdStr == "e97857cd-aa4e-4af4-9448-6f6973fd6755" & AR_OD_MA_T_f_RMD$Time == 3 & AR_OD_MA_T_f_RMD$tmz_concentration_in_microM_numeric == 50] <- 23.5
AR_OD_MA_T_f_RMD$SD.1[AR_OD_MA_T_f_RMD$StudyIdStr == "65b98083-98fc-469f-bd94-4a2074894498" & AR_OD_MA_T_f_RMD$Time == 3] <- 2.191
AR_OD_MA_T_f_RMD$SD.1[AR_OD_MA_T_f_RMD$StudyIdStr == "35c0a928-2ad6-4912-919e-46eab215fafb"] <- 2


### Transform all cell viability means and errors to a percentage of the control viability mean:

# Change all viability means and errors to presentation as a percentage of the control cell viability mean:
AR_OD_MA_T_f_RMD$mean_tmz <- AR_OD_MA_T_f_RMD$Average / AR_OD_MA_T_f_RMD$Average.1
AR_OD_MA_T_f_RMD$mean_control <- AR_OD_MA_T_f_RMD$Average.1 / AR_OD_MA_T_f_RMD$Average.1
AR_OD_MA_T_f_RMD$sd_tmz <- AR_OD_MA_T_f_RMD$SD / AR_OD_MA_T_f_RMD$Average.1
AR_OD_MA_T_f_RMD$sd_control <- AR_OD_MA_T_f_RMD$SD.1 / AR_OD_MA_T_f_RMD$Average.1

# Adjustments for Control viability = 0:
AR_OD_MA_T_f_RMD$mean_control <- ifelse(AR_OD_MA_T_f_RMD$Average.1 == 0, 1, AR_OD_MA_T_f_RMD$mean_control)
AR_OD_MA_T_f_RMD$mean_tmz <- ifelse(AR_OD_MA_T_f_RMD$Average.1 == 0, AR_OD_MA_T_f_RMD$Average / 100, AR_OD_MA_T_f_RMD$mean_tmz)
AR_OD_MA_T_f_RMD$sd_tmz <- ifelse(AR_OD_MA_T_f_RMD$Average.1 == 0, AR_OD_MA_T_f_RMD$SD / 100, AR_OD_MA_T_f_RMD$sd_tmz)
AR_OD_MA_T_f_RMD$sd_control <- ifelse(AR_OD_MA_T_f_RMD$Average.1 == 0, AR_OD_MA_T_f_RMD$SD.1 / 100, AR_OD_MA_T_f_RMD$sd_control)

# Adjustments for greater is worse = FALSE (data presented as growth inhibition ratio):
AR_OD_MA_T_f_RMD$mean_control <- ifelse(AR_OD_MA_T_f_RMD$Greater.is.worse. == FALSE, 1, AR_OD_MA_T_f_RMD$mean_control)
AR_OD_MA_T_f_RMD$mean_tmz <- ifelse(AR_OD_MA_T_f_RMD$Greater.is.worse. == FALSE, 1 - AR_OD_MA_T_f_RMD$Average/(100 - AR_OD_MA_T_f_RMD$Average.1), AR_OD_MA_T_f_RMD$mean_tmz)
AR_OD_MA_T_f_RMD$sd_tmz <- ifelse(AR_OD_MA_T_f_RMD$Greater.is.worse. == FALSE, AR_OD_MA_T_f_RMD$SD/(100 - AR_OD_MA_T_f_RMD$Average.1), AR_OD_MA_T_f_RMD$sd_tmz)
AR_OD_MA_T_f_RMD$sd_control <- ifelse(AR_OD_MA_T_f_RMD$Greater.is.worse. == FALSE, AR_OD_MA_T_f_RMD$SD.1/(100 - AR_OD_MA_T_f_RMD$Average.1), AR_OD_MA_T_f_RMD$sd_control)

# Adjustment for some infinite sd_control values (where growth inhibition ratios were presented):
AR_OD_MA_T_f_RMD$sd_control <- ifelse(AR_OD_MA_T_f_RMD$sd_control == Inf, AR_OD_MA_T_f_RMD$SD.1/100, AR_OD_MA_T_f_RMD$sd_control) 


### calculating raw mean difference effect sizes and their variances with the escalc function from metafor while creating a new dataframe:
escalcnew <- escalc(measure = "MD",n1i = AR_OD_MA_T_f_RMD$Number.of.animals.in.cohort.1, n2i = AR_OD_MA_T_f_RMD$Number.of.animals.in.cohort, m1i = AR_OD_MA_T_f_RMD$mean_control, m2i = AR_OD_MA_T_f_RMD$mean_tmz, sd1i = AR_OD_MA_T_f_RMD$sd_control, sd2i = AR_OD_MA_T_f_RMD$sd_tmz, data = AR_OD_MA_T_f_RMD)



# Excluding data points with treatment time = 0 (as no effect of the intervention can be expected at this baseline time):
escalcnew <- escalcnew %>% filter(Time != 0)




#### Addition of possible moderators to the outcome values dataframe

# Using the "rn" column that shows the row number of the data in the AR dataframe:

### Moderators from AR:
escalcnew$authentification <- AR_updated[escalcnew$rn, 19]
escalcnew$mycoplasma <- AR_updated[escalcnew$rn, 20]
escalcnew$glucose <- AR_updated[escalcnew$rn, 38]



### Combined antibiotics:

# changing the antibiotics from "other" to "not reported" because in these articles they did not name the antibiotics used (manually checked):
AR_updated$Which.antibiotic.did.they.used.in.the.U87.cell.line.medium.[AR_updated$StudyIdStr == "3c41250d-a6ef-47cd-bdce-1323921aa009" |
                                                                         AR_updated$StudyIdStr == "828c1097-b791-43fe-9bf8-3cc8b020bdd7" |
                                                                         AR_updated$StudyIdStr == "cc7b32c6-5485-4ed8-90cf-562706fa6955"] <- "not reported"

# Adding it to AR_OD_MA_T_f_RMD
escalcnew$antibiotics <- AR_updated[escalcnew$rn,36]



### U87 age (passages):

# Adding it to AR_OD_MA_T_f_RMD
escalcnew$U87_age_passages <- AR_updated[escalcnew$rn, 18]

# replace empty cells with NA:
escalcnew$U87_age_passages[escalcnew$U87_age_passages == ""] <- NA

# Assuming the center of a range if the age was reported as range:
escalcnew$U87_age_passages[escalcnew$joined_StudyID_EL == "65e0c64b-fd79-4f96-b170-ef785d6aa038 e"] <- 10

# converting the age into numeric values:
escalcnew$U87_age_passages <- as.numeric(escalcnew$U87_age_passages)



### Type of control:

# Loading a manually written dataframe with the control parameters of the articles (also extracted by two independent reviewers, reconciliation by a third independent reviewer):
ARSCT1R <- read.csv("ARSCT1R after adding and removing studies on 25th february.csv") 

# Adding parameters for a study that had missing control parameters:
ARSCT1R[131,1] <- paste("f08c5376-874b-4753-b626-98f3d9870491 TRUE")
ARSCT1R[131,2] <- TRUE
ARSCT1R[131,3] <- paste("not reported")
ARSCT1R[131,5] <- FALSE
ARSCT1R[131,7] <- paste("f08c5376-874b-4753-b626-98f3d9870491")

# checking if all studies in the meta-analysis are an element of ARSCT1R:
escalcnew %>% filter(is.element(StudyIdStr, ARSCT1R$StudyIdStr)) %>% nrow() == nrow(escalcnew)

# Which studies have no control parameters?:
escalcnew %>% filter(!is.element(StudyIdStr, ARSCT1R$StudyIdStr)) %>% select_(1) %>% unique()

# Adding the control parameters for this study to the ARSCT1R dataframe after looking again into the original article
ARSCT1R[132,7] <- "31268e7c-5b50-47cb-81cd-65a249bedf43"
ARSCT1R[132,3] <- "DMSO"

# again checking if all studies in the meta-analysis are an element of ARSCT1R:
escalcnew %>% filter(is.element(StudyIdStr, ARSCT1R$StudyIdStr)) %>% nrow() == nrow(escalcnew)
# if it is TRUE we can continue!

# Adding a column to ARSCT1R with the combined info of the type of control including others:
ARSCT1R$Control_type_combined <- ifelse(ARSCT1R$What.is.the.type.of.control.treatment. != "other" & ARSCT1R$What.is.the.type.of.control.treatment. != "media & PBS",
                                        ARSCT1R$What.is.the.type.of.control.treatment.,
                                        ifelse(ARSCT1R$What.is.the.type.of.control.treatment. == "media & PBS",
                                               "medium without DMSO",
                                               ifelse(!is.na(ARSCT1R$What.other.control.did.they.use.), "medium without DMSO", NA)))

# Changing control types that were obtained through contacting authors:
ARSCT1R$Control_type_combined[ARSCT1R$StudyIdStr == "64a785ad-4714-4ec6-b731-4d66128602f7" | ARSCT1R$StudyIdStr == "4b98515f-057a-4d17-af4a-527e9b58ec36"] <- "medium without DMSO"
ARSCT1R$Control_type_combined[ARSCT1R$StudyIdStr == "86d1eac3-4f54-4d81-bd3d-673f115de939" |
                                ARSCT1R$StudyIdStr == "e8c88cd2-830a-4b62-95e8-18fdeaa25347" |
                                ARSCT1R$StudyIdStr == "6f3a6613-897c-4bb6-9f5a-a510bf1ec6ad" |
                                ARSCT1R$StudyIdStr == "bb056e51-f833-4bbc-a8f1-eecfe62bf9c9"] <- "DMSO"

# Adding the type of control using the escalcnew dataframe:
for (ctparameter in 1:nrow(escalcnew)) {
  escalcnew[ctparameter, 32] <- ARSCT1R %>% filter(StudyIdStr == escalcnew[ctparameter, 1]) %>% select_(8) 
}



### Concentration of cells plated:

# Adding a joined StudyID and ExperimentLabel column to AR:
AR_updated$joined_StudyID_EL <- paste(AR$StudyIdStr, AR$ExperimentLabel)


## Creating a new column with the combined cell concentrations to AR:

# For the data rows that have a defined concentration of cells reported:
AR_updated$combined_cell_concentration <- ifelse(AR_updated$Is.a.defined.concentration.or.a.range.of.conc..of.cells.plated.reported. == "defined concentration", AR_updated$What.is.the.defined.concentration.of.cells.plated., NA)

# Transforming the volume the cells are plated in into numerical values:
AR_updated$What.is.the.volume.the.cells.are.plated.in. <- as.numeric(AR_updated$What.is.the.volume.the.cells.are.plated.in.)

# For the data rows that have the number and the volume of the plated cells reported we can calculate the concentration:
AR_updated$combined_cell_concentration <-  ifelse(
  AR_updated$Is.the.number...concentration...volume.the.U87.cells.are.plated.reported. == "c(\"number of cells plated is reported\", \"volume the cells are plated in is reported\")" &
    AR_updated$Is.a.defined.number.or.a.range.of.cells.plated.reported. == "defined number",
  AR_updated$What.is.the.number.of.cells.plated.per.well. / AR_updated$What.is.the.volume.the.cells.are.plated.in., 
  AR_updated$combined_cell_concentration)
# data rows with ranges of numbers of cells are not considered in the meta-analysis!
# manually checked, it was successful!

# adding it to AR_OD_MA_T_f_RMD:
escalcnew$cell_concentration <- AR_updated[escalcnew$rn, 69]

# converting it to numerical values:
escalcnew$cell_concentration <- as.numeric(escalcnew$cell_concentration)



### Decadic logarithm of the TMZ conc.:
escalcnew$tmz_conc_microM_log10 <- log(escalcnew$tmz_concentration_in_microM_numeric, 10)


## Decadic logarithm of the Time (treatment duration) to AR_OD_MA_T_f_RMD:
escalcnew$time_log10 <- log(escalcnew$Time, 10)



### Log-logistic TMZ concentration: 

# For the optimal fit dose-response-fitting with the Hill-equation (4 parameters) was applied:

# Building a 4-parameter hill-equation with the help of the drc package:
library(drc)
model_hill <- drm(yi ~ tmz_concentration_in_microM_numeric, fct = LL.4(), data = escalcnew)
summary(model_hill)
# Some coefficients are not significant (p > 0). However, we use this model as an possible moderator of the effect to test its fit.

# adding a column with the 4-parameter-log-logistic tmz concentration predictor:
escalcnew$loglog4pred <- model_hill$coefficients[3] + ((model_hill$coefficients[2] - model_hill$coefficients[3]) / (1 + (escalcnew$tmz_concentration_in_microM_numeric / model_hill$coefficients[4])^model_hill$coefficients[1]))



### Log-logistic Time (treatment duration):

# For the optimal fit time-response-fitting with the Hill-equation (4 parameters) was applied:

# building a 4-parameter hill-equation with the help of the drc package:
model_hill_time <- drm(yi ~ Time, fct = LL.4(), data = escalcnew)
summary(model_hill_time)

# adding a column with the 4-parameter-log-logistic time predictor:
escalcnew$loglog4predtime <- model_hill_time$coefficients[3] + ((model_hill_time$coefficients[2] - model_hill_time$coefficients[3]) / (1 + (escalcnew$Time / model_hill_time$coefficients[4])^model_hill_time$coefficients[1]))



### Number of experiments:
escalcnew$number_of_experiments <- AR[escalcnew$rn, 7]


### Articles reporting quality:
for (ctparameter3 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter3, 39] <- ARRQSP %>% filter(ARIP1R.StudyIdStr == escalcnew[ctparameter3, 1]) %>% select_(23)
}


### Year of publication:
for (ctparameter4 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter4, 40] <- ARRQSP %>% filter(ARIP1R.StudyIdStr == escalcnew[ctparameter4, 1]) %>% select_(24)
}


### Journal Impact Factor:
for (ctparameter5 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter5, 41] <- ARRQSP %>% filter(ARIP1R.StudyIdStr == escalcnew[ctparameter5, 1]) %>% select_(25)
}


### Combined U87 Source:

# As the combined U87 sources were already assessed in AR_updated we can use this dataframe
for (ctparameter6 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter6, 42] <- AR_updated %>% filter(StudyIdStr == escalcnew[ctparameter6, 1]) %>% select_(61) %>% unique()
}


### Combined cell viability assays:

# As the combined cell viability assays were already assessed in AR_updated we can use this dataframe
for (ctparameter7 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter7, 43] <- AR_updated %>% filter(StudyIdStr == escalcnew[ctparameter7, 1]) %>% select_(63) %>% unique()
}


### Combined FBS sources:

# As the combined FBS sources were already assessed in AR_updated we can use this dataframe
for (ctparameter8 in 1:nrow(escalcnew)) {
  escalcnew[ctparameter8, 44] <- AR_updated %>% filter(StudyIdStr == escalcnew[ctparameter8, 1]) %>% select_(62) %>% unique()
}


### Detailed cell passaging criteria:

# Create a dataframe with one row per article:
AR_updated_distinct <- AR_updated %>% distinct(StudyIdStr, .keep_all = T)

# overview about the cell passaging criteria
table(AR_updated_distinct$How.ist.the.reporting.on.cell.passaging.criteria.)

# create a dataframe with the detailed cell passaging criteria:
cell_passaging_criteria <- AR_updated_distinct[ , c(1,27,22,25)]
cell_passaging_criteria <- cell_passaging_criteria %>% filter(How.ist.the.reporting.on.cell.passaging.criteria. != "not reported")
cell_passaging_criteria <- cell_passaging_criteria[order(cell_passaging_criteria$How.ist.the.reporting.on.cell.passaging.criteria.), ]

# Writing a csv table for manual adjustments:
# write.table(cell_passaging_criteria, file = "cell_passaging_criteria.csv", row.names = F)

# Loading the csv with the manually reconciled numerical values:
cell_passaging_criteria <- read.csv("cell_passaging_criteria_finished.csv")
# If ranges were reported the mean got assumed for calculations:

# adding it to the AR_updated dataframe:
for(cpc in 1:nrow(AR_updated)) {
  AR_updated[cpc,68] <- ifelse(!is.element(AR_updated[cpc, 1], cell_passaging_criteria$StudyIdStr),
                               NA,
                               ifelse(
                                 !is.na(cell_passaging_criteria$confluence.numerical[cell_passaging_criteria$StudyIdStr == AR_updated[cpc,1]]),
                                 cell_passaging_criteria$confluence.numerical[cell_passaging_criteria$StudyIdStr == AR_updated[cpc,1]],
                                 NA
                               ))
}

for(cpc in 1:nrow(AR_updated)) {
  AR_updated[cpc,67] <- ifelse(!is.element(AR_updated[cpc, 1], cell_passaging_criteria$StudyIdStr),
                               NA,
                               ifelse(
                                 !is.na(cell_passaging_criteria$time.numerical[cell_passaging_criteria$StudyIdStr == AR_updated[cpc,1]]),
                                 cell_passaging_criteria$time.numerical[cell_passaging_criteria$StudyIdStr == AR_updated[cpc,1]],
                                 NA
                               ))
}

# Changing column names:
colnames(AR_updated)[67] <- paste("cpc_time_num")
colnames(AR_updated)[68] <- paste("cpc_confluence_num")

# Adding the numerical confluence based cell passaging criteria to AR_OD_MA_T_f_OD (time based criteria not because they are to few articles):
escalcnew$cpc_confluence_num <- AR_updated[escalcnew$rn, 68]

# Converting it into numerical values:
escalcnew$cpc_confluence_num <- as.numeric(escalcnew$cpc_confluence_num)





####  How many different experiments are included in the systematic review and meta-analysis, respectively?

### Number of experiments in systematic review:
AR_OD %>% distinct(StudyIdStr, ExperimentLabel) %>% nrow()

### Number of experiments in meta-analysis:
AR_OD_MA_T_f_RMD %>% distinct(StudyIdStr, ExperimentLabel) %>% nrow()

### Number of drug concentrations in systematic review:
AR_updated %>% filter(tmz_concentration_in_microM_numeric > 0) %>% nrow()

### Number of unique drug concentrations in systematic review:
AR_updated %>% filter(tmz_concentration_in_microM_numeric > 0) %>% select_(58) %>% unique() %>% nrow()

### Number of treatment durations in systematic review:
AR_updated %>% filter(tmz_concentration_in_microM_numeric > 0) %>% filter(Time > 0) %>% nrow()

### Number of unique treatment durations in systematic review:
AR_updated %>% filter(tmz_concentration_in_microM_numeric > 0) %>% filter(Time > 0) %>% select_(9) %>% unique() %>% nrow()




##### Three-level meta-analysis

# Random-effects-model
# Tau-estimator: REML (resticted maximum likelihood)
# Robust-variance-estimation

# Level 3: Each particular article
# Level 2: Every effect size
# Level 1: Raw data of every independent experiment the effect were calculated with (not available as the authors of the included articles mainly did not published these; but they are not relevant for our meta-analytical objectives)




#### Meta-analysis without moderators:

mlmodel_escalc <- rma.mv(yi, vi, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist = TRUE, test = "t")
summary(mlmodel_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation) :
coef_test(mlmodel_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(mlmodel_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# Prediction interval for the effect:
coef(mlmodel_escalc) - 1.96 * sqrt(mlmodel_escalc$sigma2[2])
coef(mlmodel_escalc) + 1.96 * sqrt(mlmodel_escalc$sigma2[2])

# I^2 for level two and three:
i2_ml(mlmodel_escalc)

# confidence intervals for tau^2:
confint.rma.mv(mlmodel_escalc)




#### Univariable meta-analysis with moderators:

### Moderator: glucose level
Glucose_escalc <- rma.mv(yi, vi, mods = ~ factor(glucose), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Glucose_escalc) 

# Confidence interval for the effects estimate (respecting robust-variance-estimation) :
coef_test(Glucose_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(Glucose_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# Prediction interval for the effect:
coef(Glucose_escalc) - 1.96 * sqrt(Glucose_escalc$sigma2[2])
coef(Glucose_escalc) + 1.96 * sqrt(Glucose_escalc$sigma2[2])

# I^2 for level two and three:
i2_ml(Glucose_escalc)

# confidence intervals for tau^2:
confint.rma.mv(Glucose_escalc)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(Glucose_escalc) 

# Explained between-articles-variance:
1 - (Glucose_escalc$sigma2[2] / mlmodel_escalc$sigma2[2])

# Number of effects per phenotype:
table(escalcnew$glucose)

# Number of articles per phenotype:
escalcnew %>% distinct(StudyIdStr,.keep_all = TRUE) %>% select_(29) %>% table()



### Moderator: U-87 MG source

U87_source_escalc <- rma.mv(yi,vi, mods = ~ factor(combined_U87_source), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(U87_source_escalc) 
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: U-87 MG authentication
U87_authentification_escalc <- rma.mv(yi, vi, mods = ~ factor(authentification), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(U87_authentification_escalc) 
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: U-87 MG age
U87_age_escalc <- rma.mv(yi, vi, mods = ~ U87_age_passages, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(U87_age_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Cell concentration
Cell_conc_escalc <- rma.mv(yi, vi, mods = ~ cell_concentration, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Cell_conc_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Confluency of cell culture at the time of cell passaging
CPC_escalc <- rma.mv(yi, vi, mods = ~ cpc_confluence_num, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(CPC_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Mycoplasma exclusion
Mycoplasma_escalc <- rma.mv(yi, vi, mods = ~ factor(mycoplasma), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Mycoplasma_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Supplemented antibiotics:
Antibiotics_escalc <- rma.mv(yi, vi, mods = ~ factor(antibiotics), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Antibiotics_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: FBS source
FBS_source_escalc <- rma.mv(yi, vi, mods = ~ factor(combined_FBS_source), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(FBS_source_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Type of untreated control
Control_Type_escalc <- rma.mv(yi,vi, mods = ~ factor(Control_type_combined), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Control_Type_escalc)
# As this parameter shows no significant result no further parameters are calculated.



### Moderator: Articles reporting quality:
Rep_ratio_escalc <- rma.mv(yi, vi, mods = ~ Reporting_ratio, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(Rep_ratio_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation) :
coef_test(Rep_ratio_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(Rep_ratio_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(Rep_ratio_escalc) 

# I^2
i2_ml(Rep_ratio_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(Rep_ratio_escalc)

# Explained between-articles-variance:
1 - (Rep_ratio_escalc$sigma2[2] / mlmodel_escalc$sigma2[2])




### Moderator: TMZ concentration (4-parameter log-logistic)
LogLog_TMZ_conc_escalc <- rma.mv(yi, vi, mods = ~ loglog4pred, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_TMZ_conc_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation) :
coef_test(LogLog_TMZ_conc_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_TMZ_conc_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_TMZ_conc_escalc) 

# I^2
i2_ml(LogLog_TMZ_conc_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_TMZ_conc_escalc)

# Explained within-articles-variance:
1 - (LogLog_TMZ_conc_escalc$sigma2[1] / mlmodel_escalc$sigma2[1])




### Moderator: Treatment duration (4-parameter log-logistic)
LogLog_Time_escalc <- rma.mv(yi, vi, mods = ~ loglog4predtime, random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_Time_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation) :
coef_test(LogLog_Time_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_Time_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_Time_escalc) 

# I^2
i2_ml(LogLog_Time_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_Time_escalc)

# Explained within-articles-variance:
1 - (LogLog_Time_escalc$sigma2[1] / mlmodel_escalc$sigma2[1])




#### Multivariable meta-analysis:

### Moderators: TMZ concentration, Treatment duration (both 4-parameter log-logistic)
LogLog_TMZ_conc_Time_escalc <- rma.mv(yi, vi, mods = ~ cbind(loglog4pred,loglog4predtime), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_TMZ_conc_Time_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation):
coef_test(LogLog_TMZ_conc_Time_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_TMZ_conc_Time_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_TMZ_conc_Time_escalc) 

# I^2
i2_ml(LogLog_TMZ_conc_Time_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_TMZ_conc_Time_escalc)



### Moderators: TMZ concentration, Treatment duration, Glucose level (first two are 4-parameter log-logistic)
LogLog_TMZ_conc_Time_Glucose_escalc <- rma.mv(yi, vi, mods = ~ cbind(loglog4pred,loglog4predtime,factor(glucose)), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_TMZ_conc_Time_Glucose_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation):
coef_test(LogLog_TMZ_conc_Time_Glucose_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_TMZ_conc_Time_Glucose_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_TMZ_conc_Time_Glucose_escalc) 

# I^2
i2_ml(LogLog_TMZ_conc_Time_Glucose_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_TMZ_conc_Time_Glucose_escalc)



### Moderators: TMZ concentration, Treatment duration, Articles reporting quality (first two are 4-parameter log-logistic)
LogLog_TMZ_conc_Time_RQ_escalc <- rma.mv(yi,vi, mods = ~ cbind(loglog4pred,loglog4predtime,Reporting_ratio), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_TMZ_conc_Time_RQ_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation):
coef_test(LogLog_TMZ_conc_Time_RQ_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_TMZ_conc_Time_RQ_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_TMZ_conc_Time_RQ_escalc) 

# I^2
i2_ml(LogLog_TMZ_conc_Time_RQ_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_TMZ_conc_Time_RQ_escalc)



### Moderators: TMZ concentration, Treatment duration, Glucose level, Articles reporting quality (first two are 4-parameter log-logistic)
LogLog_TMZ_conc_Time_Glucose_RQ_escalc <- rma.mv(yi,vi, mods = ~ cbind(loglog4pred,loglog4predtime,factor(glucose),Reporting_ratio), random = list(~ 1 | rn, ~ 1 | StudyIdStr), data = escalcnew, method = "REML", tdist=TRUE, test = "t")
summary(LogLog_TMZ_conc_Time_Glucose_RQ_escalc)

# Confidence interval for the effects estimate (respecting robust-variance-estimation):
coef_test(LogLog_TMZ_conc_Time_Glucose_RQ_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)
conf_int(LogLog_TMZ_conc_Time_Glucose_RQ_escalc, vcov = "CR2", cluster = escalcnew$StudyIdStr)

# R2 (using the method described by Nakagawa and Schielzeth, 2012 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x))
r2_ml(LogLog_TMZ_conc_Time_Glucose_RQ_escalc) 

# I^2
i2_ml(LogLog_TMZ_conc_Time_Glucose_RQ_escalc)

# Confidence intervals for tau^2 etc.
confint.rma.mv(LogLog_TMZ_conc_Time_Glucose_RQ_escalc)






###### Risk of Bias analysis:

### Preparations:

# Loading the Risk of Bias (ROB) project StudyIDs (which are different SudyIDs in comparison to the original Syrf project StudyIds):
ROB_StudyIDs <- read.csv("Risk of Bias Assessment - Systematic Review and Meta-Analysis of reporting on quality control standards in in-vitro glioblastoma research_screeningsLearn.csv")[,1:2]

# Excluding later excluded articles from the ROB_studylist:
ROB_StudyIDs <- ROB_StudyIDs %>% filter(Title != "EGFR Inhibition in Glioma Cells Modulates Rho Signaling to Inhibit Cell Motility and Invasion and Cooperates with Temozolomide to Reduce Cell Growth" &
                                          Title != "Combination of lentivirus-mediated silencing of PPM1D and temozolomide chemotherapy eradicates malignant glioma through cell apoptosis and cell cycle arrest" &
                                          Title != "The Fanconi anemia (FA) pathway confers glioma resistance to DNA alkylating agents" &
                                          PublicationNumber != "552afe89-a702-4001-8e57-1936c83048d3")

# Loading the Annotation project StudyIDs (which are the "correct" StudyIDs):
Annotation_StudyIDs <- included_articles[,1:2]

# Ordering both by the title:
ROB_StudyIDs <- ROB_StudyIDs[order(ROB_StudyIDs$Title),]
Annotation_StudyIDs <- Annotation_StudyIDs[order(Annotation_StudyIDs$Title),]

# Checking if the titles of both are identical:
ROB_StudyIDs$Title == Annotation_StudyIDs$Title
# They are not identical

# It seems like the Study "Beta-Elemene inhibits proliferation through crosstalk between glia maturation factor ? and extracellular signal-regulated kinase 1/2 and impairs drug resistance to temozolomide in glioblastoma cells" causes problems with the title.
# Fixing this:
# Finding the Study with the key unique word "crosstalk" from the title of the studies and redefining the correct title:
ROB_StudyIDs$Title[
  str_detect(ROB_StudyIDs$Title, "crosstalk") == TRUE
] <- "Beta-Elemene inhibits proliferation through crosstalk between glia maturation factor ? and extracellular signal-regulated kinase 1/2 and impairs drug resistance to temozolomide in glioblastoma cells"
Annotation_StudyIDs$Title[
  str_detect(Annotation_StudyIDs$Title, "crosstalk") == TRUE
] <- "Beta-Elemene inhibits proliferation through crosstalk between glia maturation factor ? and extracellular signal-regulated kinase 1/2 and impairs drug resistance to temozolomide in glioblastoma cells"

# Ordering both by the title again:
ROB_StudyIDs <- ROB_StudyIDs[order(ROB_StudyIDs$Title),]
Annotation_StudyIDs <- Annotation_StudyIDs[order(Annotation_StudyIDs$Title),]

# Checking if the titles of both are identical after fixing the "beta-Elemene inhibits..." study title:
all(ROB_StudyIDs$Title == Annotation_StudyIDs$Title)
# As all values are TRUE we can continue!

# Creating a dataframe with each study title and both StudyIDs:
Both_StudyIDs <- merge(Annotation_StudyIDs, ROB_StudyIDs, by = "Title")

# Changing the column names:
colnames(Both_StudyIDs) <- c("Title", "StudyID Annotation", "StudyID ROB")

# Ordering by title again:
Both_StudyIDs <- Both_StudyIDs[order(Both_StudyIDs$Title),]

# Final dataframe with both StudyIDs for each study title:
# write.table(Both_StudyIDs, file = "StudyID Comparison.csv", row.names = FALSE)


## Changing the StudyIDs of the ROB assessment results:

# Loading the Risks of Bias assessment results:
ROB <- read.csv("Risk of Bias Assessment - Systematic Review and Meta-Analysis of reporting on quality control standards in in-vitro glioblastoma research_studyAnnotationsSpread.csv")

# Excluding not included articles' ROB results:
ROB <- ROB %>% filter(is.element(StudyId, Both_StudyIDs$`StudyID ROB`))


## Checking whether the reconciliation was only conducted if differences between Timo's and Joly's results were existing:
ROB_r_new <- data.frame()
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,1] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(5) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(5), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(5) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(5), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,2] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(8) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(8), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(8) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(8), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,3] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(11) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(11), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(11) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(11), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,4] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(14) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(14), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(14) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(14), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,5] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(17) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(17), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(17) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(17), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,6] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(20) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(20), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(20) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(20), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,7] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(23) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(23), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(23) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(23), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,8] <- ifelse(ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(26) != 
                                    ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "e368d6d0-eeb5-42aa-8aee-4de0c8530da3") %>% select_(26), TRUE, ifelse(
                                      ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "a92f6308-5fb6-4754-a80f-5fcead380da6") %>% select_(26) ==
                                        ROB %>% filter(StudyId == Both_StudyIDs[robstudy,3] & AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b") %>% select_(26), TRUE, FALSE))
}
for (robstudy in 1:nrow(Both_StudyIDs)) {
  ROB_r_new[robstudy,9] <- Both_StudyIDs[robstudy,3] 
}

# Filtering data where the reconciled ROB was different then both ROBs from Joly Ghanawhi and Timo Sander:
ROB_r_new %>% filter(V1 == FALSE)
ROB_r_new %>% filter(V2 == FALSE)
ROB_r_new %>% filter(V3 == FALSE)
ROB_r_new %>% filter(V4 == FALSE)
ROB_r_new %>% filter(V5 == FALSE)
ROB_r_new %>% filter(V6 == FALSE)
ROB_r_new %>% filter(V7 == FALSE)
ROB_r_new %>% filter(V8 == FALSE)
# For only one ROB parameter of one article this was the case. This articles ROB parameter was changed manually.

# Getting the reconciled data via filtering for Emmas Annotator ID:
ROB_r <- ROB %>% filter(AnnotatorIdStr == "857cd48e-9080-43f8-9a91-f6ea13e1ce0b")

# Excluding later excluded articles:
ROB_r <- ROB_r %>% filter(StudyIdStr != "8d176ce8-1e76-49b8-b100-e1bc1e3f3b48" &
                            StudyIdStr != "552afe89-a702-4001-8e57-1936c83048d3" &
                            StudyIdStr != "faca9bd2-760e-4c53-9640-2771cc30c576" &
                            StudyIdStr != "9c16f6bd-a1d6-41da-829a-30465f870365")

# Ordering the reconciled ROB results and the Both_StudyIDs after the ROB StudyIDs:
ROB_r <- ROB_r[order(ROB_r$StudyId),]
Both_StudyIDs_ordered_by_ROB <- Both_StudyIDs[order(Both_StudyIDs$`StudyID ROB`),]

# Checking if the ROB StudyIDs match:
all(ROB_r$StudyId == Both_StudyIDs_ordered_by_ROB$`StudyID ROB`)
# As all values are TRUE we can continue!

# Changing the StudyIDs in the reconciled ROB results dataframe to the "correct" Annotation StudyIDs:
ROB_r$StudyId <- Both_StudyIDs_ordered_by_ROB$`StudyID Annotation`

# Cave: The StudyIDs are only changed for the reconciled ROB results!

# Changing column names:
colnames(ROB_r)[1] <- paste("StudyID annotations")
colnames(ROB_r)[3] <- paste("StudyID RoB")



### Analysis of risks of bias:

# Creating a new dataframe ROB_r_ov ("Risk of Bias reconciled overview") with the reconciled Risk of Bias:
ROB_r_ov <- as.data.frame(ROB_r[,c(1,3,seq(5,26,3))])

# Changing column names:
colnames(ROB_r_ov)[seq(3,10,1)] <- paste(c("repeats and replications","blinded outcome assessment","sample size calculation","protocol","all data",
                                           "random group allocation","all data points","mean/error calculation"))

# converting to numerical values (1 means no risk and 0 means a risk is present):
ROB_r_ov[,3:10] <- ifelse(ROB_r_ov[,3:10] == TRUE, 1, 0)

# Getting a summary score for each article:
ROB_r_ov$sum_of_risks_factors <- 8 - rowSums(ROB_r_ov[,3:10])

# Getting a summary score for each parameter:
ListROB <- list()
for (ROBp in 3:10) {
  ROBpsum <- sum(ROB_r_ov[, ROBp], na.rm = TRUE)
  ListROB[[length(ListROB)+1]] = ROBpsum
}
aROB <- data.frame(ListROB)
colnames(aROB) <- colnames(ROB_r_ov[,3:10])
bROB <- data.frame(t(aROB))
bROB$risk_prevalence <- (137 - bROB$t.aROB.) / 137
bROB_ordered <- bROB[order(bROB$risk_prevalence, decreasing = TRUE),]

# Creating a dataframe with the prevalence of each risk of bias parameter:
RoB_parameters_prevalence_in_percent <- bROB_ordered
RoB_parameters_prevalence_in_percent$in_percent <- RoB_parameters_prevalence_in_percent$risk_prevalence * 100

# Mean number of risks of bias factors present per article:
mean(ROB_r_ov$sum_of_risks_factors)






####### Additional analyses: Clinically relevant TMZ concentrations and treatment durations

### Clinically relevant TMZ concentrations

# How many effects with TMZ concentrations above clinically relevant concentrations (> 5 microM and > 10 microM)?
AR_OD_MA_T_f_RMD %>% filter(tmz_concentration_in_microM_numeric > 5) %>% nrow()
AR_OD_MA_T_f_RMD %>% filter(tmz_concentration_in_microM_numeric > 10) %>% nrow()

# What is the share of studies that included no clinically relevant concentration of Temozolomide (< 10 microM)?
Studies_clin_rel_conc <- AR_OD_MA_T_f_RMD %>% filter(tmz_concentration_in_microM_numeric <= 10) %>% select_(1) %>% unique()
Studies_clin_irrel_conc <- AR_OD_MA_T_f_RMD %>% filter(tmz_concentration_in_microM_numeric > 10) %>% select_(1) %>% unique()
Studies_clin_irrel_conc_only <- Studies_clin_irrel_conc %>% filter(!is.element(StudyIdStr, Studies_clin_rel_conc$StudyIdStr))
nrow(Studies_clin_irrel_conc_only) / 101



### Clinically relevant treatment durations

# How many treatment durations are lower than the doubling time of U87 (34 h)?
AR_OD_MA_T_f_RMD %>% filter(Time < (1.5*34/24)) %>% nrow()

# What is the share of studies that had only treatment durations that are lower than the 1.5 times the doubing time of U87 (34 h)?
Studies_clin_rel_time <- AR_OD_MA_T_f_RMD %>% filter(Time >= (1.5*34/24)) %>% select_(1) %>% unique()
Studies_clin_irrel_time <- AR_OD_MA_T_f_RMD %>% filter(Time < (1.5*34/24)) %>% select_(1) %>% unique()
Studies_clin_irrel_time_only <- Studies_clin_irrel_time %>% filter(!is.element(StudyIdStr, Studies_clin_rel_time$StudyIdStr))
nrow(Studies_clin_irrel_time_only) / 101






######## Code for creating the figures



### Reporting quality of the parameters:

## Getting the reporting quality of every parameter summarized over all studies:

# Creating a  dataframe with the calculated reporting quality for every parameter:
Listpara <- list()
for (para in 2:20) {
  parasum <- sum(ARRQSP[, para], na.rm = TRUE)
  Listpara[[length(Listpara)+1]] = parasum
}

apara <- data.frame(Listpara)
colnames(apara) <- colnames(ARRQSP[,2:20])
bpara <- data.frame(t(apara))
bpara[1:10,2] <- c("Mycoplasma exclusion", "U-87 MG age", "Source of U-87 MG", "Cell passaging criteria", "Antibiotics", "Source of FBS", "U-87 MG authentication", "Cell growth assay", "Glucose level","Source of TMZ")
bpara[11:19,2] <- c("Volume of added TMZ", "Type of control", "Volume of added control", "U-87 MG concentration", "Treatment duration", "Concentration of TMZ", "Errortype", "Number of experiments", "Conflicts of interests")
bpara$rep_percentage <- (137 - bpara$t.apara.) / 137
bpara_ordered <- bpara[order(bpara$rep_percentage, decreasing = TRUE),]
RQparameters <- bpara_ordered

# Remove the parameter cell growth assay (as reporting is 100% because it's an inclusion criterion):
RQparameters <- RQparameters[2:19,]

# Rounding the reporting percentage values:
RQparameters$rep_percentage_rounded <- round(RQparameters$rep_percentage,3)

# Plotting:
par(mar=c(15,5,2,1))
barplot_parameter_reporting_without_assay <- barplot(RQparameters$rep_percentage * 100,
                                                     main = "2a) Parameters reporting quality",
                                                     cex.main = 1.6,
                                                     col = ifelse(RQparameters$rep_percentage >= 0.5, "grey70", "grey90"),
                                                     ylab = "Percentage of reporting",
                                                     ylim = c(0,100),
                                                     cex.lab=1.4,
                                                     names.arg = RQparameters$V2, las = 2, cex.names = 1.2)
grid(NA, NULL, col = "black")

# Adding values on top of the bars:
text(x = barplot_parameter_reporting_without_assay, y = RQparameters$rep_percentage_rounded * 100, label = RQparameters$rep_percentage_rounded * 100, pos = ifelse((RQparameters$rep_percentage_rounded * 100) > 50,1,3), cex = 1, col = "black")



### Articles reporting quality per year:

# Building a inear regression model:
model_RQ_years <- lm((ARRQSP$Reporting_ratio *100) ~ ARRQSP$Year)
summary(model_RQ_years)

# Workaround: Replacing the year "2003 - 2011" with the numerical value 2011 to make the regression plotting possible. Label gets changed later.
RQ_Years[1,1] <- 2011
RQ_Years$Years <- as.numeric(RQ_Years$Years)

## Plotting it with the number of articles per year:
ggplot(data = RQ_Years, aes(x = Years, y = `Mean reporting ratio` * 100)) +
  geom_line(group = 1, size = 1) +
  geom_errorbar(aes(ymin = (`Mean reporting ratio` - `Standard deviation`) * 100, ymax = (`Mean reporting ratio` + `Standard deviation`) * 100), size = 0.7, color = "black", width = 0.2) +
  xlab("Year") + ylab("Percentage of reported parameters") +
  theme_classic() +
  theme(axis.title = element_text(size = 18)) +
  theme(axis.text = element_text(size = 16)) +
  ggtitle("2b) Articles reporting quality over time") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  geom_hline(yintercept=seq(40,70,10), linetype = "dotted", size = 0.7) +
  geom_segment(x = 2012, xend = 2020, y = coef(model_RQ_years)[1] + coef(model_RQ_years)[2] * 2012, yend = coef(model_RQ_years)[1] + coef(model_RQ_years)[2] * 2020, size = 1.2, linetype = "twodash") +
  scale_x_continuous(limits = c(2010.8,2020.2), breaks = seq(2011,2020,1), labels = c("2003 - 2011",2012:2020)) +
  scale_y_continuous(limits = c(35,72)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  geom_text(x=2015, y=67, label="Linear regression model estimate\n+0.635% reporting score per year\np = .011", size = 6) +
  geom_text(x = seq(2011,2020,1), y = 34.5, label = RQ_Years[1:10,4], size = 5.5) +
  geom_text(x = 2012.11, y = 35.5, label = expression(underline("Number of articles per year:")), size = 5.5)



### Articles reporting quality depending on the journal impact factors of the publishing journals:

# Building a linear model:
model_RQ_JIF <- lm(Reporting_ratio~JIF, data = ARRQSP, na.action = na.omit)
summary(model_RQ_JIF)

# Plotting it:
ggplot(ARRQSP, aes(x=JIF, y=Reporting_ratio*100)) + 
  geom_point(color="black", size = 2) + 
  geom_smooth(method=lm, color = "black") +
  xlab("Journal Impact Factor (JIF)") +
  ylab("Percentage of reported parameters") +
  ggtitle("2c) Relationship between the articles reporting quality\nand the JIF of publishing journal") +
  theme_classic() +
  theme(axis.title = element_text(size = 18)) +
  theme(axis.text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  scale_x_continuous(limits = c(0,10.2), breaks = seq(0,10,1)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  scale_y_continuous(limits = c(20,90), breaks = seq(20,90,10)) +
  geom_hline(yintercept=seq(20,90,10), linetype = "dotted", size = 0.7) +
  geom_text(x=1.5, y=80, label="Linear regression model\nEstimate: +1.74% per JIF-unit\np < .001", size = 6)



### TMZ concentrations:

# Extracting the TMZ concentrations:
TMZ_concentrations <- as.data.frame(AR_updated$tmz_concentration_in_microM_numeric)

# Excluding untreated controls:
TMZ_concentrations <- TMZ_concentrations %>% filter(`AR_updated$tmz_concentration_in_microM_numeric` != 0)

# calculating the decadic logarithmic TMZ concentrations:
TMZ_concentrations$dec_log_tmz_conc <- log(TMZ_concentrations$`AR_updated$tmz_concentration_in_microM_numeric`,10)

# Plotting a histogramm:
ggplot(aes(x = TMZ_concentrations$dec_log_tmz_conc), data = TMZ_concentrations) +
  geom_histogram(color="black", fill="grey") +
  theme_classic() +
  xlab("Temozolomide concentration (log. scale)") + ylab("Number of experiments") +
  scale_x_continuous(limits = c(-2,5), breaks = seq(-2,5,1), labels = c("0.01 µM","0.1 µM","1 µM", "10 µM", "100 µM", "1 000 µM", "10 000 µM", "100 000 µM")) +
  ggtitle("Temozolomide concentrations") +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.title = element_text(size = 18)) +
  theme(axis.text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold")) +
  geom_hline(yintercept=seq(50,150,50), linetype = "dotted", size = 0.7) +
  geom_vline(xintercept = median(TMZ_concentrations$dec_log_tmz_conc), size = 1.1) +
  geom_text(x=-0.3, y=100, label="Median Temozolomide concentration:\n100 µM\nLowest: 0.01 µM\nHighest: 16 000 µM", size = 6)

# What is the total number of experiments with Temozolomide concentration given?
TMZ_concentrations %>% filter(!is.na(`AR_updated$tmz_concentration_in_microM_numeric`)) %>% nrow()

# What is the number of unique Temozolomide concentrations?
TMZ_concentrations %>% select_(1) %>% unique() %>% nrow()

# Minimum TMz concentration:
min(TMZ_concentrations$`AR_updated$tmz_concentration_in_microM_numeric`)

# Maximum TMz concentration:
max(TMZ_concentrations$`AR_updated$tmz_concentration_in_microM_numeric`)

# Median TMz concentration:
median(TMZ_concentrations$`AR_updated$tmz_concentration_in_microM_numeric`)



### Treatment duration:

# Extracting the Times:
Treatment_times <- as.data.frame(AR_updated$Time)

# Adding the TMZ concentrations:
Treatment_times$TMZ_conc <- AR_updated$tmz_concentration_in_microM_numeric

# Excluding untreated control values:
Treatment_times <- Treatment_times %>% filter(TMZ_conc != 0)

# Excluding time = 0:
Treatment_times <- Treatment_times %>% filter(`AR_updated$Time` > 0)

# Plotting a histogram:
ggplot(aes(x = Treatment_times$`AR_updated$Time`), data = Treatment_times) +
  geom_histogram(color="black", fill="grey", binwidth = 1) +
  theme_classic() +
  xlab("Treatment duration [days]") + ylab("Number of experiments") +
  scale_x_continuous(limits = c(-0,12.5), breaks = seq(0,12,1)) +
  scale_y_continuous(limits = c(0,300), breaks = seq(50,300,50)) +
  ggtitle("Treatment durations") +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.title = element_text(size = 18)) +
  theme(axis.text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold")) +
  geom_hline(yintercept=seq(50,300,50), linetype = "dotted", size = 0.7) +
  geom_vline(xintercept = median(Treatment_times$`AR_updated$Time`, na.rm = T), size = 1.1) +
  geom_text(x=7, y=200, label="Median treatment duration:\n3 days\nShortest: 4 hours\nLongest: 12 days", size = 6)

# What is the total number of experiments with Temozolomide?
Treatment_times %>% filter(`AR_updated$Time` > 0) %>% nrow()

# What is the number of unique Treatment durations?
Treatment_times %>% select_(1) %>% unique() %>% nrow()

# Minimum treatment duration:
min(Treatment_times$`AR_updated$Time`)

# Maximum treatment duration:
max(Treatment_times$`AR_updated$Time`)

# Median treatment duration:
median(Treatment_times$`AR_updated$Time`)





######## THE END ####################
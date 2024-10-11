#==========================================================================================
# Circulating amino acids in CAD patients vs. healthy controls (table1)
# Author: Justin Carrard
#==========================================================================================
#------------------------------------------------------------------------------------------
# Reminder 
#------------------------------------------------------------------------------------------

# sex = 1 -> male, sex = 0 -> female

#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------
#install.packages("htmlTable")
library(htmlTable)
#install.packages("table1")
library(table1)
library(magrittr)
library(boot) 
library(knitr)
library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(transplantr)

#install.packages("rmarkdown")
library(rmarkdown)

#install.packages("writexl")
library("writexl")

#install.packages("data.table")
library(data.table)

#install.packages("MatchIt")
library(MatchIt)

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#------------------------------------------------------------------------------------------
# Define paths
#------------------------------------------------------------------------------------------

# Paths
setwd("/Volumes/nuxyde32/Uni-Basel/DSBG/Forschung/COmPLETE Health/Metabolomics/MSc Luisa/Table 1") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

#------------------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------------------

#import_data
dat <- as.data.frame(read_excel(paste0(data_path, "/", "CHL_CHR_complete.xlsx")))

head(dat)
str(dat)

#------------------------------------------------------------------------------------------
# Categorise sampling time
#------------------------------------------------------------------------------------------

#remove wrong date from sampling time column
dat$"Sampling_time" = gsub("1899-12-31","",dat$"Sampling_time")

#Group sampling time into 5 categories
dat$Sampling_time_cat <- cut(
  chron::times(dat$"Sampling_time")
  , breaks = chron::times(c(
    "08:00:00"
    , "10:00:00"
    , "12:00:00"
    , "14:00:00"
    , "16:00:00"
    , "23:59:59"
  ))
  , labels = c("08-09.59","10-11.59","12-13.59","14-15.59","16-18.59")
  , include.lowest = TRUE
  , right = TRUE
)

dat[, c("Sampling_time", "Sampling_time_cat")]

#------------------------------------------------------------------------------------------
# Convert categorical data to factors
#------------------------------------------------------------------------------------------

#Convert Sex and CAD to factor
dat$Sex <- factor(dat$Sex)
dat$Smoking_status <- factor(dat$Smoking_status)
dat$CAD <- factor(dat$CAD)
dat$LVEF <- factor(dat$LVEF)
dat$Number_of_vessels <- factor(dat$Number_of_vessels)
dat$Sampling_time_cat <- factor(dat$Sampling_time_cat)

#------------------------------------------------------------------------------------------
# Add eGFR based on creatinine values in micromol/L
#------------------------------------------------------------------------------------------

dat$eGFR <- ckd_epi(creat = dat$creatinine, age = dat$Age, sex = dat$Sex, eth = rep("non-black", nrow(dat)))

#------------------------------------------------------------------------------------------
# Add a stage of kidney function (1 to 5) according to eGFR
#------------------------------------------------------------------------------------------

dat$eGFR_stage <- NA

dat$eGFR_stage[dat$eGFR > 90] <- "G1"
dat$eGFR_stage[dat$eGFR <= 90 & dat$eGFR > 60] <- "G2"
dat$eGFR_stage[dat$eGFR <= 60 & dat$eGFR > 45] <- "G3a"
dat$eGFR_stage[dat$eGFR <= 45 & dat$eGFR > 30] <- "G3b"
dat$eGFR_stage[dat$eGFR <= 30 & dat$eGFR > 15] <- "G4"
dat$eGFR_stage[dat$eGFR<= 15] <- "G5"

dat$eGFR_stage <- factor(dat$eGFR_stage, levels = c("G1", "G2", "G3a", "G3b", "G4", "G5"))

#Convert eGFR_stage to factor
dat$eGFR_stage <- factor(dat$eGFR_stage)

#count the frequency of each stages
table(dat$eGFR_stage)

#------------------------------------------------------------------------------------------
# Check for missing data
#------------------------------------------------------------------------------------------

#Check missing data
colSums(is.na(dat))

#------------------------------------------------------------------------------------------
# Match CAD patients with healthy controls
#------------------------------------------------------------------------------------------

# No matching; constructing a pre-match matchit object
m.out0 <- matchit(CAD ~ Age + Sex, data = dat,
                  method = NULL, distance = "glm")

# Checking balance prior to matching
summary(m.out0)

#plot data prior matching
library(ggplot2)
p <- ggplot(dat, aes(fill = factor(CAD))) + 
  geom_bar(position = "dodge") + 
  scale_fill_discrete("CAD")
p + aes(x = Age)
p + aes(x = Sex)

# Exact matching on a probit PS
m.out1 <- matchit(CAD ~ Age + Sex, data = dat, estimand = "ATE",
                  method = "cem")
m.out1

#Checking balance after matching by numbers
summary(m.out1)
plot(summary(m.out1))

#Estimate treatment effect by numbers
m.data1 <- match.data(m.out1)

head(m.data1)

#Create matched data frame
dat.match <- match.data(m.out1)

#Check missing data
colSums(is.na(dat.match))

#------------------------------------------------------------------------------------------
# Compute table 1
#------------------------------------------------------------------------------------------

#Specify categorical variables
dat.match$Sex <- 
  factor(dat.match$Sex, levels=c(0,1),
         labels=c("Female", 
                  "Male"))

dat.match$Smoking_status <- 
  factor(dat.match$Smoking_status, levels=c(0, 1, 2, 3),
         labels=c("Never smoked", 
                  "Former smoker, stopped > 10 years", 
                  "Former smoker, stopped < 10 years", 
                  "Active smokers"))

dat.match$CAD <- 
  factor(dat.match$CAD, levels=c(0,1),
         labels=c("Healthy individuals", 
                  "CAD patients"))

dat.match$LVEF <- 
  factor(dat.match$LVEF, levels=c(0, 1, 2, 3),
         labels=c("Normal",
                  "Preserved (≥50%)", 
                  "Mildly reduced (40-49%)", 
                  "Reduced (≤40%)"))

dat.match$Number_of_vessels <- 
  factor(dat.match$Number_of_vessels, levels=c(0, 1, 2, 3),
         labels=c("None",
                  "Single-vessel CAD", 
                  "Double-vessel CAD", 
                  "Triple-vessel CAD"))

dat.match$eGFR_stage <- 
  factor(dat.match$eGFR_stage, levels=c("G1", "G2", "G3a", "G3b", "G4", "G5"),
         labels=c("Normal or high", 
                  "Mildly decreased", 
                  "Mildly to moderately decreased", 
                  "Moderately to severely decreased",
                  "Severely decreased",
                  "Kidney failure"))


#Specify units
units(dat.match$Age) <- "years"
units(dat.match$BMI) <- "kg/m^2"
units(dat.match$PBF) <- "%"
units(dat.match$Smoking_status) <- "n(%)"
units(dat.match$VO2peak_ml_min_kg) <- "mL/min/kg"
units(dat.match$Total_MVPA) <- "min/day"
units(dat.match$Total_PA) <- "min/day"
units(dat.match$SBP) <- "mmHg"
units(dat.match$DBP) <- "mmHg"
units(dat.match$Chol_mmol_l) <- "mmol/L"
units(dat.match$HDL_mmol_l) <- "mmol/L"
units(dat.match$LDL_mmol_l) <- "mmol/L"
units(dat.match$TG_mmol_l) <- "mmol/L"
units(dat.match$HbA1c) <- "%"
units(dat.match$eGFR_stage) <- "ml/min/1.73m2"

#Rename variables
label(dat.match$Smoking_status) <- "Smoking status"
label(dat.match$BMI) <- "Body mass index"
label(dat.match$PBF) <- "Percentage of body fat"
label(dat.match$SBP) <- "Systolic blood pressure"
label(dat.match$DBP) <- "Diastolic blood pressure"
label(dat.match$VO2peak_ml_min_kg) <- "VO2peak"
label(dat.match$Total_MVPA) <- "Daily moderate to vigorous PA"
label(dat.match$Total_PA) <- "Daily total PA"
label(dat.match$Chol_mmol_l) <- "Total cholesterol"
label(dat.match$HDL_mmol_l) <- "HDL-cholesterol"
label(dat.match$LDL_mmol_l) <- "LDL-cholesterol"
label(dat.match$TG_mmol_l) <- "Triglycerides"
label(dat.match$eGFR_stage) <- "Glomerular filtration rate categories"

#Compute Table1
Table1 <- table1(~ Age + PBF + BMI + SBP + DBP + VO2peak_ml_min_kg + Total_MVPA + Total_PA +Smoking_status + Fasting + Chol_mmol_l + HDL_mmol_l + LDL_mmol_l + TG_mmol_l + HbA1c + eGFR_stage + Number_of_vessels + LVEF | CAD*Sex, data=dat.match, overall = TRUE)

my.render.cont <- function(x) {
    with(stats.default(x), 
         c("",
           
          "Mean (SD)" = sprintf("%s (%s)",
                                round_pad(MEAN, 1),
                                round_pad(SD, 1)),
         
          "Median (Min, Max)" = sprintf("%s (%s, %s)",
                                       round_pad(MEDIAN, 1), 
                                       round_pad(MIN, 1), 
                                       round_pad(MAX, 1)))
    )
}

#Compute Table1
Table1 <- table1(~ Age + PBF + BMI + SBP + DBP + VO2peak_ml_min_kg + Total_MVPA + Total_PA +Smoking_status + Fasting + Chol_mmol_l + HDL_mmol_l + LDL_mmol_l + TG_mmol_l + HbA1c + eGFR_stage + Number_of_vessels + LVEF | CAD*Sex, data=dat.match, overall = TRUE)

# Convert to HTML
Table1_html <- htmlTable(as.matrix(Table1))

# Save to HTML file
writeLines(Table1_html, con = "./output/text/Table1.html")


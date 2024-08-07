#==========================================================================================
# Circulating amino acids in CAD patients vs. healthy controls
# Author: Luisa, Prechtl, Justin Carrard & Denis Infanger
#==========================================================================================
#------------------------------------------------------------------------------------------
# Reminder 
#------------------------------------------------------------------------------------------

# sex = 1 -> male, sex = 0 -> female

#------------------------------------------------------------------------------------------
# Load packages
# Import library
#------------------------------------------------------------------------------------------
library(readxl)
library(car)
library(broom)
library(dplyr)
library(tidyr)
library(ggdendro)
library(egg)
library(ggplot2)
library(reshape2)
library(lmtest)
library(sandwich)

#install.packages("multcomp")
library(multcomp)  # for glht

#install.packages("writexl")
library("writexl")

#install.packages("data.table")
library(data.table)

#install.packages("ccoptimalmatch")
library(ccoptimalmatch)

#install.packages("MatchIt")
library(MatchIt)

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#------------------------------------------------------------------------------------------
# Define paths
#------------------------------------------------------------------------------------------

# Paths
setwd("") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

#------------------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------------------

#import_data
dat <- as.data.frame(read_excel(paste0(data_path, "/", "xxx.xlsx")))

# View(dat)
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
dat$CAD <- factor(dat$CAD)
dat$Sampling_time_cat <- factor(dat$Sampling_time_cat)

#------------------------------------------------------------------------------------------
# Check for missing data
#------------------------------------------------------------------------------------------

#Check missing data
colSums(is.na(dat))

#------------------------------------------------------------------------------------------
# Log2 transform data
#------------------------------------------------------------------------------------------

#log2_transformation_AminoAcids
(variables_to_transform <- names(dat)[37:68]) # Name of variables to transform with log2

for (i in variables_to_transform) {
  dat[, paste0(i, "_log2")] <- log2(dat[, i])
}

names(dat)

head(dat)

#------------------------------------------------------------------------------------------
# Descriptive scatterplots
#------------------------------------------------------------------------------------------

# Create scatterplot matrix to inspect correlation among clinical data
png(paste("scatterplot_matrix_age.png", sep = ""), width = 9*2.5, height = 9*2.5, units = "in", res = 300)

car::scatterplotMatrix(~Age + Sex + VO2peak_ml_min_kg + Total_PA,  data = dat
                       , diagonal = list(method = "boxplot")
                       , smooth = list(method = loessLine, spread = FALSE, col.smooth = "red")
                       , by.groups = TRUE
                       , use = "pairwise.complete.obs"
                       , regLine  = FALSE
                       # , lwd = 2
                       , col = c("#00ABCE", "violet")
                       # , pch = c(15, 1)
                       # , cex = 1.5
)

dev.off()

#------------------------------------------------------------------------------------------
# Data standardisation
#------------------------------------------------------------------------------------------

#z_standardisation_dependent_variables
(variables_to_standardize <- names(dat)[c(2, 21:30, 69, 71:102)])

for (i in variables_to_standardize) {
  dat[, paste0(i, "_std")] <- scale(dat[, i], center = TRUE, scale = TRUE)
}

names(dat)

head(dat)

#------------------------------------------------------------------------------------------
# Plot descriptive graphics
#------------------------------------------------------------------------------------------

# Some descriptive graphics
(vars_to_plot <- names(dat)[c(114:146)])

# Reshape to long
dat_long <- reshape2::melt(
  dat
  , id.vars = c("Sex", "Age", "CAD")
  , measure.vars = vars_to_plot
)

#theme_CAD(theme_bw())
p <- ggplot(data = dat_long, aes(x = CAD, y = value, group = CAD)) +
  geom_violin(aes(fill = CAD), trim = TRUE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3) +
  xlab("Sex") +
  ylab("Value") +
  scale_fill_manual(breaks = c(0, 1), values = c("#00ACEE", "#F4A749")) +
  facet_wrap(~variable, scales = "free_y") +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    legend.position="right",
    legend.text=element_text(size=12.5),
    legend.key=element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text.x=element_text(size=15)
  )

p

ggsave(paste(graphics_path, paste("violin_cad.png", sep = ""), sep = "/"), p, width = 18*1, height = 14*1, units = "in", dpi = 300)

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
# Multiple linear regressions
#------------------------------------------------------------------------------------------
my_lms <- lapply(dat.match[,c(114:146)], function(x) lm(x ~ CAD + Age_std * Sex + Fasting_std + Sampling_time_cat + Total_PA_std, data = dat.match, weights = weights))

coef_results <- lapply(my_lms, FUN = function(x) coeftest(x, vcov = vcovCL(x, cluster = dat.match$subclass))) # cluster-robust standard errors

#------------------------------------------------------------------------------------------
# Summaries & plot
#------------------------------------------------------------------------------------------
#summaries
lapply(my_lms, summary)

#------------------------------------------------------------------------------------------
# Check residuals
#------------------------------------------------------------------------------------------

# Get resid for all models
list_resid <- lapply(my_lms, resid)

# If you want them in a data.frame instead of list
df_resid <- do.call(cbind.data.frame, list_resid)

#plot residuals
p <- par(mfrow = c(2, 2))
lapply(my_lms, plot)
p

#------------------------------------------------------------------------------------------
# Obtain 95% confidence intervalls
#------------------------------------------------------------------------------------------

# Use "coefci" to calculate robust confidence intervals
conf_intervals <- lapply(my_lms, FUN = function(x) coefci(x, vcov. = vcovCL(x, cluster = dat.match$subclass))) # cluster-robust standard errors

# Print the results
print(conf_intervals)

# Assuming each confidence interval is a matrix with appropriate row names
convert_to_dataframe <- function(ci, model_index) {
  # Check if ci is a matrix and convert to dataframe
  if (is.matrix(ci)) {
    ci_df <- as.data.frame(ci)
  } else {
    ci_df <- ci  # if it's already a dataframe
  }
  
  # Ensure that rownames can be converted to a column
  ci_df$Term <- rownames(ci_df)
  rownames(ci_df) <- NULL  # Clean up row names
  
  # Name the bounds columns if not named
  names(ci_df)[1:2] <- c("2.5%", "97.5%")
  
  # Add a column for the model index
  ci_df$Model <- sprintf("Model %d", model_index)

  return(ci_df)
}

# Apply this function to each element in the list and combine into one dataframe
df_conf_intervals <- do.call(rbind, lapply(seq_along(conf_intervals), function(i) {
  convert_to_dataframe(conf_intervals[[i]], i)
}))

# View the first few rows to verify
head(df_conf_intervals, 20)

# Select rows of interest
df_conf_intervals2 <- subset(df_conf_intervals, Term %in% c("CAD1", "Age_std", "Sex1", "Age_std:Sex1", "Total_PA_std"))

# Replace Model name by effective amino acid names
df_conf_intervals3 <- df_conf_intervals2 %>%
  mutate(Model = case_when(
    Model == "Model 1" ~ "glutamine/glutamate",
    Model == "Model 2" ~ "creatinine",
    Model == "Model 3" ~ "phenylalanine",
    Model == "Model 4" ~ "leucine",
    Model == "Model 5" ~ "kynurenine",
    Model == "Model 6" ~ "tryptophan",
    Model == "Model 7" ~ "isoleucine",
    Model == "Model 8" ~ "methionine",
    Model == "Model 9" ~ "taurine",
    Model == "Model 10" ~ "valine",
    Model == "Model 11" ~ "proline",
    Model == "Model 12" ~ "pipecolate",
    Model == "Model 13" ~ "tyrosine",
    Model == "Model 14" ~ "alpha-aminobutyrate",
    Model == "Model 15" ~ "beta-alanine",
    Model == "Model 16" ~ "creatine",
    Model == "Model 17" ~ "alanine",
    Model == "Model 18" ~ "hydroxyproline",
    Model == "Model 19" ~ "guanidinoacetate",
    Model == "Model 20" ~ "threonine",
    Model == "Model 21" ~ "2-aminoadipate",
    Model == "Model 22" ~ "glycine",
    Model == "Model 23" ~ "glutamate",
    Model == "Model 24" ~ "serine",
    Model == "Model 25" ~ "glutamine",
    Model == "Model 26" ~ "asparagine",
    Model == "Model 27" ~ "homocitrulline",
    Model == "Model 28" ~ "citrulline",
    Model == "Model 29" ~ "aspartate",
    Model == "Model 30" ~ "arginine",
    Model == "Model 31" ~ "histidine",
    Model == "Model 32" ~ "lysine",
    Model == "Model 33" ~ "ornithine",
    TRUE ~ Model  # Keep other values unchanged
  ))
# View the resulting dataset
print(df_conf_intervals3)

# Replace term names by effective names
df_conf_intervals4 <- df_conf_intervals3 %>%
  mutate(Term = case_when(
    Term == "CAD1" ~ "CAD",
    Term == "Sex1" ~ "Male",
    Term == "Age_std" ~ "Age",
    Term == "Age_std:Sex1" ~ "Interaction age:male",
    Term == "Total_PA_std" ~ "Daily physical activity",
       TRUE ~ Term  # Keep other values unchanged
  ))
# View the resulting dataset
print(df_conf_intervals4)

#------------------------------------------------------------------------------------------
# Extract variable from my_lms
#------------------------------------------------------------------------------------------

#extract CAD
res_frame_CAD <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  res_frame_CAD$estimate[which(res_frame_CAD$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["CAD1", "Estimate"]
  res_frame_CAD$std.error[which(res_frame_CAD$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["CAD1", "Std. Error"]
  res_frame_CAD$statistic[which(res_frame_CAD$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["CAD1", "t value"]
  res_frame_CAD$p.value[which(res_frame_CAD$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["CAD1", "Pr(>|t|)"]
  res_frame_CAD$response[which(res_frame_CAD$term %in% names(my_lms)[i])] <- "CAD"

}  

#extract Age
res_frame_Age <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  res_frame_Age$estimate[which(res_frame_Age$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std", "Estimate"]
  res_frame_Age$std.error[which(res_frame_Age$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std", "Std. Error"]
  res_frame_Age$statistic[which(res_frame_Age$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std", "t value"]
  res_frame_Age$p.value[which(res_frame_Age$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std", "Pr(>|t|)"]
  res_frame_Age$response[which(res_frame_Age$term %in% names(my_lms)[i])] <- "Age"
  
}  

#extract Sex
res_frame_Sex <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  res_frame_Sex$estimate[which(res_frame_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Sex1", "Estimate"]
  res_frame_Sex$std.error[which(res_frame_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Sex1", "Std. Error"]
  res_frame_Sex$statistic[which(res_frame_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Sex1", "t value"]
  res_frame_Sex$p.value[which(res_frame_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Sex1", "Pr(>|t|)"]
  res_frame_Sex$response[which(res_frame_Sex$term %in% names(my_lms)[i])] <- "Male"
  
}  

#extract Age:Sex results
res_frame_Age_Sex <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {

  res_frame_Age_Sex$estimate[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std:Sex1", "Estimate"]
  res_frame_Age_Sex$std.error[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std:Sex1", "Std. Error"]
  res_frame_Age_Sex$statistic[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std:Sex1", "t value"]
  res_frame_Age_Sex$p.value[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Age_std:Sex1", "Pr(>|t|)"]
  res_frame_Age_Sex$response[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- "Interaction age:male"
  
}

#extract Total_PA results
res_frame_Total_PA <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {

  res_frame_Total_PA$estimate[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Total_PA_std", "Estimate"]
  res_frame_Total_PA$std.error[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Total_PA_std", "Std. Error"]
  res_frame_Total_PA$statistic[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Total_PA_std", "t value"]
  res_frame_Total_PA$p.value[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- coef_results[[names(my_lms)[i]]]["Total_PA_std", "Pr(>|t|)"]
  res_frame_Total_PA$response[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- "Daily physical activity"
  
}

#------------------------------------------------------------------------------------------
# Combine extracted data into a data frame
#------------------------------------------------------------------------------------------

#Combine_data_frames
Overall <- rbind(
   res_frame_CAD
  , res_frame_Age
  , res_frame_Sex
  , res_frame_Age_Sex
  , res_frame_Total_PA
)

#------------------------------------------------------------------------------------------
# Adjust p-values for multiple testing
#------------------------------------------------------------------------------------------

Overall <- Overall[order(Overall$p.value),]

Overall$BH <- p.adjust(Overall$p.value, method = "BH")

#Categorise BH p-values
Overall$BH_cat <- NA

Overall$BH_cat[Overall$BH > 0.05] <- "> 0.05"
Overall$BH_cat[Overall$BH <= 0.05 & Overall$BH > 0.01] <- "≤ 0.05"
Overall$BH_cat[Overall$BH <= 0.01 & Overall$BH > 0.001] <- "≤ 0.01"
Overall$BH_cat[Overall$BH <= 0.001 & Overall$BH > 0.0001] <- "≤ 0.001"
Overall$BH_cat[Overall$BH <= 0.0001] <- "≤ 0.0001"

Overall$BH_cat <- factor(Overall$BH_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

#Order by decreasing estimates
Overall <- Overall[order(-Overall$estimate),]

#remove suffixes from amino acids name
Overall$term=gsub("_std","",Overall$term)
Overall$term=gsub("_log2","",Overall$term)

#preparing merging : rename columns of df_conf_intervals4 to match Overall columns
names(df_conf_intervals4)[3] <- "response"
names(df_conf_intervals4)[4] <- "term"

#merge statistical results with 95% confidence intervals
Overall_final <- merge(Overall, df_conf_intervals4, by = c("term", "response"))  

#rename columns
names(Overall_final)[1] <- "Dependent variable"
names(Overall_final)[2] <- "Independent variables"
names(Overall_final)[3] <- "β coefficient"
names(Overall_final)[4] <- "standard error"
names(Overall_final)[6] <- "p-value"
names(Overall_final)[7] <- "BH p-value"
names(Overall_final)[8] <- "Categorical BH p-value"
names(Overall_final)[9] <- "95% CI lower bound for β coefficient"
names(Overall_final)[10] <- "95% CI higher bound for β coefficient"

#reorder columns order to have 95% CI lower and higher bounds directly after β coefficient
Overall_final <- Overall_final[, c(1, 2, 3, 9, 10, 4, 5, 6, 7, 8)]
head(Overall_final)

write_xlsx(Overall_final, text_path, "Regressions_CAD.xlsx")

#------------------------------------------------------------------------------------------
# Start rain plot
#------------------------------------------------------------------------------------------

#Import data
plot_data <- Overall

# Define theme and palette

## Palette

palette <-
  # Blue
  c("#053061",
    "#313695",
    "#4575b4",
    "#74add1",
    "#abd9e9",
    "#e0f3f8",
    "#fee090",
    "#fdae61",
    "#f46d43",
    "#d73027",
    "#a50026",
    '#67001f')
# Red

# Calculate symmetric limits based on most extreme value
max_abs_estimate <- max(abs(plot_data$estimate))

max_lim <- max_abs_estimate
min_lim = -1 * max_lim

## theme

thm <-
  # Good starting theme + set text size
  theme_light(base_size = 7) +
  theme(
    # Remove axis ticks and titles
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # Remove gridlines and boxes
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.key = element_blank(),
    
    # White backgrounds
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white'),
    legend.background = element_rect(fill = 'white'),
    
    # Angle text
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )

#------------------------------------------------------------------------------------------
# Bare-bones rainplot
#------------------------------------------------------------------------------------------

plot_data$response <- factor(plot_data$response, levels = c("Age", "Male", "Interaction age:male", "Daily physical activity", "CAD"))

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))

print(rainplot)

#------------------------------------------------------------------------------------------
# Basic Rainplot
#------------------------------------------------------------------------------------------

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))  +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  scale_color_gradientn(
    'Effect Size\n(β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)

#------------------------------------------------------------------------------------------
# Ordering by Cluster
#------------------------------------------------------------------------------------------

# Convert to matrix and reshape for clustering.
cluster_data <-
  plot_data %>%
  dplyr::select(response, term, estimate) %>%
  spread(response, estimate)

rnms <-
  cluster_data$term

cluster_data <-
  cluster_data %>%
  dplyr::select(-term) %>%
  as.matrix()

rownames(cluster_data) <- rnms

# cluster dependent variable terms
clust <- hclust(dist(cluster_data), method = 'ward.D2')

# `clust$order` orders `term` into clusters
term_order <-
  clust$labels[clust$order]

# Convert term to a factor, ordered by `term_order`
plot_data_clo <-
  plot_data %>%
  mutate(term = factor(term, levels = term_order))


rainplot <-
  # Use cluter ordered data
  ggplot(plot_data_clo) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat)) +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  scale_color_gradientn(
    'Effect Size\n(β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)

#------------------------------------------------------------------------------------------
# Adding dendrograms
#------------------------------------------------------------------------------------------

dendro_dat <- segment(dendro_data(clust))


dendro <-
  # Empty ggplot with same y-scale as rainplot
  ggplot() +
  geom_blank(aes(y = term), data = plot_data) +
  theme_dendro() +
  # 'expand' controls whitespace around the dendrogram. The non-zero argument
  # may need to be increasesed if the line thickness of the dendrogram is
  # increased to make sure the entire dendrogram is plotted
  scale_x_discrete(position = 'top', expand = c(0, 0.03, 0, 0)) +
  # Draw dendrogram
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend),
               colour = 'black',
               data = dendro_dat)


p <- ggarrange(dendro, rainplot, ncol = 2, widths = c(1, 5))

p

ggsave(paste(graphics_path, paste("rain_plot_CAD.png", sep = ""), sep = "/"), p, width = 4*1, height = 9*1, units = "in", dpi = 300)


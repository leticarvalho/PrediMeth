##########################################################
#
# Methylation clocks (discovery cohort) - Extra
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, we investigate if there is association
# between methylation clocks and adherence to MedDiet 
# assessed by predimed_score as continuous variable,
# in the discovery cohort. 
# We also stratify the analysis by incident CVD. 

# INPUT --------------------------------------------------
# 	Discovery cohort phenotype data (diet adherence + covariables)
# 	Calculated methylation clocks for discovery cohort

# OUTPUT -------------------------------------------------
# 	Forest plots
#   Forest plot (stratified incident CVD)


##########################################################

#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(ggplot2)

## Input paths -------------------------------------------

path_to_discovery_pheno <- "" #.R
path_to_meth_clocks_discovery <- "" #.csv

# To assess incident CVD cases (for stratification)
questionario_dir <- ""

## Output
results_folder <- ""
dir.create(results_folder)


##########################################################

#### LOAD DATA ####

# Discovery phenotype data (.R)
load(path_to_discovery_pheno)

# Methylation clocks (discovery)
clocks <- read.csv(path_to_meth_clocks_discovery)

#### PREPROCESSING DATA ####

## Preprocessing phenotype data frame
str(pheno)
pheno_sub <- pheno[, c("sample_id", "entity_id", 
                       "predimed_score", "predimed_high", 
                       "age", "sex", "smoking_type", "bmi", "batch",
                       "Bcell", "CD4T", "CD8T", "Mono", "NK")]
str(pheno_sub)

## Preprocessing methylation clocks
clocks$X <- NULL
clocks <- clocks %>% rename(sample_id = id)
names(clocks)
clocks <- clocks[, c("sample_id", "age", 
                     "Levine", "ageAcc.Levine", "ageAcc3.Levine",
                     "Horvath", "ageAcc.Horvath", "ageAcc3.Horvath",
                     "Hannum", "ageAcc.Hannum", "ageAcc2.Hannum",
                     "BLUP", "ageAcc.BLUP", "ageAcc3.BLUP",
                     "EN", "ageAcc.EN", "ageAcc3.EN")]

## Merge
clock_merged <- merge(pheno_sub, clocks, by = "sample_id")
names(clock_merged)
str(clock_merged)

##########################################################

#### STUDYING ASSOCIATIONS ####

# MODEL: AccAge ~ PREDIMED score + SEX + SMOKING_HABIT + BATCH + CELL COMPOSITION #

#### Levine #### 
modLevine_cont <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                       Bcell + CD4T + CD8T + Mono + NK, 
                     data = clock_merged)
summary(modLevine_cont)

#### Horvath #### 
modHorvath_cont <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                        Bcell + CD4T + CD8T + Mono + NK, 
                      data = clock_merged)
summary(modHorvath_cont)

#### Hannum #### 
modHannum_cont <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                       Bcell + CD4T + CD8T + Mono + NK, 
                     data = clock_merged)
summary(modHannum_cont)

#### BLUP #### 
modBLUP_cont <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                     Bcell + CD4T + CD8T + Mono + NK, 
                   data = clock_merged)
summary(modBLUP_cont)

#### EN #### 
modEN_cont <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                   Bcell + CD4T + CD8T + Mono + NK, 
                 data = clock_merged)
summary(modEN_cont)


###############################################################

#### STORE MODEL RESULTS ####

## List with methylation clocks names
clocks <- c("Levine", "Horvath", "Hannum", "BLUP", "EN")

## List with models names
clocks_cont <- sapply(clocks, function(clock) paste0("mod", clock, "_cont"))

## Empty variables to store models results
beta <- c()
se <- c()
pval <- c()
ci_low <- c()
ci_high <- c()

## For loop in all models to store results
for (c in 1:length(clocks_cont)) {
  mod <- get(clocks_cont[c])
  beta[c] <- summary(mod)$coefficients["predimed_score", "Estimate"]
  se[c] <- summary(mod)$coefficients["predimed_score", "Std. Error"]
  pval[c] <- summary(mod)$coefficients["predimed_score", "Pr(>|t|)"]
  ci_low[c]  <- beta[c] - 1.96 * se[c]
  ci_high[c] <- beta[c] + 1.96 * se[c]
}

## Store results in a data frame
results_cont <- data.frame(clock = clocks, beta, se, pval, ci_low, ci_high)

## Add significance and direction labels
results_cont <- results_cont %>%
  mutate(sig = case_when(pval < 0.001 ~ "***", pval < 0.01  ~ "**", pval < 0.05  ~ "*", TRUE ~ "ns"),
         direction = ifelse(beta < 0, "Lower AccAge (protective)", "Higher AccAge (detrimental)"))


#################################################################

#### FOREST PLOT CONSTRUCTION ####

## Color by significance
sig_colors <- c("***" = "#d73027", "**" = "#f46d43", "*" = "#fdae61", "ns" = "#878787")

## Forest plot (continuous score)
forest_plot_cont <- ggplot(results_cont, aes(x = beta, y = clock, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, linewidth = 0.9) +
  geom_point(size = 4, shape = 18) +
  geom_text(aes(x = max(ci_high) + 0.3, label = sprintf("β = %.2f %s", beta, sig)), 
            hjust = 0, size = 3.2, color = "grey20", lineheight = 1.3) +
  scale_color_manual(values = sig_colors, name = "Significance", 
                     breaks = c("***", "**", "*", "ns"),
                     labels = c("p < 0.001", "p < 0.01", "p < 0.05", "n.s.")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(title = "Association of High Mediterranean Diet Adherence\nwith Epigenetic Age Acceleration (Methylation clocks)", 
       subtitle = "Linear regression: AccAge ~ predimed_score (continuous) + covariates",
       x = "Beta Coefficient (years)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(15, 80, 15, 15))

print(forest_plot_cont)


###############################################################

#### DEFINING INCIDENT CVD VARIABLE - IN DISCOVERY COHORT ####

# Diseases
conditions_1 <- fread(file.path(questionario_dir, "conditions_nivell1.csv"))
conditions_fu <- fread(file.path(questionario_dir, "fu_conditions_nivell1.csv"))

# Diseases (baseline)
conditions1_subset <- subset(conditions_1, select = c("entity_id", "DISEASE", "N"))
conditions1_subset$DISEASE_baseline <- conditions1_subset$DISEASE
conditions1_subset$DISEASE <- NULL
conditions1_subset$Num_DISEASE_baseline <- conditions1_subset$N
conditions1_subset$N <- NULL

# Diseases (follow-up)
conditions_fu_subset <- subset(conditions_fu, select = c("entity_id", "DISEASE", "EDAD"))
conditions_fu_subset$DISEASE_followup <- conditions_fu_subset$DISEASE
conditions_fu_subset$DISEASE <- NULL
conditions_fu_subset$EDAD_disease_followup <- conditions_fu_subset$EDAD
conditions_fu_subset$EDAD <- NULL

## Information about diseases (1 row per disease)
malaltias <- merge(conditions1_subset, conditions_fu_subset, by = "entity_id") 
malaltias$entity_id <- stringr::str_remove(malaltias$entity_id, "=*")

## Subset malaltias to only see samples in pheno discovery
malaltias <- malaltias[malaltias$entity_id %in% pheno_sub$entity_id,]

## Incident CVD in discovery cohort
incidentCVD <- c("HTA", "ICTUS", "ANGINA", "DIABETES", "INFARTO", "HIPERCOLESTEROLEMIA")
incidentCVD_ids <- c()

## For loop that evaluates each sample from the discovery group (pheno data frame)
for(r in 1:nrow(pheno_sub)) {
  id <- pheno_sub[r,]$entity_id # string with individual ID
  if(id %in% malaltias$entity_id) { # if this ID is in the table, this person has a disease 
    malaltias_id <- malaltias[malaltias$entity_id==id,] # data frame with all diseases related to this person (ID)
    for(m in 1:nrow(malaltias_id)) { # for each disease related to this person...
      if(malaltias_id[m,]$DISEASE_followup %in% incidentCVD) { # if this disease is in the group of incident CVD, 
        incidentCVD_ids <- append(incidentCVD_ids, id) # append this ID in the list
      }
    }
  }
}
length(incidentCVD_ids) # how many IDs were added (there can be individuals added more than once, for having more than 1 incident CVD)
incidentCVD_ids <- unique(incidentCVD_ids) # unique IDs
length(incidentCVD_ids) # 115/974 have any incident CVD

## Add this information as variable, to classify incident CVD
clock_merged$incidentCVD <- ifelse(clock_merged$entity_id %in% incidentCVD_ids, 
                                "yes", # have anyone of the incident CVDs
                                "no") # doesn't have any incident CVD
table(clock_merged$incidentCVD) # check

################################################################

#### STRATIFICATION BY INCIDENT CVD - DATA PREPARATION #### 

## Split sample by incident CVD
data_cvd_yes <- subset(clock_merged, incidentCVD == "yes")
data_cvd_no  <- subset(clock_merged, incidentCVD == "no")

## Run models separately

# Levine
Levine_cont_CVDyes <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                       Bcell + CD4T + CD8T + Mono + NK, 
                     data = data_cvd_yes)
Levine_cont_CVDno <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                           Bcell + CD4T + CD8T + Mono + NK, 
                         data = data_cvd_no)

# Hannum 
Hannum_cont_CVDyes <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                            Bcell + CD4T + CD8T + Mono + NK, data = data_cvd_yes)
Hannum_cont_CVDno <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                           Bcell + CD4T + CD8T + Mono + NK, data = data_cvd_no)

# Horvath
Horvath_cont_CVDyes <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                        Bcell + CD4T + CD8T + Mono + NK, data = data_cvd_yes)
Horvath_cont_CVDno <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                         Bcell + CD4T + CD8T + Mono + NK, data = data_cvd_no)

# EN 
EN_cont_CVDyes <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                      Bcell + CD4T + CD8T + Mono + NK, 
                    data = data_cvd_yes)
EN_cont_CVDno <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                   Bcell + CD4T + CD8T + Mono + NK, 
                 data = data_cvd_no)

# BLUP
BLUP_cont_CVDyes <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                     Bcell + CD4T + CD8T + Mono + NK, 
                   data = data_cvd_yes)
BLUP_cont_CVDno <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                         Bcell + CD4T + CD8T + Mono + NK, 
                       data = data_cvd_no)

################################################################

#### STRATIFIED BY INCIDENT CVD - FOREST PLOT CONSTRUCTION ####

strata <- c("CVDyes", "CVDno")

## Tracking models
clocks_cont_CVDyes <- sapply(clocks, function(clock) paste0(clock, "_cont_CVDyes"))
clocks_cont_CVDno  <- sapply(clocks, function(clock) paste0(clock, "_cont_CVDno"))

## For loop for each stratified group (stratum)
for (stratum in strata) {
  model_list <- get(paste0("clocks_cont_", stratum))
  # Empty vectors to store data later
  beta <- c()
  se <- c()
  pval <- c()
  ci_low <- c()
  ci_high <- c()
  # Store results for each model
  for (c in 1:length(model_list)) {
    mod <- get(model_list[c])
    beta[c] <- summary(mod)$coefficients["predimed_score", "Estimate"]
    se[c] <- summary(mod)$coefficients["predimed_score", "Std. Error"]
    pval[c] <- summary(mod)$coefficients["predimed_score", "Pr(>|t|)"]
    ci_low[c] <- beta[c] - 1.96 * se[c]
    ci_high[c] <- beta[c] + 1.96 * se[c]
  }
  # Join results in a data frame
  df <- data.frame(clock = clocks, beta, se, pval, ci_low, ci_high) %>%
    mutate(sig = case_when(pval < 0.001 ~ "***", pval < 0.01  ~ "**", pval < 0.05  ~ "*", TRUE ~ "ns"), # add significance label
           direction = ifelse(beta < 0, "Lower AccAge (protective)", "Higher AccAge (detrimental)"), # add direction label
           group = stratum) # add group
  # Name results after stratum 
  assign(paste0("results_cont_", stratum), df)
}

## Combine overall results + stratified results
results_dodged <- bind_rows(results_cont %>% mutate(group = "Overall"),
                            results_cont_CVDyes %>% mutate(group = "CVD: Yes"),
                            results_cont_CVDno  %>% mutate(group = "CVD: No")) %>%
  mutate(clock = factor(clock, levels = rev(clocks)), # rev() so Levine plots at top
         group = factor(group, levels = c("Overall", "CVD: Yes", "CVD: No")))

## Adding a grey band per clock (to make it more understandable)
band_data <- data.frame(clock = levels(results_dodged$clock)) %>%
  mutate(y = as.numeric(factor(clock, levels = levels(results_dodged$clock))),shade = y %% 2 == 0) %>% filter(shade)

## Different colour per stratum 
group_colors <- c("Overall"  = "#2166ac", "CVD: Yes" = "#d73027", "CVD: No"  = "#1a9850")

## Constructing forest plot (stratified)
pd <- position_dodge(width = 0.55) # dodge width controls vertical spacing
forest_plot_dodged <- ggplot(results_dodged, aes(x = beta, y = clock, color = group, group = group)) +
  geom_rect(data = band_data, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf), fill = "grey92", color = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7) + 
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.75, position = pd) +
  geom_point(size = 3, shape = 16, position = pd) +
  scale_color_manual(values = group_colors, name = "Stratum") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Association of High Mediterranean Diet Adherence\nwith Epigenetic Age Acceleration (Methylation clocks)",
       subtitle = "Linear regression: AccAge ~ predimed_score (continuous) + covariates | Stratified by incident CVD",
       x = "Beta Coefficient (years)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(15, 20, 15, 15))

forest_plot_dodged

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
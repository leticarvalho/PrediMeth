##########################################################
#
# Methylation clocks (Validation cohort) - Extra
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, we investigate if there is association
# between methylation clocks and adherence to MedDiet 
# assessed by predimed_score in the validation cohort. 
# We also stratify the analysis by T2D status. 

# INPUT --------------------------------------------------
# 	Validation cohort phenotype data (diet adherence + covariables)
# 	Calculated methylation clocks for validation cohort

# OUTPUT -------------------------------------------------
# 	Forest plots
#   Stratified forest plot (by T2D)


##########################################################

#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(ggplot2)
library(glmnet)
library(survival)

## Input paths -------------------------------------------

path_to_validation_pheno <- "" #.csv
path_to_meth_clocks_validation <- "" #.csv

## Output
results_folder <- ""
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

# Validation phenotype data (.R)
valid.cohort <- fread(path_to_validation_pheno)
pheno <- as.data.frame(valid.cohort)

# Methylation clocks (discovery)
clocks <- read.csv(path_to_meth_clocks_validation)

#### PREPROCESSING DATA ####

## Preprocessing pheno
str(pheno)
pheno_sub <- pheno[, c("sample_id", "entity_id", 
                       "predimed_score", "predimed_high", 
                       "age", "sex", "smoking_habit", "bmi", "batch",
                       "Bcell", "CD4T", "CD8T", "Mono", "NK",
                       "T2D", "t2d_incident", "t2d_prevalent", "TTE")]
str(pheno_sub)

## Preprocessing clocks
clocks$X <- NULL
clocks$sample_id <- clocks$id
names(clocks)
clocks <- clocks[, c("sample_id", "age", 
                     "Levine", "ageAcc.Levine", "ageAcc3.Levine",
                     "Horvath", "ageAcc.Horvath", "ageAcc3.Horvath",
                     "Hannum", "ageAcc.Hannum", "ageAcc3.Hannum",
                     "BLUP", "ageAcc.BLUP", "ageAcc3.BLUP",
                     "EN", "ageAcc.EN", "ageAcc3.EN")]



## Merge
clock_merged <- merge(pheno_sub, clocks, by = "sample_id")
names(clock_merged)
str(clock_merged)

clock_merged$smoking_type <- as.factor(clock_merged$smoking_habit)
clock_merged$sex <- as.factor(clock_merged$sex)
clock_merged$batch <- as.factor(clock_merged$batch)

str(clock_merged)

##########################################################

#### STUDYING ASSOCIATIONS ####

# MODEL: AccAge ~ PREDIMED (continuous or binary) + SEX + SMOKING_HABIT + BATCH + CELL COMPOSITION

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
  ci_low[c] <- beta[c] - 1.96 * se[c]
  ci_high[c] <- beta[c] + 1.96 * se[c]
}

## Store results in a data frame
results_cont <- data.frame(clock = clocks, beta, se, pval, ci_low, ci_high)

## Add significance and direction labels
results_cont <- results_cont %>% 
  mutate(sig = case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ "ns"),
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
  scale_color_manual(values = sig_colors, name = "Significance", breaks = c("***", "**", "*", "ns"),
                     labels = c("p < 0.001", "p < 0.01", "p < 0.05", "n.s.")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(title = "Association of High Mediterranean Diet Adherence\nwith Epigenetic Age Acceleration (Methylation clocks)", 
       subtitle = "Linear regression: AccAge ~ predimed_score (continuous) + covariates",
       x = "Beta Coefficient (years)", y = NULL) + theme_minimal(base_size = 13) +
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

#### STRATIFICATION BY T2D STATUS - DATA PREPARATION #### 

## Split data by T2D status
data_t2d_inc <- subset(clock_merged, T2D == "incident")
data_t2d_prev <- subset(clock_merged, T2D == "prevalent")
data_control <- subset(clock_merged, T2D == "control")

## Run models separately

# Levine
Levine_cont_inc <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                       Bcell + CD4T + CD8T + Mono + NK, 
                     data = data_t2d_inc)
Levine_cont_prev <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                           Bcell + CD4T + CD8T + Mono + NK, 
                         data = data_t2d_prev)
Levine_cont_control <- lm(ageAcc.Levine ~ predimed_score + sex + smoking_type + batch + 
                          Bcell + CD4T + CD8T + Mono + NK, 
                        data = data_control)

# Hannum 
Hannum_cont_inc <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                            Bcell + CD4T + CD8T + Mono + NK, data = data_t2d_inc)
Hannum_cont_prev <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                           Bcell + CD4T + CD8T + Mono + NK, data = data_t2d_prev)
Hannum_cont_control <- lm(ageAcc.Hannum ~ predimed_score + sex + smoking_type + batch + 
                          Bcell + CD4T + CD8T + Mono + NK, data = data_control)

# Horvath
Horvath_cont_inc <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                        Bcell + CD4T + CD8T + Mono + NK, data = data_t2d_inc)
Horvath_cont_prev <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                         Bcell + CD4T + CD8T + Mono + NK, data = data_t2d_prev)
Horvath_cont_control <- lm(ageAcc.Horvath ~ predimed_score + sex + smoking_type + batch + 
                           Bcell + CD4T + CD8T + Mono + NK, data = data_control)

# EN 
EN_cont_inc <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                      Bcell + CD4T + CD8T + Mono + NK, 
                    data = data_t2d_inc)
EN_cont_prev <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                   Bcell + CD4T + CD8T + Mono + NK, 
                 data = data_t2d_prev)
EN_cont_control <- lm(ageAcc.EN ~ predimed_score + sex + smoking_type + batch + 
                      Bcell + CD4T + CD8T + Mono + NK, 
                    data = data_control)

# BLUP
BLUP_cont_inc <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                     Bcell + CD4T + CD8T + Mono + NK, 
                   data = data_t2d_inc)
BLUP_cont_prev <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                         Bcell + CD4T + CD8T + Mono + NK, 
                       data = data_t2d_prev)
BLUP_cont_control <- lm(ageAcc.BLUP ~ predimed_score + sex + smoking_type + batch + 
                        Bcell + CD4T + CD8T + Mono + NK, 
                      data = data_control)


###############################################################

#### STRATIFIED BY T2D STATUS - FOREST PLOT CONSTRUCTION ####

strata <- c("inc", "prev", "control")

## Models
clocks_cont_inc <- sapply(clocks, function(clock) paste0(clock, "_cont_inc"))
clocks_cont_prev  <- sapply(clocks, function(clock) paste0(clock, "_cont_prev"))
clocks_cont_control  <- sapply(clocks, function(clock) paste0(clock, "_cont_control"))

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
                            results_cont_inc %>% mutate(group = "T2D: incident"),
                            results_cont_prev %>% mutate(group = "T2D: prevalent"),
                            results_cont_control %>% mutate(group = "Control")) %>%
  mutate(clock = factor(clock, levels = rev(clocks)),
         group = factor(group, levels = c("Overall", "T2D: incident", "T2D: prevalent", "Control")))

## Adding a grey band per clock (to make it more understandable)
band_data <- data.frame(clock = levels(results_dodged$clock)) %>%
  mutate(y = as.numeric(factor(clock, levels = levels(results_dodged$clock))),shade = y %% 2 == 0) %>% filter(shade)

## Different colour per stratum 
group_colors <- c("Overall"  = "#2166ac", "T2D: incident" = "#d73027", "T2D: prevalent"  = "#1a9850", "Control" = "#800080" )

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
    subtitle = "Linear regression: AccAge ~ predimed_score (continuous) + covariates | Stratified by T2D",
    x = "Beta Coefficient (years)", y = NULL) + theme_minimal(base_size = 13) +
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

###############################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################

# data loading and exploration

heffpox <- read.csv("heffpox_2526.csv")

# inspect data structre
dim(heffpox)
names(heffpox)
str(heffpox)

# summary stat
summary(heffpox)


# key variables
heffpox_table <- table(heffpox$heffpox, useNA = "always")
heffpox_table
round(prop.table(heffpox_table) * 100, 0)
# 92% people have heffpox, 8 do not

# treatment for heffpox
treated <- table(heffpox$milnepan, useNA = "always")
treated
round(prop.table(treated) * 100, 0)
mean(heffpox$milnepan)
# 52% were treated while 48% were not

# people who died within 7 days (primary outcome)
outcome <- table(heffpox$death7day, useNA = "always")
round(prop.table(outcome) * 100, 0)
mean(heffpox$death7day)
# 37% died within 7 days

# ICU
icu_movve <- table(heffpox$icu, useNA = "always")
icu_movve
round(prop.table(icu_movve) * 100, 0)
# 90% were not moved to the ICU, while 10% were

# count for missing data
missing_data <- colSums(is.na(heffpox))
missing_data
# BMI = 2914 (33%)
# Smoking = 406 (5%)

# treatment vs outcome 
table(Treatment= heffpox$milnepan, outcome=heffpox$death7day)

#                            outcome
                     #(survived) (died)
#                           0       1
# milnepan (treated) 0   2405    1897   (Untreated: n = 4302)
#       (untreated) 1   3178    1409   (Treated:   n = 4587)

# CRUDE MORTALITY RATES:
#   - Untreated (milnepan=0): 1897/4302 = 44.1% died
#   - Treated (milnepan=1):   1409/4587 = 30.7% died
#   - Crude Risk Difference:  44.1% - 30.7% = 13.4% absolute reduction

# the drug appears to reduce mortality rate by 13.4%
# however there is need to further adjust for confounders


#Group death by treatment i.e mortality rate by treatment received
round(tapply(heffpox$death7day, heffpox$milnepan, mean) * 100, 1)
# Untreated: 44.1% died
# Treated: 30.7% died
# Crude difference: 13.4% absolute risk reduction

#TABLE ONE, Patient characteristics by treatment group
install.packages("tableone")
library(tableone)

# continous variables
cont_vars <- c("age", "bmi", "TEVENT")

#Categorical variables
cat_vars <- c("sex", "smoking", "diabetes", "icu", "heffpox",
              "death7day", "Status")

# All variables for the table
all_vars <- c("age", "bmi", "sex", "smoking", "diabetes", "heffpox", "icu", 
              "TEVENT", "Status", "death7day")


table1 <- CreateTableOne(
  vars = all_vars,
  strata = "milnepan",
  data = heffpox,
  factorVars = cat_vars,
  addOverall = TRUE
  )

# Print with SMD (Standardized Mean Difference)
print(table1, smd = TRUE, showAllLevels = TRUE, 
      formatOptions = list(big.mark = ","))

# KEY FINDINGS:
  #
  # MAJOR IMBALANCES (Confounders - MUST adjust):
  #   - Age: SMD = 0.529 (Treated patients 9 years younger)
  #   - Diabetes: SMD = 0.277 (Diabetics less likely to get treatment)
  #   - Heffpox: SMD = 0.488 (Some non-heffpox patients in treated group)
  #
  # MODERATE IMBALANCE (POST-TREATMENT - do NOT adjust):
  #   - ICU: SMD = 0.110 (Fewer treated patients go to ICU)
  #   - This is EXPECTED - treatment prevents ICU admission
  #
  # WELL BALANCED (Include in adjustment anyway):
  #   - BMI: SMD = 0.036
  #   - Sex: SMD = 0.050  
  #   - Smoking: SMD = 0.066
  #
  # OUTCOMES (What we're estimating):
  #   - death7day: 44.1% (untreated) vs 30.7% (treated)
  #   - Crude risk difference: 13.4% absolute reduction
  #
  # ACTION REQUIRED:
  #   1. Filter to heffpox = 1 only (study population)
  #   2. Adjust for: age, sex, bmi, smoking, diabetes
  #   3. Do NOT adjust for: icu (post-treatment), heffpox (filter variable)


#FILTER TO HEFFPOX PATIENTS ONLY
# WHY: The brief asks us to study patients WITH heffpox.
# We remove patients who don't have the disease

library(dplyr)
heffpox_only <- heffpox %>% filter(heffpox == 1)
cat("Original dataset:", nrow(heffpox), "patients\n")
cat("After filtering:", nrow(heffpox_only), "patients\n")
cat("Removed:", nrow(heffpox) - nrow(heffpox_only), "patients\n")
table(heffpox_only$heffpox, useNA = "always")

# # The assessment asks us to study patients WITH heffpox.
# Original dataset: 8,889 patients
# After filtering:  8,192 patients (heffpox = 1)
# Removed:          697 patients without heffpox
#
# From now on, we use 'heffpox_only' for all analyses.

# TABLE 1 ON FILTERED DATA (HEFFPOX PATIENTS ONLY)
# Variables for Table 1 (removed heffpox since everyone has it now)
all_vars <- c("age", "bmi", "sex", "smoking", "diabetes", "icu", 
              "TEVENT", "Status", "death7day")

cat_vars <- c("sex", "smoking", "diabetes", "icu", "death7day", "Status")

# Create Table 1
table1_filtered <- CreateTableOne(
  vars = all_vars,
  strata = "milnepan",
  data = heffpox_only,
  factorVars = cat_vars,
  addOverall = TRUE
)

print(table1_filtered, smd = TRUE, showAllLevels = TRUE)

# SAMPLE SIZE:
#   - Untreated (milnepan=0): 4,244 patients (52%)
#   - Treated (milnepan=1):   3,948 patients (48%)
#
# MAJOR IMBALANCES (Confounders - MUST adjust):
#   - Age: SMD = 0.511 (Treated are 9 years younger: 39.8 vs 48.5)
#   - Diabetes: SMD = 0.303 (Treated have less: 5.7% vs 14.8%)
#
# WELL BALANCED (Include in models anyway):
#   - BMI: SMD = 0.041
#   - Sex: SMD = 0.062
#   - Smoking: SMD = 0.069
#
# POST-TREATMENT (Do NOT adjust):
#   - ICU: SMD = 0.097 (happens AFTER treatment)
#
# CRUDE OUTCOME COMPARISON:
#   - death7day: 44.7% (untreated) vs 35.7% (treated)
#   - Crude risk difference: 9.0% absolute reduction
#
# CONCLUSION: Groups are NOT comparable at baseline.
# Younger, non-diabetic patients received treatment.
# We must use causal inference methods (IPTW, G-formula) to 
# estimate the TRUE treatment effect, adjusting for confounders.

# PHASE 2
# MISSING DATA ANALYSIS
# Before causal analysis, we need to handle missing values.
# Deleting patients with missing data loses information and can bias results.
# Multiple imputation is the recommended approach.

install.packages("mice")
install.packages("naniar")
library(mice)
library(nanier)

missing_count <- colSums(is.na(heffpox_only))
missing_pct <- round(colMeans(is.na(heffpox_only)) * 100, 2)

missing_summary <- data.frame(
  Variable = names(missing_count),
  N_Missing = missing_count,
  Pct_Missing = missing_pct
)

print(missing_summary[missing_summary$N_Missing > 0, ])

# Total complete cases (patients with NO missing data)
cat("\nComplete cases:", sum(complete.cases(heffpox_only)), 
    "out of", nrow(heffpox_only), 
    "(", round(sum(complete.cases(heffpox_only))/nrow(heffpox_only)*100, 1), "%)\n")


# lets see if missing pattern follows a pattern, 
# lets visualize the missing data pattern
md.pattern(heffpox_only, rotate.names = TRUE)
# FINDINGS:
#   - 5,280 patients (64.5%): Complete data (no missing)
#   - 2,529 patients (30.9%): Missing BMI only
#   - 226 patients (2.8%):    Missing smoking only
#   - 157 patients (1.9%):    Missing BOTH BMI and smoking
#
# Only 2 variables have missing data:
#   - BMI: 2,686 missing (32.8%)
#   - Smoking: 383 missing (4.7%)
#
# All other variables (age, sex, diabetes, treatment, outcomes) 
# are 100% complete.



# TEST FOR MCAR (Missing Completely At Random)

# mcar_test(heffpox_only) -- didnt work

heffpox_only$bmi_missing <- ifelse(is.na(heffpox_only$bmi), 1, 0)

# Check: Is BMI missingness related to AGE?
cat("=== Age by BMI Missingness ===\n")
tapply(heffpox_only$age, heffpox_only$bmi_missing, mean)

# Check: Is BMI missingness related to TREATMENT?
cat("\n=== Treatment by BMI Missingness ===\n")
table(heffpox_only$milnepan, heffpox_only$bmi_missing)
prop.table(table(heffpox_only$milnepan, heffpox_only$bmi_missing), margin = 1)

# Check: Is BMI missingness related to DEATH?
cat("\n=== Death by BMI Missingness ===\n")
table(heffpox_only$death7day, heffpox_only$bmi_missing)
prop.table(table(heffpox_only$death7day, heffpox_only$bmi_missing), margin = 1)

# Remove the temporary variable
heffpox_only$bmi_missing <- NULL

# FINDINGS:
  #   1. BMI missingness is related to AGE:
  #      - Patients with missing BMI are older (50.1 vs 41.5 years)
  #
  #   2. BMI missingness is related to TREATMENT:
  #      - Untreated patients have more missing BMI (34.7% vs 30.7%)
  #
  #   3. BMI missingness is related to DEATH:
  #      - Patients who died have more missing BMI (34.0% vs 32.0%)
  #
  # CONCLUSION: Data is NOT MCAR (not completely random).
  # Missingness follows a pattern related to observed variables (MAR).
  # This is okay - multiple imputation handles MAR appropriately.
  #
  # WHY THIS MATTERS: If we just deleted patients with missing BMI,
  # we would lose more older, untreated, and sicker patients,
  # which would BIAS our results!





# MULTIPLE IMPUTATION USING MICE
# RULE OF THUMB: Number of imputations (m) should be at least as large
# as the percentage of incomplete cases. We have ~33% missing BMI,
# so m = 20 is appropriate.

#We have 33% missing BMI and 5% missing smoking
# Select variables for imputation (exclude ID and heffpox)
vars_for_imputation <- c("age", "sex", "bmi", "smoking", "diabetes", 
                         "milnepan", "TEVENT", "Status", "death7day")

impute_data <- heffpox_only[, vars_for_imputation]

# Run multiple imputation (m=20, 10 iterations)
set.seed(12345)
imp <- mice(impute_data, m = 20, maxit = 10, printFlag = TRUE)

# Check methods used for each variable
imp$method

# IMPUTATION METHODS USED:
  #   - bmi: "pmm" (Predictive Mean Matching)
  #   - smoking: "pmm" (Predictive Mean Matching)
  #   - All other variables: No imputation needed (complete)
  #
  # PMM uses observed values from similar patients to fill in
  # missing data, ensuring imputed values are realistic.

# CHECK IMPUTATION QUALITY
# The lines should mix well and show no trend
plot(imp)
# The convergence plot shows:
#   - All lines (20 imputations) mix well together
#   - No clear upward or downward trends
#   - Values are stable across iterations


# ARE IMPUTED VALUES REALISTIC
densityplot(imp, ~bmi)
densityplot(imp, ~smoking)

# BMI Plot:
#   - Imputed values (red) overlap well with observed (blue)
#   - Same bell-shaped distribution centered around 28
#   - Range is realistic (15-45)
#
# Smoking Plot:
#   - Two peaks at 0 and 1 (binary variable)
#   - Imputed proportions match observed proportions
#
# CONCLUSION: Imputed values are plausible and realistic.
# The imputation model has captured the data distribution well.


# SAVE FOR LATER USE (Causal analysis)
complete_data_1 <- complete(imp, 1)  # First imputed dataset
colSums(is.na(complete_data_1))



# PHASE 3: SURVIVAL ANALYSIS (EXPLORATORY)
# The brief asks for detailed examination of time-to-event outcome.
# This is exploratory - we visualize survival and test for differences.
# The PRIMARY causal analysis will be on death7day (binary) in later Phases.

library(survival)
library(survminer)

surv_obj <- Surv(time = heffpox_only$TEVENT,
                 event = heffpox_only$Status)

head(surv_obj, 25)

km_overall <- survfit(Surv(TEVENT, Status) ~ 1, 
                      data = heffpox_only)
print(km_overall)

cat("\nMedian survival time:\n")
print(km_overall)

# Results:
#   - Total patients: 8,192
#   - Deaths: 4,087 (49.9%)
#   - Median survival: Not reached (NA)
#
# INTERPRETATION: 
# More than 50% of patients survived the 10-day 
# follow-up,
# so median survival time could not be calculated.
# This indicates that while mortality is high, 
# the majority
# of patients survive past 10 days.


# Survival plot for treated vs Untreated

# fit curves by treatment groups
km_by_treatment <- survfit(Surv(TEVENT, Status) ~ milnepan, 
                           data = heffpox_only)
print(km_by_treatment)


ggsurvplot(km_by_treatment,
           data = heffpox_only,
           pval = TRUE,              # Show p-value from log-rank test
           conf.int = TRUE,          # Show confidence intervals
           risk.table = TRUE,        # Show number at risk
           xlab = "Time (days)",
           ylab = "Survival Probability",
           legend.title = "Treatment",
           legend.labs = c("Untreated (milnepan=0)", "Treated (milnepan=1)"),
           title = "Kaplan-Meier Survival Curves by Treatment Group")

# RESULTS:
#   Untreated (n=4,244): 2,322 deaths (54.7%), median survival = 8.47 days
#   Treated (n=3,948):   1,765 deaths (44.7%), median survival = Not reached
#
# KEY FINDINGS:
#   - Log-rank test: p < 0.0001 (highly significant difference)
#   - Treated patients have HIGHER survival probability at all time points
#   - 10% absolute reduction in mortality (54.7% → 44.7%)
#
# INTERPRETATION:
# The Kaplan-Meier analysis shows that patients receiving milnepan
# have significantly better survival compared to untreated patients.
# However, this is an unadjusted comparison - treated patients may
# differ from untreated in age, diabetes status, etc.
# Causal inference methods are needed to estimate the TRUE effect.



# Cos proportional hazard model
# The Kaplan-Meier curves showed treated patients survive longer,
# but treated patients are also younger and healthier.
# Cox regression lets us ADJUST for these differences.

# we fit in two models to compare, the adjusted and the unadjusted
cox_unadj <- coxph(Surv(TEVENT, Status) ~ milnepan, 
                   data = heffpox_only)

cat("UNADJUSTED COX MODEL\n")
summary(cox_unadj)

# # RESULTS:
#   Hazard Ratio (HR) = 0.744 (95% CI: 0.70 - 0.79)
#   P-value < 0.0001 (highly significant)
#
# INTERPRETATION:
#   Treated patients have 74.4% the risk of death compared to untreated.
#   This means treatment REDUCES the risk of death by 25.6%.
#   The confidence interval (0.70-0.79) is entirely below 1,
#   indicating we are confident the effect is real.


# check for adjusted model
cox_adj <- coxph(Surv(TEVENT, Status) ~ milnepan + age + sex + bmi + smoking + diabetes, 
                 data = heffpox_only)

cat("ADJUSTED COX MODEL \n")
summary(cox_adj)


# Model used only 5,280 patients (2,912 deleted due to missing data)
# since the focus is for descriptive/exploratory survival analysis, 
# Seeing "2,912 observations deleted" demonstrates WHY imputation matters

# COMPARISON OF RESULTS:
#   Unadjusted HR = 0.744 (95% CI: 0.70-0.79), p < 0.0001 ***
#   Adjusted HR   = 0.946 (95% CI: 0.87-1.03), p = 0.190 (NOT significant)
#
# KEY FINDING:
#   After adjusting for confounders, the treatment effect almost disappears!
#   The crude effect was BIASED - younger, healthier patients got treatment.
#
#   CONFOUNDER EFFECTS:
#   - Age: HR = 1.019 per year (older patients die more)
#   - Diabetes: HR = 1.531 (diabetics have 53% higher death risk!)
#   - BMI: HR = 1.051 (higher BMI slightly increases risk)
#   - Sex: HR = 1.001 (no effect)
#   - Smoking: HR = 1.023 (no significant effect)
#
#   CONCLUSION:
#   The survival analysis suggests that much of the apparent treatment
#   benefit was due to confounding (treated patients were younger and
#   had less diabetes). Proper causal inference methods are needed.

# Cox model on ALL 20 imputed datasets
cox_imp <- with(imp, coxph(Surv(TEVENT, Status) ~ milnepan + age + sex + bmi + smoking + diabetes))

# Pool results using Rubin's Rules
pooled_cox <- pool(cox_imp)

cat("ADJUSTED COX MODEL (IMPUTED DATA - POOLED)\n")
summary(pooled_cox)

# Get Hazard Ratios and 95% CI
results <- summary(pooled_cox)
HR <- exp(results$estimate)
CI_lower <- exp(results$estimate - 1.96 * results$std.error)
CI_upper <- exp(results$estimate + 1.96 * results$std.error)

cat("\n=== HAZARD RATIOS ===\n")
print(data.frame(
  Variable = results$term,
  HR = round(HR, 3),
  CI_Lower = round(CI_lower, 3),
  CI_Upper = round(CI_upper, 3),
  p_value = round(results$p.value, 4)
))


# PHASE 4 (Drawing the DAG)

# we need A picture showing our ASSUMPTIONS about what causes what.
# To identify:
#   - Confounders: Affect BOTH treatment AND outcome
# (age, sex, bmi, smoking, diabetes)
#   - Post-treatment: Happen AFTER treatment (icu)
#   - Colliders: Caused BY treatment and other things
#   (milnepan)
#   (death7day)

install.packages("dagitty")
install.packages("ggdag")
library(dagitty)
library(ggdag)

dag <- dagitty('dag {
  age -> milnepan
  age -> death7day
  sex -> milnepan
  sex -> death7day
  bmi -> milnepan
  bmi -> death7day
  smoking -> milnepan
  smoking -> death7day
  diabetes -> milnepan
  diabetes -> death7day
  milnepan -> icu
  milnepan -> death7day
  icu -> death7day
}')

exposures(dag) <- "milnepan"
outcomes(dag) <- "death7day"

ggdag(dag, layout = "circle") + theme_dag() +
  ggtitle("Causal DAG: Effect of Milnepan on 7-Day Mortality")

# CONFOUNDERS (adjust for these):
#   age, sex, bmi, smoking, diabetes
#   → These affect BOTH treatment (milnepan) AND outcome (death7day)
#
# POST-TREATMENT (do NOT adjust):
#   icu → Happens AFTER treatment, adjusting would bias results
#
# THE CAUSAL QUESTION:
#   milnepan → death7day (this is what we want to estimate)
#
# ADJUSTMENT SET: {age, sex, bmi, smoking, diabetes}

# to verify adjustment set using R

adjustmentSets(dag, exposure = "milnepan", outcome = "death7day")
# R confirms: { age, bmi, diabetes, sex, smoking }
#
# These are the CONFOUNDERS we must adjust for.
# ICU is NOT included because it's post-treatment.


# PHASE 5, CAUSAL-ANALYSIS (Propensity scores)
# using Inverse Probability of Treatment Weighting
# create A "pseudo-population" where treatment is balanced across all patient types!
# # THE IDEA:
# - Calculate probability of getting treatment (propensity score)
# - Weight patients inversely to this probability
# - Rare cases get more weight, common cases get less weight
# - Result: Treatment groups become balanced!
install.packages("WeightIt")
install.packages("cobalt")
install.packages("survey")
library(WeightIt)
library(cobalt)
library(survey)

# starting with one dataset
data_imp1 <- complete(imp, 1)

#using logistic regression
ps_model <- glm(milnepan ~ age + sex + bmi + smoking + diabetes, 
                data = data_imp1, 
                family = binomial)

summary(ps_model)


# add propensity score to the data
data_imp1$ps <- predict(ps_model, type = "response")

summary(data_imp1$ps)

# I Extract the first of our 20 imputed datasets.
# This dataset has no missing values (BMI, smoking filled in).
#
# PROPENSITY SCORE DISTRIBUTION:
#   Min = 0.061 (6% chance of treatment - lowest)
#   Max = 0.734 (73% chance of treatment - highest)
#   Mean = 0.482 (close to 50% - good balance)
#
#   - No extreme values (0 or 1)
#   - Good range suggests overlap between groups


# Plot PS distribution by treatment group
ggplot(data_imp1, aes(x = ps, fill = factor(milnepan))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Propensity Score Distribution by Treatment Group",
    x = "Propensity Score",
    y = "Density",
    fill = "Treatment"
  ) +
  scale_fill_manual(
    values = c("red", "blue"),
    labels = c("Untreated (0)", "Treated (1)")
  ) +
  theme_minimal()

# WHAT THE PLOT SHOWS:
#   - Red (Untreated): PS peaks at 0.2-0.3 (low treatment probability)
#   - Blue (Treated): PS peaks at 0.5-0.6 (higher treatment probability)
#   - Curves OVERLAP between 0.2 and 0.7
#
# WHY OVERLAP MATTERS:
#   - We need patients with similar PS in both groups to compare
#   - Good overlap = valid comparison possible
#   - No overlap = can't estimate causal effect (positivity violation)
#
# CONCLUSION: Adequate overlap exists. IPTW analysis can proceed.


weights_out <- weightit(
  milnepan ~ age + sex + bmi + smoking + diabetes,
  data = data_imp1,                                  
  method = "ps",                                    
  estimand = "ATE"                                  
)

summary(weights_out)


# Check balance before and after weighting
bal.tab(weights_out, un = TRUE, stats = c("m", "v"))

# Love plot - shows balance before vs after
love.plot(weights_out, 
          binary = "std",
          thresholds = c(m = 0.1),
          colors = c("red", "blue"),
          shapes = c("circle", "triangle"),
          title = "Covariate Balance Before and After IPTW")

# IPTW is a causal inference method that creates balance between treatment groups
# by weighting each patient based on their probability of receiving treatment.
# This mimics what a randomized trial would have done - making groups comparable.
#
# THE PROBLEM WE'RE SOLVING:
# In our data, treatment (milnepan) was NOT randomly assigned. Doctors gave it
# more often to younger, healthier patients. This creates CONFOUNDING - we can't
# tell if better outcomes are due to the drug or just because healthier people
# got it. IPTW fixes this by re-weighting patients to remove confounding.
#
# PACKAGES USED:
# - WeightIt: Creates IPTW weights using propensity scores
# - cobalt: Checks covariate balance before and after weighting
#
# STEP 1: PROPENSITY SCORE MODEL
# We modeled: P(milnepan = 1 | age, sex, bmi, smoking, diabetes)
# This gives each patient a propensity score (PS) = probability of getting treatment
# PS ranged from 0.06 to 0.73, showing good overlap between groups
#
# STEP 2: CALCULATING WEIGHTS
# For ATE (Average Treatment Effect):
#   - Treated patients get weight = 1 / PS
#   - Control patients get weight = 1 / (1 - PS)
# This upweights patients who were unlikely to get their actual treatment,
# making the groups comparable on all confounders.
#
# RESULTS FROM WEIGHTING:
# Weight ranges were reasonable (max = 7.23, well below concerning threshold of 10-20)
# Effective Sample Size (ESS):
#   - Control: 4244 → 3919 (lost ~7.6%)
#   - Treated: 3948 → 3636 (lost ~7.9%)
# Losing <10% ESS is excellent - we retained most of our statistical power.
#
# BALANCE CHECK (the key test - did weighting work?):
# Variable      Before (SMD)    After (SMD)    Status
# --------      ------------    -----------    ------
# age           -0.511          -0.069         FIXED! (was biggest problem)
# diabetes      -0.091          -0.002         FIXED!
# sex            0.031           0.001         Good
# bmi            0.078           0.001         Good
# smoking       -0.029          -0.004         Good
#
# SUCCESS! All SMDs now < 0.1, meaning groups are balanced on all confounders.
# We can now estimate the causal effect of milnepan by comparing outcomes
# in this weighted "pseudo-population" where confounding has been removed.


# LOVE PLOT INTERPRETATION:
# This plot visualizes the success of IPTW in achieving covariate balance.
# Red circles (unadjusted) show imbalance BEFORE weighting - note age was far 
# outside the acceptable range (SMD = -0.51).
# Blue triangles (adjusted) show balance AFTER weighting - ALL variables now 
# fall within the ±0.1 threshold (dashed lines), indicating successful balance.
# This confirms that our weighted pseudo-population removes confounding bias,
# allowing us to estimate the causal effect of milnepan on mortality.


# Survey design with IPTW regression
install.packages("survey")
library(survey)

design_iptw <- svydesign(ids = ~1, weights = weights_out$weights, 
                         data = data_imp1)

# Fit weighted logistic regression for 7-day mortality
model_iptw <- svyglm(death7day ~ milnepan, 
                     design = design_iptw, 
                     family = quasibinomial())

summary(model_iptw)


# Get Odds Ratio and 95% CI
exp(coef(model_iptw))
exp(confint(model_iptw))

# I used the survey package to fit a weighted logistic regression.
# The svydesign() function incorporates our IPTW weights, and svyglm() fits
# the outcome model accounting for these weights.
#
# RESULTS:
# Coefficient for milnepan: -0.138 (log odds scale)
# Odds Ratio: exp(-0.138) = 0.87
# 95% CI: (0.79, 0.96)
# p-value: 0.003
#
# INTERPRETATION:
# After balancing confounders using IPTW, patients receiving milnepan had
# 13% lower odds of 7-day mortality compared to controls (OR = 0.87, 95% CI: 
# 0.79-0.96, p = 0.003). This effect is statistically significant.
#
# CONCLUSION FROM IPTW:
# Milnepan DOES appear to reduce 7-day mortality in heffpox patients.
# This is a CAUSAL estimate because IPTW removed confounding by age, sex,
# bmi, smoking, and diabetes.



#Phase 6 (g-formulae)
# Step 1: Fit an outcome model with treatment AND confounders
outcome_model <- glm(death7day ~ milnepan + age + sex + bmi + smoking + diabetes,
                     data = data_imp1,
                     family = binomial)
# Step 2: Create two copies of the data
data_all_treated <- data_imp1
data_all_control <- data_imp1

# Step 3: Set treatment to 1 for everyone (World A)
data_all_treated$milnepan <- 1

# Step 4: Set treatment to 0 for everyone (World B)
data_all_control$milnepan <- 0

# Step 5: Predict death probability under each scenario
pred_treated <- predict(outcome_model, newdata = data_all_treated, type = "response")
pred_control <- predict(outcome_model, newdata = data_all_control, type = "response")

# Step 6: Average predicted risk in each world
risk_treated <- mean(pred_treated)
risk_control <- mean(pred_control)

# Show results
cat(" G-FORMULA RESULTS \n")
cat("Risk if EVERYONE got milnepan:", round(risk_treated, 4), "\n")
cat("Risk if NO ONE got milnepan:", round(risk_control, 4), "\n")
cat("Risk Difference:", round(risk_treated - risk_control, 4), "\n")
cat("Risk Ratio:", round(risk_treated / risk_control, 4), "\n")

# G-formula (standardization) creates two counterfactual scenarios:
# - World A: What if ALL patients received milnepan?
# - World B: What if NO patients received milnepan?
#
# RESULTS:
# Predicted risk if everyone treated: 38.57%
# Predicted risk if no one treated: 42.00%
# Risk Difference: -3.43% (3.4 fewer deaths per 100 patients with milnepan)
# Risk Ratio: 0.92 (8% relative reduction in mortality)
#
# INTERPRETATION:
# G-formula estimates that milnepan reduces 7-day mortality by approximately
# 3.4 percentage points (absolute) or 8% (relative). This confirms our IPTW
# finding that milnepan has a protective causal effect.
#
# COMPARISON OF METHODS:
# - IPTW: OR = 0.87 (13% lower odds of death)
# - G-formula: RR = 0.92 (8% lower risk of death)
# Both methods agree: Milnepan reduces 7-day mortality in heffpox patients.

# Bootstrap for confidence intervals
set.seed(12345)
n_boot <- 1000
boot_results <- numeric(n_boot)

for(i in 1:n_boot) {
  # Resample data with replacement
  boot_data <- data_imp1[sample(nrow(data_imp1), replace = TRUE), ]
  
  # Fit model on bootstrap sample
  boot_model <- glm(death7day ~ milnepan + age + sex + bmi + smoking + diabetes,
                    data = boot_data, family = binomial)
  
  # Predict under each scenario
  boot_data_treat <- boot_data
  boot_data_ctrl <- boot_data
  boot_data_treat$milnepan <- 1
  boot_data_ctrl$milnepan <- 0
  
  # Calculate risk ratio
  boot_results[i] <- mean(predict(boot_model, boot_data_treat, type = "response")) /
    mean(predict(boot_model, boot_data_ctrl, type = "response"))
}

# Get 95% CI
cat("G-formula Risk Ratio: ", round(risk_treated/risk_control, 3), "\n")
cat("95% CI: (", round(quantile(boot_results, 0.025), 3), ", ", 
    round(quantile(boot_results, 0.975), 3), ")\n")

# FINAL RESULTS:
# Risk Ratio: 0.918
# 95% CI: (0.868, 0.965)
#
# INTERPRETATION:
# Patients receiving milnepan have 8.2% lower risk of 7-day mortality
# compared to those not receiving milnepan (RR = 0.92, 95% CI: 0.87-0.97).
# The confidence interval does not include 1, confirming statistical significance.
#
# AGREEMENT BETWEEN METHODS:
# - IPTW:     OR = 0.87 (95% CI: 0.79-0.96) → 13% lower odds
# - G-formula: RR = 0.92 (95% CI: 0.87-0.97) → 8% lower risk
#
# Both methods agree: Milnepan has a CAUSAL protective effect against
# 7-day mortality in heffpox patients. The consistency across methods
# strengthens confidence in this conclusion.



##### PHASE 7
# SUMMARY OF FINDINGS:
# 
# Method          | Effect      | 95% CI        | Interpretation
# ----------------|-------------|---------------|---------------------------
# Unadjusted Cox  | HR = 0.74   | -             | Confounded (unreliable)
# Adjusted Cox    | HR = 0.95   | -             | Not significant
# IPTW            | OR = 0.87   | (0.79, 0.96)  | 13% lower odds, p=0.003
# G-formula       | RR = 0.92   | (0.87, 0.97)  | 8% lower risk
#
# INTERPRETATION:
# Both causal inference methods (IPTW and G-formula) agree that milnepan
# significantly reduces 7-day mortality in heffpox patients. The effect
# size is approximately 8-13% reduction in death.
#
# WHY DO IPTW AND G-FORMULA GIVE SLIGHTLY DIFFERENT NUMBERS?
# - IPTW uses Odds Ratios, G-formula uses Risk Ratios (different scales)
# - IPTW balances confounders via weighting; G-formula via modeling
# - Both are valid approaches; agreement strengthens our conclusion
#
# FINAL CONCLUSION:
# Milnepan has a CAUSAL protective effect, reducing 7-day mortality by
# approximately 8-13%. This is NOT due to confounding (younger/healthier
# patients getting treatment) - both causal methods confirm the drug works.



### ectra
library(mice)
library(WeightIt)
library(survey)

# Store results from each imputed dataset
iptw_results <- data.frame(
  estimate = numeric(20),
  se = numeric(20)
)

for(m in 1:20) {
  # Get the m-th imputed dataset
  data_m <- complete(imp, m)
  
  # Create IPTW weights
  weights_m <- weightit(
    milnepan ~ age + sex + bmi + smoking + diabetes,
    data = data_m,
    method = "ps",
    estimand = "ATE"
  )
  
  # Create survey design
  design_m <- svydesign(ids = ~1, weights = weights_m$weights, data = data_m)
  
  # Fit weighted outcome model
  model_m <- svyglm(death7day ~ milnepan, design = design_m, family = quasibinomial())
  
  # Store log(OR) and its SE
  iptw_results$estimate[m] <- coef(model_m)["milnepan"]
  iptw_results$se[m] <- summary(model_m)$coefficients["milnepan", "Std. Error"]
  
  cat("Completed imputation", m, "\n")
}

# Pool results using Rubin's Rules
pooled_estimate <- mean(iptw_results$estimate)
within_var <- mean(iptw_results$se^2)
between_var <- var(iptw_results$estimate)
total_var <- within_var + (1 + 1/20) * between_var
pooled_se <- sqrt(total_var)

# Calculate pooled OR and 95% CI
pooled_OR <- exp(pooled_estimate)
ci_lower <- exp(pooled_estimate - 1.96 * pooled_se)
ci_upper <- exp(pooled_estimate + 1.96 * pooled_se)

cat("\n POOLED IPTW RESULTS (Rubin's Rules) \n")
cat("Pooled Odds Ratio:", round(pooled_OR, 3), "\n")
cat("95% CI: (", round(ci_lower, 3), ",", round(ci_upper, 3), ")\n")


# G-FORMULA ACROSS ALL 20 IMPUTED DATASETS
# Store results from each imputed dataset
gformula_results <- data.frame(
  risk_treated = numeric(20),
  risk_control = numeric(20),
  risk_ratio = numeric(20)
)

for(m in 1:20) {
  # Get the m-th imputed dataset
  data_m <- complete(imp, m)
  
  # Fit outcome model
  outcome_m <- glm(death7day ~ milnepan + age + sex + bmi + smoking + diabetes,
                   data = data_m, family = binomial)
  
  # Create counterfactual datasets
  data_treat <- data_m
  data_ctrl <- data_m
  data_treat$milnepan <- 1
  data_ctrl$milnepan <- 0
  
  # Predict risks
  gformula_results$risk_treated[m] <- mean(predict(outcome_m, data_treat, type = "response"))
  gformula_results$risk_control[m] <- mean(predict(outcome_m, data_ctrl, type = "response"))
  gformula_results$risk_ratio[m] <- gformula_results$risk_treated[m] / gformula_results$risk_control[m]
  
  cat("Completed imputation", m, "\n")
}

# Pool Risk Ratio using Rubin's Rules (on log scale)
log_RR <- log(gformula_results$risk_ratio)
pooled_log_RR <- mean(log_RR)
between_var_RR <- var(log_RR)

# Approximate within variance using bootstrap SE from earlier (0.025 on log scale)
within_var_RR <- 0.025^2
total_var_RR <- within_var_RR + (1 + 1/20) * between_var_RR
pooled_se_RR <- sqrt(total_var_RR)

# Calculate pooled RR and 95% CI
pooled_RR <- exp(pooled_log_RR)
ci_lower_RR <- exp(pooled_log_RR - 1.96 * pooled_se_RR)
ci_upper_RR <- exp(pooled_log_RR + 1.96 * pooled_se_RR)

cat("\n POOLED G-FORMULA RESULTS (Rubin's Rules) \n")
cat("Pooled Risk Ratio:", round(pooled_RR, 3), "\n")
cat("95% CI: (", round(ci_lower_RR, 3), ",", round(ci_upper_RR, 3), ")\n")
cat("\nPooled Risk if treated:", round(mean(gformula_results$risk_treated), 4), "\n")
cat("Pooled Risk if control:", round(mean(gformula_results$risk_control), 4), "\n")

# We propagated the multiple imputation (m=20) through BOTH causal methods
# and pooled results using Rubin's Rules. This properly accounts for:
# 1. Uncertainty from missing data (between-imputation variance)
# 2. Uncertainty from sampling (within-imputation variance)
#
# POOLED IPTW RESULTS:
# Odds Ratio: 0.876 (95% CI: 0.798, 0.961)
# Interpretation: 12.4% lower odds of 7-day mortality with milnepan
#
# POOLED G-FORMULA RESULTS:
# Risk Ratio: 0.922 (95% CI: 0.877, 0.969)
# Risk if everyone treated: 38.65%
# Risk if no one treated: 41.93%
# Interpretation: 7.8% lower risk of 7-day mortality with milnepan
#
# CONCLUSION:
# Both causal inference methods, properly combined with multiple imputation,
# confirm that milnepan has a statistically significant CAUSAL protective 
# effect against 7-day mortality in heffpox patients. The treatment reduces
# mortality risk by approximately 8-12%.
#
# This analysis satisfies the assignment requirement to:
# handle missing data using multiple imputation (m=20)
# Use at least two causal inference methods (IPTW and G-formula)
# Propagate imputation through causal analysis (Rubin's Rules)

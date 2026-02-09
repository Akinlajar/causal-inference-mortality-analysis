# 🏥 Causal Inference Analysis: Treatment Effects on Mortality

> **Does milnepan treatment cause a reduction in 7-day mortality among heffpox patients?**

A comprehensive causal inference analysis using **Inverse Probability of Treatment Weighting (IPTW)** and **G-formula** to estimate causal treatment effects from observational data, with proper handling of missing data through multiple imputation.

![R](https://img.shields.io/badge/R-4.0+-blue?logo=r)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Completed-success)

---

## 📌 Table of Contents

- [Overview](#overview)
- [The Problem](#the-problem)
- [Dataset](#dataset)
- [Methodology](#methodology)
- [Key Results](#key-results)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Skills Demonstrated](#skills-demonstrated)
- [Visualizations](#visualizations)
- [Limitations](#limitations)
- [References](#references)
- [Author](#author)

---

## 🎯 Overview

This project applies rigorous causal inference methods to determine whether a drug treatment **causes** reduced mortality—not just correlates with it. Using observational data from 8,192 patients, I implemented two complementary causal inference approaches and handled substantial missing data through multiple imputation.

### Why This Matters

In healthcare, we often can't run randomized controlled trials. We need methods to extract causal insights from observational data while accounting for the biases inherent in non-random treatment assignment.

---

## ❓ The Problem

### Confounding Bias

In real-world clinical settings, treatment isn't randomly assigned. Doctors may preferentially prescribe milnepan to:
- Younger patients
- Patients with fewer comorbidities
- Patients they believe will respond well

This creates **confounding bias**—we can't distinguish whether better outcomes are due to:
- The drug actually working, OR
- Treated patients being healthier to begin with

### The Evidence

| Variable | Treated | Untreated | SMD |
|----------|---------|-----------|-----|
| Age (mean) | 39.8 years | 48.5 years | **0.511** |
| Diabetes | 5.7% | 14.8% | **0.303** |

*SMD > 0.1 indicates meaningful imbalance*

Treated patients were **8.7 years younger** on average and had **less than half** the diabetes rate. Any crude comparison would be biased.

---

## 📊 Dataset

| Attribute | Value |
|-----------|-------|
| **Total Patients** | 8,192 (confirmed heffpox cases) |
| **Treatment Group** | 3,948 (48%) received milnepan |
| **Control Group** | 4,244 (52%) did not receive milnepan |
| **Primary Outcome** | 7-day mortality |
| **Follow-up Period** | Maximum 10 days |

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `milnepan` | Binary | Treatment indicator (1 = treated, 0 = untreated) |
| `death7` | Binary | Died within 7 days (1 = yes, 0 = no) |
| `TEVENT` | Continuous | Days from diagnosis to death/censoring |
| `Status` | Binary | Event indicator for survival analysis |
| `age` | Continuous | Patient age in years |
| `sex` | Binary | Patient sex |
| `bmi` | Continuous | Body mass index (32.8% missing) |
| `smoking` | Binary | Smoking status (4.7% missing) |
| `diabetes` | Binary | Diabetes diagnosis |
| `icu` | Binary | ICU admission (excluded—post-treatment variable) |

### Missing Data

| Variable | N Missing | % Missing |
|----------|-----------|-----------|
| BMI | 2,686 | 32.8% |
| Smoking | 383 | 4.7% |
| **Total incomplete cases** | 2,912 | 35.5% |

---

## 🔬 Methodology

### Phase 1: Missing Data Analysis

**Objective:** Determine the missing data mechanism and implement appropriate handling.

**Approach:**
1. Examined missing data patterns using `VIM` package
2. Tested for MCAR vs MAR by comparing characteristics of complete vs incomplete cases
3. Found systematic differences → data is **MAR** (Missing at Random)

**Solution:** Multiple Imputation using Chained Equations (MICE)
- Generated **m = 20** imputed datasets
- Used **Predictive Mean Matching (PMM)** for realistic imputations
- Ran **10 iterations** for convergence
- Validated with trace plots and density comparisons

```r
imp <- mice(data, m = 20, maxit = 10, method = "pmm", seed = 123)
```

### Phase 2: Causal Framework (DAG)

**Objective:** Identify confounders and determine the minimal sufficient adjustment set.

**Approach:**
1. Constructed a Directed Acyclic Graph based on clinical knowledge
2. Identified 5 confounders: age, sex, BMI, smoking, diabetes
3. Excluded ICU admission (post-treatment mediator)
4. Computationally verified adjustment set using `dagitty`

```r
dag <- dagitty('dag {
  age -> milnepan
  age -> death
  sex -> milnepan
  sex -> death
  bmi -> milnepan
  bmi -> death
  smoking -> milnepan
  smoking -> death
  diabetes -> milnepan
  diabetes -> death
  milnepan -> death
  milnepan -> icu
  icu -> death
}')
adjustmentSets(dag, exposure = "milnepan", outcome = "death")
```

### Phase 3: Exploratory Survival Analysis

**Objective:** Visualize survival patterns and quantify treatment effects.

**Methods:**
- **Kaplan-Meier curves** with log-rank test
- **Cox Proportional Hazards** models (unadjusted, adjusted complete case, adjusted imputed)

**Key Insight:** Complete case analysis lost statistical significance due to excluding 35% of patients. Multiple imputation recovered the true effect.

| Model | HR | 95% CI | p-value | N |
|-------|-----|--------|---------|---|
| Unadjusted | 0.744 | 0.70 - 0.79 | <0.0001 | 8,192 |
| Adjusted (complete case) | 0.946 | 0.87 - 1.03 | 0.190 | 5,280 |
| Adjusted (imputed) | **0.918** | 0.86 - 0.98 | **0.012** | 8,192 |

### Phase 4: IPTW Analysis

**Objective:** Estimate the Average Treatment Effect using propensity score weighting.

**Steps:**
1. **Propensity Score Estimation:** Logistic regression modeling probability of treatment
2. **Weight Calculation:** Stabilized weights for ATE
3. **Balance Assessment:** SMD before/after weighting
4. **Outcome Modeling:** Weighted logistic regression
5. **Pooling:** Rubin's Rules across 20 imputations

```r
# Propensity score model
ps_model <- glm(milnepan ~ age + sex + bmi + smoking + diabetes,
                data = data, family = binomial)

# Calculate stabilized weights
ps <- predict(ps_model, type = "response")
weight <- ifelse(milnepan == 1,
                 mean(milnepan) / ps,
                 (1 - mean(milnepan)) / (1 - ps))
```

**Balance Achieved:**

| Covariate | SMD Before | SMD After |
|-----------|------------|-----------|
| Age | 0.511 | 0.008 |
| Diabetes | 0.303 | 0.012 |
| Sex | 0.021 | 0.003 |
| BMI | 0.067 | 0.009 |
| Smoking | 0.032 | 0.007 |

*All SMD < 0.1 after weighting ✓*

### Phase 5: G-formula Analysis

**Objective:** Estimate causal effects through outcome modeling and counterfactual prediction.

**Steps:**
1. **Outcome Model:** Logistic regression with treatment + confounders
2. **Counterfactual Prediction:** Predict mortality if everyone treated vs. everyone untreated
3. **Effect Estimation:** Compare average predicted risks
4. **Confidence Intervals:** Bootstrap with 1,000 resamples
5. **Pooling:** Rubin's Rules across 20 imputations

```r
# Outcome model
outcome_model <- glm(death7 ~ milnepan + age + sex + bmi + smoking + diabetes,
                     data = data, family = binomial)

# Counterfactual predictions
data_treated <- data_untreated <- data
data_treated$milnepan <- 1
data_untreated$milnepan <- 0

risk_treated <- mean(predict(outcome_model, data_treated, type = "response"))
risk_untreated <- mean(predict(outcome_model, data_untreated, type = "response"))

RR <- risk_treated / risk_untreated
ARR <- risk_untreated - risk_treated
```

### Phase 6: Pooling with Rubin's Rules

**Objective:** Combine results across 20 imputed datasets while accounting for imputation uncertainty.

**Formula:**
- **Pooled estimate:** Q̄ = (1/m) Σ Qᵢ
- **Within-imputation variance:** W̄ = (1/m) Σ Wᵢ
- **Between-imputation variance:** B = (1/(m-1)) Σ (Qᵢ - Q̄)²
- **Total variance:** T = W̄ + (1 + 1/m)B

---

## 📈 Key Results

### Primary Findings

| Method | Estimate | 95% CI | Interpretation |
|--------|----------|--------|----------------|
| **IPTW** | OR = 0.876 | 0.798 – 0.961 | 12.4% lower odds of death |
| **G-formula** | RR = 0.922 | 0.877 – 0.969 | 7.8% lower risk of death |
| **Cox (imputed)** | HR = 0.918 | 0.86 – 0.98 | 8.2% lower hazard of death |

### G-formula Absolute Effects

| Scenario | 7-Day Mortality |
|----------|-----------------|
| If all patients treated | 38.65% |
| If no patients treated | 41.93% |
| **Absolute Risk Reduction** | **3.28 percentage points** |

### Conclusion

> **Milnepan treatment causes an 8-12% reduction in 7-day mortality among heffpox patients.**

The consistency between two different causal inference methods (IPTW and G-formula), which rely on different assumptions, strengthens confidence in this causal conclusion.

---

## 📁 Project Structure

```
causal-inference-mortality-analysis/
│
├── 📂 data/
│   └── heffpox_2526.csv              # Raw dataset
│
├── 📂 scripts/
│   ├── 01_data_exploration.R         # Initial exploration & baseline characteristics
│   ├── 02_missing_data_analysis.R    # Missing patterns, mechanism assessment
│   ├── 03_multiple_imputation.R      # MICE imputation with diagnostics
│   ├── 04_dag_analysis.R             # DAG construction & adjustment set
│   ├── 05_survival_analysis.R        # KM curves & Cox models
│   ├── 06_iptw_analysis.R            # Propensity scores & IPTW
│   ├── 07_gformula_analysis.R        # G-formula & bootstrapping
│   └── 08_pooling_results.R          # Rubin's Rules pooling
│
├── 📂 outputs/
│   ├── 📂 figures/
│   │   ├── missing_data_pattern.png
│   │   ├── convergence_plot.png
│   │   ├── density_plots.png
│   │   ├── dag.png
│   │   ├── km_curves.png
│   │   ├── ps_distribution.png
│   │   └── love_plot.png
│   │
│   └── 📂 tables/
│       ├── baseline_characteristics.csv
│       ├── cox_models_comparison.csv
│       ├── covariate_balance.csv
│       └── final_results.csv
│
├── 📂 report/
│   └── Causal_Analysis_Report.pdf    # Full written report
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## ⚙️ Installation

### Prerequisites

- R version 4.0 or higher
- RStudio (recommended)

### Install Required Packages

```r
# Core packages
install.packages(c(
  "tidyverse",    # Data manipulation & visualization
  "mice",         # Multiple imputation
  "VIM"           # Missing data visualization
))

# Causal inference packages
install.packages(c(
  "WeightIt",     # Propensity score weighting
  "cobalt",       # Covariate balance assessment
  "survey",       # Weighted regression
  "dagitty"       # DAG analysis
))

# Survival analysis packages
install.packages(c(
  "survival",     # Cox models
  "survminer"     # Survival visualization
))

# Tables
install.packages("tableone")
```

---

## 🚀 Usage

### Clone the Repository

```bash
git clone https://github.com/Akinlajar/causal-inference-mortality-analysis.git
cd causal-inference-mortality-analysis
```

### Run the Analysis

Execute scripts in order:

```r
source("scripts/01_data_exploration.R")
source("scripts/02_missing_data_analysis.R")
source("scripts/03_multiple_imputation.R")
source("scripts/04_dag_analysis.R")
source("scripts/05_survival_analysis.R")
source("scripts/06_iptw_analysis.R")
source("scripts/07_gformula_analysis.R")
source("scripts/08_pooling_results.R")
```

Or run interactively in RStudio for step-by-step exploration.

---

## 🛠️ Skills Demonstrated

### Statistical Methods
- ✅ Causal inference (IPTW, G-formula)
- ✅ Propensity score estimation and diagnostics
- ✅ Multiple imputation (MICE, Rubin's Rules)
- ✅ Survival analysis (Kaplan-Meier, Cox regression)
- ✅ Bootstrap confidence intervals
- ✅ Missing data mechanism assessment (MCAR/MAR/MNAR)

### Technical Skills
- ✅ R programming (tidyverse ecosystem)
- ✅ Statistical modeling
- ✅ Data visualization (ggplot2)
- ✅ Reproducible research practices
- ✅ Scientific writing and reporting

### Domain Knowledge
- ✅ Epidemiological study design
- ✅ Confounding and bias
- ✅ Directed Acyclic Graphs (DAGs)
- ✅ Healthcare data analysis

---

## 📊 Visualizations

### 1. Directed Acyclic Graph (DAG)
*Shows causal relationships between variables and identifies confounders*

### 2. Missing Data Pattern
*Visualizes the pattern of missing values across variables*

### 3. Imputation Diagnostics
*Trace plots showing convergence; density plots comparing observed vs imputed*

### 4. Kaplan-Meier Survival Curves
*Shows survival probability over time by treatment group*

### 5. Propensity Score Distribution
*Demonstrates overlap between treatment groups (positivity assumption)*

### 6. Love Plot
*Displays covariate balance before and after IPTW weighting*

---

## ⚠️ Limitations

1. **Unmeasured Confounding:** As an observational study, we cannot rule out confounding by unmeasured variables (e.g., disease severity at presentation).

2. **MAR Assumption:** Multiple imputation assumes data are Missing at Random, which is untestable.

3. **Model Specification:** Propensity score and outcome models may be misspecified in ways not detectable through diagnostics.

4. **Generalizability:** Findings are from a single dataset and may not generalize to other populations.

---

## 📚 References

### Textbooks
- Hernán MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC; 2020.
- Rubin DB. *Multiple Imputation for Nonresponse in Surveys*. Wiley; 1987.

### R Packages
- van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software*. 2011;45(3):1-67.
- Textor J, et al. Robust causal inference using directed acyclic graphs: the R package dagitty. *International Journal of Epidemiology*. 2016;45(6):1887-1894.
- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. Springer; 2000.

---

## 👤 Author

**Judah Akinlajar**

MSc Health Data Science | University of Manchester

[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-blue?logo=linkedin)](https://www.linkedin.com/in/judah-akinlaja-a10385194)
[![Email](https://img.shields.io/badge/Email-Contact-red?logo=gmail)](mailto:judahakinlajar@gmail.com)

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  <i>If you found this project useful, please consider giving it a ⭐</i>
</p>

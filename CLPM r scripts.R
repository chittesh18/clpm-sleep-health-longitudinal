install.packages(c("haven", "lavaan", "tidyverse"))

library(haven)
library(lavaan)
library(tidyverse)

setwd("E:/Thesis/Understanding_society/UKDA-6614-stata/stata/stata14_se/ukhls")
clean_neg <- function(x) ifelse(x < 0, NA, x)

# Load each wave
w1_raw <- read_dta("a_indresp.dta")
w1 <- data.frame(
  pidp  = w1_raw$pidp,
  ghq1  = clean_neg(as.numeric(w1_raw$a_scghq1_dv)),
  slp1  = clean_neg(as.numeric(w1_raw$a_scslp_qual)),
  sex   = clean_neg(as.numeric(w1_raw$a_sex_dv)),
  age   = clean_neg(as.numeric(w1_raw$a_age_dv))
)

w2_raw <- read_dta("d_indresp.dta")
w2 <- data.frame(
  pidp  = w2_raw$pidp,
  ghq2  = clean_neg(as.numeric(w2_raw$d_scghq1_dv)),
  slp2  = clean_neg(as.numeric(w2_raw$d_slp_qual))
)

w3_raw <- read_dta("g_indresp.dta")
w3 <- data.frame(
  pidp  = w3_raw$pidp,
  ghq3  = clean_neg(as.numeric(w3_raw$g_scghq1_dv)),
  slp3  = clean_neg(as.numeric(w3_raw$g_slp_qual))
)

w4_raw <- read_dta("j_indresp.dta")
w4 <- data.frame(
  pidp  = w4_raw$pidp,
  ghq4  = clean_neg(as.numeric(w4_raw$j_scghq1_dv)),
  slp4  = clean_neg(as.numeric(w4_raw$j_slp_qual))
)

w5_raw <- read_dta("m_indresp.dta")
w5 <- data.frame(
  pidp  = w5_raw$pidp,
  ghq5  = clean_neg(as.numeric(w5_raw$m_scghq1_dv)),
  slp5  = clean_neg(as.numeric(w5_raw$m_slp_qual))
)

# Merge
wide_data <- merge(w1, w2, by = "pidp")
wide_data <- merge(wide_data, w3, by = "pidp")
wide_data <- merge(wide_data, w4, by = "pidp")
wide_data <- merge(wide_data, w5, by = "pidp")


#..................................................................... Standardize…………………………………………………
wide_data$ghq1s <- scale(wide_data$ghq1)[,1]
wide_data$ghq2s <- scale(wide_data$ghq2)[,1]
wide_data$ghq3s <- scale(wide_data$ghq3)[,1]
wide_data$ghq4s <- scale(wide_data$ghq4)[,1]
wide_data$ghq5s <- scale(wide_data$ghq5)[,1]

wide_data$slp1s <- scale(wide_data$slp1)[,1]
wide_data$slp2s <- scale(wide_data$slp2)[,1]
wide_data$slp3s <- scale(wide_data$slp3)[,1]
wide_data$slp4s <- scale(wide_data$slp4)[,1]
wide_data$slp5s <- scale(wide_data$slp5)[,1]

# Extract covariates from wave 1
w1_cov <- data.frame(
  pidp    = w1_raw$pidp,
  edu     = clean_neg(as.numeric(w1_raw$a_hiqual_dv)),
  employ  = clean_neg(as.numeric(w1_raw$a_jbstat)),
  nssec   = clean_neg(as.numeric(w1_raw$a_jbnssec5_dv)),
  sex     = clean_neg(as.numeric(w1_raw$a_sex_dv)),
  age     = clean_neg(as.numeric(w1_raw$a_age_dv))
)

# Merge into wide_data
wide_data <- merge(wide_data, w1_cov, by = "pidp", all.x = TRUE)

# Education — recode to 3 levels
# hiqual_dv: 1=degree+, 2=other higher, 3=A-level, 4=GCSE, 5=other qual, 9=no qual
wide_data$edu3 <- ifelse(wide_data$edu %in% c(1,2), 1,      # Higher education
                         ifelse(wide_data$edu %in% c(3,4), 2,      # Secondary/A-level
                                ifelse(wide_data$edu %in% c(5,9), 3, NA))) # Low/no qualification
cat("Education 3-level:\n"); print(table(wide_data$edu3, useNA="always"))

# Employment — recode to working vs not working
# jbstat: 1=employed, 2=self-employed, 3=unemployed, 4=retired, 5=maternity,
#         6=family care, 7=student, 8=sick/disabled, 9=other
wide_data$working <- ifelse(wide_data$employ %in% c(1,2), 1,  # Working
                            ifelse(is.na(wide_data$employ), NA, 0))   # Not working
cat("Working:\n"); print(table(wide_data$working, useNA="always"))

#...........................................CLPM without covariates………………………………………………..
clpm_model <- '
  # Autoregressive and cross-lagged paths
  ghq2s ~ a*ghq1s + c*slp1s 
  ghq3s ~ a*ghq2s + c*slp2s
  ghq4s ~ a*ghq3s + c*slp3s
  ghq5s ~ a*ghq4s + c*slp4s

  slp2s ~ d*slp1s + b*ghq1s 
  slp3s ~ d*slp2s + b*ghq2s
  slp4s ~ d*slp3s + b*ghq3s
  slp5s ~ d*slp4s + b*ghq4s

  # Residual covariances
  ghq2s ~~ slp2s
  ghq3s ~~ slp3s
  ghq4s ~~ slp4s
  ghq5s ~~ slp5s

  # Baseline covariance
  ghq1s ~~ slp1s
'

fit_clpm <- sem(
  clpm_model,
  data      = wide_data,
  estimator = "ML",
  missing   = "FIML"
)

summary(fit_clpm,
        fit.measures = TRUE,
        standardized = TRUE,
        ci           = TRUE)

fitMeasures(fit_clpm, c("cfi", "rmsea", "srmr", "tli"))

parameterEstimates(fit_clpm, standardized = TRUE) |>
  subset(op == "~") |>
  (\(x) x[, c("lhs","rhs","est","se","pvalue","ci.lower","ci.upper","std.all")])() |>
  (\(x) { x[] <- lapply(x, function(v) if(is.numeric(v)) round(v,3) else v); x })()

#...........................................CLPM with covariates………………………………………………………..
clpm_model_with_covariates <- '
  ghq2s ~ a*ghq1s + c*slp1s + age + working
  ghq3s ~ a*ghq2s + c*slp2s
  ghq4s ~ a*ghq3s + c*slp3s
  ghq5s ~ a*ghq4s + c*slp4s

  slp2s ~ d*slp1s + b*ghq1s + edu3 + age + working
  slp3s ~ d*slp2s + b*ghq2s
  slp4s ~ d*slp3s + b*ghq3s
  slp5s ~ d*slp4s + b*ghq4s

  ghq2s ~~ slp2s
  ghq3s ~~ slp3s
  ghq4s ~~ slp4s
  ghq5s ~~ slp5s

  ghq1s ~~ slp1s
'

fit_final <- sem(
  clpm_model_with_covariates,
  data      = wide_data,
  estimator = "ML",
  missing   = "FIML"
)

summary(fit_final, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
fitMeasures(fit_final, c("cfi", "rmsea", "srmr", "tli"))

# Key paths only
parameterEstimates(fit_final, standardized = TRUE) |>
  subset(op == "~") |>
  (\(x) x[, c("lhs","rhs","est","se","pvalue","ci.lower","ci.upper","std.all")])() |>
  (\(x) { x[] <- lapply(x, function(v) if(is.numeric(v)) round(v,3) else v); x })()


# ── Create factor variables for grouping ──
wide_data$sex_group <- factor(wide_data$sex,
                              levels = c(1, 2),
                              labels = c("Male", "Female"))

wide_data$work_group <- factor(wide_data$working,
                               levels = c(0, 1),
                               labels = c("Not working", "Working"))

wide_data$edu_group <- factor(wide_data$edu3,
                              levels = c(1, 2, 3),
                              labels = c("Higher education",
                                         "Secondary level",
                                         "Low/no qualification"))

#Moderation
# ── Model for employment moderation ──
# working is now the GROUP so remove from covariates
clpm_work_mg <- '
  ghq2s ~ a*ghq1s + c*slp1s + age + sex
  ghq3s ~ a*ghq2s + c*slp2s
  ghq4s ~ a*ghq3s + c*slp3s
  ghq5s ~ a*ghq4s + c*slp4s

  slp2s ~ d*slp1s + b*ghq1s + edu3 + age + sex
  slp3s ~ d*slp2s + b*ghq2s
  slp4s ~ d*slp3s + b*ghq3s
  slp5s ~ d*slp4s + b*ghq4s

  ghq2s ~~ slp2s
  ghq3s ~~ slp3s
  ghq4s ~~ slp4s
  ghq5s ~~ slp5s
  ghq1s ~~ slp1s
'

# Configural
fit_work_config <- sem(clpm_work_mg,
                       data      = wide_data,
                       group     = "work_group",
                       estimator = "ML",
                       missing   = "FIML")

# Constrained
fit_work_con <- sem(clpm_work_mg,
                    data        = wide_data,
                    group       = "work_group",
                    group.equal = c("regressions"),
                    estimator   = "ML",
                    missing     = "FIML")

cat("\n=== EMPLOYMENT MODERATION LR TEST ===\n")
print(lavTestLRT(fit_work_config, fit_work_con))
cat("\nFit indices (configural):\n")
print(fitMeasures(fit_work_config, c("cfi","rmsea","srmr","tli")))

# Extract cross-lagged paths by employment group
cat("\n=== CROSS-LAGGED PATHS BY EMPLOYMENT ===\n")
for (g in 1:2) {
  grp_name <- c("Not working","Working")[g]
  cat("\n---", grp_name, "---\n")
  parameterEstimates(fit_work_config, standardized=TRUE) |>
    subset(op=="~" & group==g) |>
    subset(
      (lhs=="ghq2s" & rhs=="slp1s") |
        (lhs=="ghq3s" & rhs=="slp2s") |
        (lhs=="ghq4s" & rhs=="slp3s") |
        (lhs=="ghq5s" & rhs=="slp4s") |
        (lhs=="slp2s" & rhs=="ghq1s") |
        (lhs=="slp3s" & rhs=="ghq2s") |
        (lhs=="slp4s" & rhs=="ghq3s") |
        (lhs=="slp5s" & rhs=="ghq4s")
    ) |>
    (\(x) x[, c("lhs","rhs","est","se","pvalue","std.all")])() |>
    (\(x) { x[] <- lapply(x, function(v) if(is.numeric(v)) round(v,3) else v); x })() |>
    print()
}

# ── Model for education moderation ──
# education is now the GROUP so remove from covariates
clpm_edu_mg <- '
  ghq2s ~ a*ghq1s + c*slp1s + age + sex + working
  ghq3s ~ a*ghq2s + c*slp2s
  ghq4s ~ a*ghq3s + c*slp3s
  ghq5s ~ a*ghq4s + c*slp4s

  slp2s ~ d*slp1s + b*ghq1s + age + sex + working
  slp3s ~ d*slp2s + b*ghq2s
  slp4s ~ d*slp3s + b*ghq3s
  slp5s ~ d*slp4s + b*ghq4s

  ghq2s ~~ slp2s
  ghq3s ~~ slp3s
  ghq4s ~~ slp4s
  ghq5s ~~ slp5s
  ghq1s ~~ slp1s
'

# Configural
fit_edu_config <- sem(clpm_edu_mg,
                      data      = wide_data,
                      group     = "edu_group",
                      estimator = "ML",
                      missing   = "FIML")

# Constrained
fit_edu_con <- sem(clpm_edu_mg,
                   data        = wide_data,
                   group       = "edu_group",
                   group.equal = c("regressions"),
                   estimator   = "ML",
                   missing     = "FIML")

cat("\n=== EDUCATION MODERATION LR TEST ===\n")
print(lavTestLRT(fit_edu_config, fit_edu_con))
cat("\nFit indices (configural):\n")
print(fitMeasures(fit_edu_config, c("cfi","rmsea","srmr","tli")))

# Extract cross-lagged paths by education group
cat("\n=== CROSS-LAGGED PATHS BY EDUCATION ===\n")
for (g in 1:3) {
  grp_name <- c("Higher education","Secondary level","Low/no qualification")[g]
  cat("\n---", grp_name, "---\n")
  parameterEstimates(fit_edu_config, standardized=TRUE) |>
    subset(op=="~" & group==g) |>
    subset(
      (lhs=="ghq2s" & rhs=="slp1s") |
        (lhs=="ghq3s" & rhs=="slp2s") |
        (lhs=="ghq4s" & rhs=="slp3s") |
        (lhs=="ghq5s" & rhs=="slp4s") |
        (lhs=="slp2s" & rhs=="ghq1s") |
        (lhs=="slp3s" & rhs=="ghq2s") |
        (lhs=="slp4s" & rhs=="ghq3s") |
        (lhs=="slp5s" & rhs=="ghq4s")
    ) |>
    (\(x) x[, c("lhs","rhs","est","se","pvalue","std.all")])() |>
    (\(x) { x[] <- lapply(x, function(v) if(is.numeric(v)) round(v,3) else v); x })() |>
    print()
}


cat("Missing edu_group:", sum(is.na(wide_data$edu_group)), "\n")
cat("N used in edu model:", 
    sum(lavInspect(fit_edu_config, "nobs")), "\n")
cat("N in full sample:", nrow(wide_data), "\n")

lavInspect(fit_clpm, "options")$estimator


# ── Model for sex moderation ──
# sex is now the GROUP so remove from covariates
clpm_sex_mg <- '
  ghq2s ~ a*ghq1s + c*slp1s + age + edu3 + working
  ghq3s ~ a*ghq2s + c*slp2s
  ghq4s ~ a*ghq3s + c*slp3s
  ghq5s ~ a*ghq4s + c*slp4s

  slp2s ~ d*slp1s + b*ghq1s + age + edu3 + working
  slp3s ~ d*slp2s + b*ghq2s
  slp4s ~ d*slp3s + b*ghq3s
  slp5s ~ d*slp4s + b*ghq4s

  ghq2s ~~ slp2s
  ghq3s ~~ slp3s
  ghq4s ~~ slp4s
  ghq5s ~~ slp5s
  ghq1s ~~ slp1s
'

# Configural model
fit_sex_config <- sem(clpm_sex_mg,
                      data      = wide_data,
                      group     = "sex_group",
                      estimator = "ML",
                      missing   = "FIML")

# Constrained model
fit_sex_con <- sem(clpm_sex_mg,
                   data        = wide_data,
                   group       = "sex_group",
                   group.equal = c("regressions"),
                   estimator   = "ML",
                   missing     = "FIML")

cat("\n=== SEX MODERATION LR TEST ===\n")
print(lavTestLRT(fit_sex_config, fit_sex_con))

cat("\nFit indices (configural):\n")
print(fitMeasures(fit_sex_config,
                  c("cfi","rmsea","srmr","tli")))

# Extract cross-lagged paths by sex group
cat("\n=== CROSS-LAGGED PATHS BY SEX ===\n")

for (g in 1:2) {
  
  grp_name <- c("Male","Female")[g]
  
  cat("\n---", grp_name, "---\n")
  
  parameterEstimates(fit_sex_config,
                     standardized = TRUE) |>
    
    subset(op == "~" & group == g) |>
    
    subset(
      (lhs=="ghq2s" & rhs=="slp1s") |
        (lhs=="ghq3s" & rhs=="slp2s") |
        (lhs=="ghq4s" & rhs=="slp3s") |
        (lhs=="ghq5s" & rhs=="slp4s") |
        (lhs=="slp2s" & rhs=="ghq1s") |
        (lhs=="slp3s" & rhs=="ghq2s") |
        (lhs=="slp4s" & rhs=="ghq3s") |
        (lhs=="slp5s" & rhs=="ghq4s")
    ) |>
    
    (\(x) x[, c("lhs","rhs","est","se","pvalue","std.all")])() |>
    
    (\(x) {
      x[] <- lapply(x,
                    function(v)
                      if(is.numeric(v)) round(v,3) else v)
      x
    })() |>
    
    print()
}

cat("Missing sex_group:",
    sum(is.na(wide_data$sex_group)), "\n")

cat("N used in sex model:",
    sum(lavInspect(fit_sex_config, "nobs")), "\n")

cat("N in full sample:",
    nrow(wide_data), "\n")

lavInspect(fit_sex_config, "options")$estimator

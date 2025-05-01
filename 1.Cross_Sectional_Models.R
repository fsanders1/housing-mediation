library(psych)
library(relaimpo)
library(ggplot2)

# Load data

load(".../Housing_Data_with_SVs.RData")

# Descriptives

model1 <- c("depression_age30", "housing_age30", Covars_model1) 

describe(df[model1])

model2 <- c("depression_age32", "housing_age32", Covars_model2) 

describe(df[model2])

model3 <- c("depression_age33", "housing_age33", Covars_model3) 

describe(df[model3])

# Percentages for dep_fam_risk
## Calculate frequencies
freq <- table(df$dep_fam_risk)

## Calculate percentages
percentages <- prop.table(freq) * 100
percentages

# Additional descriptives for descriptives table that aren't included in models
additional <- c("age_at_depression_age30", "age_at_depression_age32", "age_at_depression_age33") 
describe(df[additional])

# Percentages for smoking
## Calculate frequencies
freq2 <- table(df$smoking_age30)

## Calculate percentages
percentages <- prop.table(freq2) * 100
percentages

# Percentages for ethnicity
## Calculate frequencies
freq3 <- table(df$ethnicity_mum)

## Calculate percentages
percentages <- prop.table(freq3) * 100
percentages

# Plots

ggplot(df, aes(x = housing_age30, y = depression_age33)) + 
  geom_point(size = 1) +
  geom_smooth(method = lm) +
  labs(y = "Depressive Symptoms (33 yrs)") +
  labs(x = "Poor Housing Quality (30 yrs)") +
  theme(
    axis.title.x = element_text(size = 20),  # Adjust x-axis title size
    axis.title.y = element_text(size = 20)   # Adjust y-axis title size
  )

ggplot(df, aes(x = housing_fg_score, y = depression_age48)) + 
  geom_point(size = 1) +
  geom_smooth(method = lm) +
  labs(y = "Depressive Symptoms (48 yrs)") +
  labs(x = "Poor Housing Quality (32-33 yrs)") +
  theme(
    axis.title.x = element_text(size = 20),  # Adjust x-axis title size
    axis.title.y = element_text(size = 20)   # Adjust y-axis title size
  )

# Standardised Regressions.

# Standardise the data

vars_to_standardize <- c(model1,model2,model3)
vars_to_standardize <- unique(vars_to_standardize)
st_df <- df
st_df[vars_to_standardize] <- scale(df[vars_to_standardize])

# Time point 1 mothers = 30 years.

formula <- as.formula(paste("depression_age30 ~ housing_age30 +", paste(Covars_model1, collapse = " + ")))
modelage30 <- lm(formula, data = st_df)

summary(modelage30)

# Time point 2 mothers = 32 years.
formula2 <- as.formula(paste("depression_age32 ~ housing_age32 +", paste(Covars_model2, collapse = " + ")))
modelage32 <- lm(formula2, data = st_df)

summary(modelage32)

# Time point 3 mothers = 33 years.
formula3 <- as.formula(paste("depression_age33 ~ housing_age33 +", paste(Covars_model3, collapse = " + ")))
modelage33 <- lm(formula3, data = st_df)

summary(modelage33)

# Calculate semi-partial correlations
rel_imp <- calc.relimp(modelage30, type = "lmg", rela = FALSE)
rel_imp2 <- calc.relimp(modelage32, type = "lmg", rela = FALSE)
rel_imp3 <- calc.relimp(modelage33, type = "lmg", rela = FALSE)

# Extract the semi-partial R-squared values (lmg values are sr^2)
sr2_values <- rel_imp$lmg
print(sr2_values)
sr2_values2 <- rel_imp2$lmg
print(sr2_values2)
sr2_values3 <- rel_imp3$lmg
print(sr2_values3)

# Confidence intervals for p values
confint(modelage30, level = 0.95)
confint(modelage32, level = 0.95)
confint(modelage33, level = 0.95)

# Unadjusted models

model1 <- lm(depression_age30 ~ housing_age30, data = st_df)
summary(model1)
model2 <- lm(depression_age32 ~ housing_age32, data = st_df)
summary(model2)
model3 <- lm(depression_age33 ~ housing_age33, data = st_df)
summary(model3)

# Confidence intervals for p values
confint(model1, level = 0.95)
confint(model2, level = 0.95)
confint(model3, level = 0.95)

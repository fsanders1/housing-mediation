library(dplyr)
library(foreign)
library(psych)
library(relaimpo)
library(corrplot)
library(ggplot2)

# Load data

load(".../Housing_Data_with_SVs.RData")


################################# Descriptives ##########################################

model4 <- c("depression_age33", "housing_age30", Covars_model4_6) 

describe(df[model4])

model5 <- c("depression_age48", "housing_fg_score", Covars_model5_7) 

describe(df[model5])


################################### Model 4 ##############################################

# Standardise data
vars_to_standardize <- c(model4,model5)
vars_to_standardize <- unique(vars_to_standardize)
st_df <- df
st_df[vars_to_standardize] <- scale(df[vars_to_standardize])

# Run model
formula <- as.formula(paste("depression_age33 ~ housing_age30 +", paste(Covars_model4_6, collapse = " + ")))
model4r <- lm(formula, data = st_df)

# Show results
summary(model4r)

#################################### Model 5 ##############################################

# Run model
formula2 <- as.formula(paste("depression_age48 ~ housing_fg_score +", paste(Covars_model5_7, collapse = " + ")))
model5r <- lm(formula2, data = st_df)

# Show results
summary(model5r)

# Calculate semi-partial correlations
rel_imp <- calc.relimp(model4r, type = "lmg", rela = FALSE)
rel_imp2 <- calc.relimp(model5r, type = "lmg", rela = FALSE)

# Extract the semi-partial R-squared values (lmg values are sr^2)
sr2_values <- rel_imp$lmg
print(sr2_values)
sr2_values2 <- rel_imp2$lmg
print(sr2_values2)

# Confidence intervals for p values
confint(model4r, level = 0.95)
confint(model5r, level = 0.95)

#################################### Unadjusted models ####################################

model1 <- lm(depression_age33 ~ housing_age30, data = st_df)
summary(model1)
model2 <- lm(depression_age48 ~ housing_fg_score, data = st_df)
summary(model2)

# Confidence intervals for p values
confint(model1, level = 0.95)
confint(model2, level = 0.95)

################################# Plots ###################################################

# Correlation Matrices
m4 = cor(st_df[model4])

title <- "Correlation Matrix for Model 4"
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(m4, method="color", col=col(200),  
         diag=FALSE, # tl.pos="d", 
         type="upper", order="hclust", 
         addCoef.col = "black",
         tl.pos = "lt", tl.col = "black")

m5 = cor(st_df[model5])

title <- "Correlation Matrix for Model 5"

corrplot(m5, method="color", col=col(200),  
         diag=FALSE, # tl.pos="d", 
         type="upper",  
         addCoef.col = "black",
         tl.pos = "lt", tl.col = "black")

# Plot housing against depressive symptoms for models 4 and 5

scatterplot1 <- ggplot(df, aes(x = housing_age30, y = depression_age33)) +
  geom_point(size = 0.25) +  # Add points
  geom_smooth(method = "lm", se = TRUE) + 
  theme_classic() +  # Use classic theme without gridlines
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("Poor Housing Quality (30 yrs)") +  # Change x-axis title
  ylab("Depressive Symptoms (33 yrs)")


scatterplot2 <- ggplot(df, aes(x = housing_fg_score, y = depression_age48)) +
  geom_point(size = 0.25) +  # Add points
  geom_smooth(method = "lm", se = TRUE) + 
  theme_classic() +  # Use classic theme without gridlines
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("Poor Housing Quality (32-33 yrs)") +  # Change x-axis title
  ylab("Depressive Symptoms (48yrs)")

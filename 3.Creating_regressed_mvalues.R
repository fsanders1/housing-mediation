# In this script, we residualise out our covariates from our data

library(dplyr)
library(psych)

# Load in data
load(".../data.RData")

################### Descriptives ###################################

for (var in c("depression_age33", "housing_age30", Covars_model4_6)) {
  cat("Descriptive statistics for", var, ":\n")
  print(describe(med.sample.data[[var]]))
  cat("\n")
}

####################################################################

## regress out technical variables and confounders from df subsetted to each adversity's complete cases -----------------
cat("Regress out tech vars and confounders from df SUBSET \n")

# Specifying our covariates
covars <-   c(Covars_model4_6, Covars_DNAm.30y)

# sort them both here to make sure in the same order for later regressions
med.sample.data <- med.sample.data %>% arrange(ID)
med.DNAm.data<- med.DNAm.data %>% arrange(ID)
if (identical(med.sample.data$ID, med.DNAm.data$ID)==FALSE){
  stop("df and beta matrices do not align by ID number.")
}

# Regression specifically for adversity's complete cases
cat("Regress out tech and confounders from df \n")

print(dim(med.DNAm.data)) 
if(identical(med.sample.data$ID, med.DNAm.data$ID)==FALSE){ 
  stop("regression var and betas.pheno.F7.adver matrices do not align by ID number.")
}
med.DNAm.data = subset(med.DNAm.data, select = -ID)# remove ID column for next step
# Remove columns from age 48 timepoint that we don't need for age 30 time point 

columns_to_remove <- c("V1_age48", "V2_age48","V3_age48", "V4_age48","V5_age48", "V6_age48",
                       "V7_age48", "V8_age48","V9_age48", "V10_age48","cellcounts.Bcell_age48", 
                       "cellcounts.CD4T_age48","cellcounts.CD8T_age48","cellcounts.Gran_age48",
                       "cellcounts.Mono_age48", "cellcounts.NK_age48")
med.sample.data <- med.sample.data[, !names(med.sample.data) %in% columns_to_remove]

# Remove cpg columns with all na's. This step is needed, despite imputing missing 
# cpg's in previous script, as it cannot impute cpg's with no data so we have to remove these.
med.DNAm.data <- med.DNAm.data[ , colSums(is.na(med.DNAm.data)) < nrow(med.DNAm.data)]


#regressing covariates
beta.corr.new <- do.call(cbind, lapply(1:ncol(med.DNAm.data), function(x) {
  if (x %% 1000 == 0) {
    print(x)
  }
  dnam <- med.DNAm.data[, x]
  my_formula <- as.formula(paste("dnam ~ ", 
                                 paste(Covars_model4_6, collapse = " + "), 
                                 "+",
                                 paste(Covars_DNAm.30y, collapse = " + ")))
  res <- residuals(lm(my_formula, data = med.sample.data,
                      na.action = na.exclude)) 
  adj.betas <- res + mean(dnam, na.rm = TRUE)
  adj.betas
}))

dim(beta.corr.new)
colnames(beta.corr.new) <- colnames(med.DNAm.data)
rownames(beta.corr.new) <- rownames(med.DNAm.data)
# beta values must stay within the 0-1 range, so changing any of those outside of that to the max/min value in the df
options(digits = 22)
beta.corr.new[which(beta.corr.new>=1)] <- max(beta.corr.new[which(beta.corr.new<1)])
beta.corr.new[which(beta.corr.new<=0)] <- min(beta.corr.new[which(beta.corr.new>0)])

# conversion of beta residuals to M values because M values have better statistical properties
cat("Conversion of beta residuals to M values \n")
adj.m <- lumi::beta2m(beta.corr.new) #convert to m-values

# check no infinite values
if (sum(is.infinite(adj.m))>0){
  stop("Error in beta to M-value conversion, infinite values found")
}

save(adj.m, file = 'Adjusted_M-value_dataframes.Rdata')


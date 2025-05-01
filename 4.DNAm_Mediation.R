suppressMessages(library(lumi))
suppressMessages(library(knitr))
suppressMessages(library(dplyr))
suppressMessages(library(regmed))
suppressMessages(library(lavaan))

############################## Part One #######################################################################

# Remove sex probes
## Load in sex probes
sex_probes <- read.csv("Annotation_XY_450.csv")
sex_probes <- sex_probes$probeID
columns_to_keep <- !names(med.DNAm.data) %in% sex_probes
med.DNAm.data <- med.DNAm.data[, columns_to_keep]

model_data <- model_data[model_data$ID %in% med.DNAm.data$ID, ]
df <- model_data

# Specify our covars
covars <-   c(Covars_model4_6, Covars_DNAm.30y)

# regmed will get you the selected mediation sites
# lavaan will get you the SEs for those estimates

# preparing outcome and exposure data for analysis ---------------------------------------------------

options(digits = 7)

# load appropriate adj.m previously created in "create_regressed_mvalues" function
load(".../Adjusted_M-value_dataframes.Rdata")

## Remove sex probes from our probes
adj.m <- as.data.frame(adj.m)
columns_to_keep <- !names(adj.m) %in% sex_probes
adj.m <- adj.m[, columns_to_keep]
rm(sex_probes)
gc()

# conditioning out confounders from analysis -----------------------------------------------------
# all model inputs must be regressed on confounders prior to analysis
# all model inputs must be standardized prior to analysis to help with later interpretation of results
#cat("conditioning out confounders from outcome and exposures \n")

## regess y.response on covariates to get residuals as y-adjusted outcome
formula <- as.formula(paste("depression_age33 ~", paste(Covars_model4_6, collapse = " + ")))
fit <- lm(formula, data = df)
y.resid <- fit$residuals

## center and scale y - standardization
y.std <- (y.resid - mean(y.resid))/sd(y.resid)

## PURIST APPROACH
# regress exposures on covariates
x <- as.numeric(as.character(df$housing_age30))
formula2 <- as.formula(paste("housing_age30 ~", paste(Covars_model4_6, collapse = " + ")))
fit.x <- lm(formula2, data = df)
x.resid <- fit.x$residuals

# center and scale x
x.std <- (x.resid - mean(x.resid))/sd(x.resid)

# standardize med before calculating correlation - this helps values converge - same matrix as before except reduced by number of sex/nonvar probes
med <- scale(adj.m)

# check if successful med standardization
# should get mean of 0 and sd of 1 - this is how you know data is standardized. These values were once normalized,
# but that's not the same thing. normalized means it was scaled to a value between 0 and 1.
colMeans(med)  # faster version of apply(scaled.dat, 2, mean)
apply(med, 2, sd) 

# Perform SIS -------------------------------------------------------------------
#cat("Perform SIS \n")
print(dim(x.std))
print(dim(med))
# finding correlations between x and meds and meds and y to identify CpGs mediation is most likely to come from
r.x <- cor(x.std, med)
r.y <- cor(y.std, med)
r.xy <- abs(r.x*r.y)

## want to choose those with higher value for r.xy, which is what rank does
rr.xy <- rank(r.xy) # with rank, 1 = lowest

## choose how many top potential mediators to analyze
nmed.test <- floor(434) # 434 = (n-4)/2 (where n = your sample size)

## logical, chooses “top” mediators by highest r.xy.
is.top <- rr.xy > (length(rr.xy) - nmed.test)

# mediators to keep
med.df <- as.data.frame(med[, is.top])

## Implementing regmed package ------------------------------------------------------------

# saving inputs for later manipulation
inputs <- cbind(x.std, y.std, med.df)
save(inputs, file = "med_step2_inputs.Rdata")

##### Due to changes in df structure have to look at ENTIRE grid again ###  ------------
# Fit grid of values to find ideal lambda.vec for this data. frac.lasso =.8 based on simulation studies id-ing this as ideal value

fit.grid <- regmed.grid(x.std, med.df, y.std, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)

# using this grid can see 0.1 is the ideal lambda value (minimum bic). 

# fit grid with lambda = 0.1 -------------------
# fitting grid of values to find ideal lambda.vec for this data. lasso =.8 based on simulation studies id-ing this as ideal value

fit.grid2 <- regmed.grid(x.std, med.df, y.std, lambda.vec= 0.1, frac.lasso=.8)
save(fit.grid2, file = "fit.grid2s_depression.Rdata")

############################## Part Two #######################################################################

meds <- as.matrix(inputs[,-c(1:2)]) 
dim(meds)
x <- inputs[,1]
y <- inputs[,2]

fit.trim2 <- regmed.grid.bestfit(fit.grid2) 
summary.regmed(fit.trim2) # how many mediators
summary(regmed.grid.bestfit(fit.grid2)) 

mediator.epsilon <- 1e-07
keep.med<-data.frame((abs(fit.trim2$alpha*fit.trim2$beta) >= mediator.epsilon))

column_to_check <- "x.std"
# Remove rows where the specified column has FALSE values
filtered_df <- keep.med[keep.med[[column_to_check]], , drop = FALSE]

which.med <- colnames(meds) %in% rownames(filtered_df)

# code to use for fit.trim2:
med.selected <- meds[, which.med]

# relaxed lasso fit - model is refit without the penalization term to get unpenalized effect estimates 
fit.lam0 <- regmed.fit(x, med.selected, y, lambda = 0, frac.lasso = 0.8)
sum.fit.lam0 <- summary.regmed(fit.lam0)
write.csv(sum.fit.lam0, file = 'summ.fit.lam0.csv')

#the above will be used to generate the p-value in the next script

## Code to use to estimate parameter SEs ---------------------------------------------------------

## choose subset of mediators that are in fit.trim
mediator.names <- dimnames(med.selected)[[2]]

## setup lavaan model
#THIS IS SOURCED FROM THE OLD REGMED PACKAGE
source(".../lavaan.model.R")


med.model <- lavaan.model(y.name = "y", x.name = "x", med.name = mediator.names, 
                          medcov = fit.trim2$MedCov)

## note: fit.trim2 is the same as fit.lam0

## setup data for lavaan
dat <- data.frame(cbind(scale(x, center = TRUE, scale = TRUE), scale(y, center = TRUE, 
                                                                     scale = TRUE), scale(med.selected, center = TRUE, scale = TRUE)))
names(dat) <- c("x", "y", dimnames(med.selected)[[2]])

## fit sem with lavaan - read lavaan documentation to interpret the output
gc() # free unused memory

fit.lavaan <- sem(model = med.model, data = dat)
sum.fit.lavaan <- summary(fit.lavaan)
write.csv(sum.fit.lavaan, file = 'sum.fit.lavaan.csv') 


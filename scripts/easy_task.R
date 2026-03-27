# ==============================
# EASY TASK: gamboostLSS Example
# ==============================

#Short-Explanation:
#================================================================
# Objective:
# This task demonstrates how to apply the gamboostLSS model
# using a Gaussian distribution on the mtcars dataset.

# Description:
# The goal is to predict the response variable 'mpg'
# (miles per gallon) using predictor variables such as
# horsepower (hp) and number of cylinders (cyl).

# Approach:
# - Load required libraries
# - Use built-in dataset (mtcars)
# - Fit GaussianLSS model using gamboostLSS
# - Apply cross-validation to find optimal mstop
# - Improve model performance and avoid overfitting

# Outcome:
# The model successfully fits the data and selects optimal
# boosting iterations using cross-validation.
# ================================================================


# Install required packages (run once)
# install.packages("gamboostLSS")
# install.packages("mlbench")

# Load libraries
library(gamboostLSS)
library(mboost)

# Load dataset (mtcars is built-in)
data("mtcars")

# Define response variable
# mpg = miles per gallon
# Using all other variables as predictors
df <- mtcars

# Convert to proper format
df$mpg <- as.numeric(df$mpg)

# ------------------------------
# Fit GaussianLSS Model
# ------------------------------

model <- gamboostLSS(
  mpg ~ wt + hp,   # fewer variables
  data = df,
  families = GaussianLSS(),
  control = boost_control(mstop = 100, nu = 0.1)
)

# ------------------------------
# Cross-validation to find mstop
# ------------------------------

# 10-fold cross-validation
cv <- cvrisk(model, folds = cv(model.weights(model), type = "kfold"))

# Plot CV results
plot(cv)

# Save plot as image
png("plots/easy_plot.png")
plot(cv)   
dev.off()

# Get optimal mstop
mstop_opt <- mstop(cv)
mstop_opt

# Apply optimal mstop
model[mstop_opt]

# ------------------------------
# Selected Variables
# ------------------------------

# Coefficients for mean (mu)
coef(model, parameter = "mu")

# Coefficients for variance (sigma)
coef(model, parameter = "sigma")

# ------------------------------
# Summary
# ------------------------------
summary(model)

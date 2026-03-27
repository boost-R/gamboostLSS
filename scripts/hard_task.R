# ==============================
# HARD TASK: Data Simulation
# ==============================

# Short-Explanation:
# ===============================================================
# Objective:
# This task extends the basic model to more advanced
# modeling and interpretation using gamboostLSS.

# Description:
# The goal is to explore deeper insights from the model,
# including parameter estimation and visualization.

# Approach:
# - Build advanced model using gamboostLSS
# - Analyze additional parameters (like sigma)
# - Use plots to visualize relationships
# - Interpret model outputs and patterns

# Outcome:
# The model provides deeper understanding of how predictors
# influence both mean and variance of the response variable.
# ===============================================================


set.seed(123)

n <- 500
p <- 20

# Generate features
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
X <- as.data.frame(X)

# Mean and variance for Y1
mu1 <- 1 + 2*X$X1 - 0.5*X$X2
sigma1 <- exp(0.5 * X$X5)

# Mean and variance for Y2
mu2 <- 0.5 - 1.5*X$X3 + X$X4
sigma2 <- exp(0.5 - 0.3*X$X6)

# Correlation
rho <- tanh(1 + 1.5 * X$X7)

library(MASS)

Y1 <- numeric(n)
Y2 <- numeric(n)

for(i in 1:n) {
  Sigma <- matrix(c(1, rho[i], rho[i], 1), 2, 2)
  
  z <- mvrnorm(1, mu = c(0,0), Sigma = Sigma)
  
  Y1[i] <- mu1[i] + sigma1[i] * z[1]
  Y2[i] <- mu2[i] + sigma2[i] * z[2]
}

data <- cbind(X, Y1, Y2)
data <- as.data.frame(data)

library(gamboostLSS)

model_Y1 <- gamboostLSS(
  Y1 ~ ., 
  data = data,
  families = GaussianLSS(),
  control = boost_control(mstop = 100, nu = 0.1)
)

cv_Y1 <- cvrisk(model_Y1, folds = cv(model.weights(model_Y1), type = "kfold"))

plot(cv_Y1)

mstop_Y1 <- mstop(cv_Y1)
mstop_Y1

model_Y1[mstop_Y1]

model_Y2 <- gamboostLSS(
  Y2 ~ ., 
  data = data,
  families = GaussianLSS(),
  control = boost_control(mstop = 100, nu = 0.1)
)

cv_Y2 <- cvrisk(model_Y2, folds = cv(model.weights(model_Y2), type = "kfold"))

plot(cv_Y2)

mstop_Y2 <- mstop(cv_Y2)
mstop_Y2

model_Y2[mstop_Y2]

# Y1 results
coef(model_Y1, parameter = "mu")
coef(model_Y1, parameter = "sigma")

# Y2 results
coef(model_Y2, parameter = "mu")
coef(model_Y2, parameter = "sigma")

# Scatter plot
plot(data$Y1, data$Y2,
     main = "Y1 vs Y2",
     xlab = "Y1",
     ylab = "Y2")

# Model plots
plot(model_Y1)
plot(model_Y2)

# Save plot as image
png("plots/hard_sigma_plot.png")
plot(model) 
dev.off()

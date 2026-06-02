require("gamboostLSS")

### check continuous_continuous_loglik() accuracy at theta = 0
theta <- 0
y <- cbind(c(-2, -1, 0, 1, 2), c(-2, -1, 0, 1, 2))
stopifnot(max(abs(gamboostLSS:::continuous_continuous_loglik(y, 
                                                             GaussianLSS(), 
                                                             GaussianLSS(), 
                                                             list(mu=0, sigma=1), 
                                                             list(mu=0, sigma=1), 
                                                             "gaussian", 0) -
                    (gamboostLSS:::get_marginal_logpdf(GaussianLSS(),
                                                       y[,1],
                                                       list(mu=0, sigma=1)) +
                       gamboostLSS:::get_marginal_logpdf(GaussianLSS(),
                                                         y[,2],
                                                         list(mu=0, sigma=1)
                                                         )))) < 1e-6)


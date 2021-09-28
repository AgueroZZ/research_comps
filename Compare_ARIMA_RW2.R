### Study volatility of second order difference
true_W <- compute_g(x)
true_gamma <- diff(true_W, differences = 2)
compute_MSE <- function(vec){
  mean((vec - true_gamma)^2)
}

## RW2
# Mean:
gamma_RW2 <- apply(samps1$samps,2, diff, differences = 2)
gamma_RW2_mean <- apply(gamma_RW2,1, mean)
gamma_RW2_upper <- apply(gamma_RW2,1,quantile, probs = 0.975)
gamma_RW2_lower <- apply(gamma_RW2,1,quantile, probs = 0.025)

plot(gamma_RW2_mean~x[-c(1,2)], type = 'l', col = 'red')
points(true_gamma ~ x[-c(1,2)], col = 'purple')
acf(gamma_RW2_mean, lag.max = 50)


# Samples
plot(gamma_RW2_mean~x[-c(1,2)], type = 'l', col = 'red', ylim = c(-0.1,0.1))
lines(true_gamma ~ x[-c(1,2)], col = 'purple')
lines(gamma_RW2_upper ~ x[-c(1,2)], col = 'orange')
lines(gamma_RW2_lower ~ x[-c(1,2)], col = 'orange')

for (i in sample.int(1000,1)) {
  lines(gamma_RW2[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

acf(gamma_RW2[,1], lag.max = 50)
acf(gamma_RW2[,2], lag.max = 50)
acf(gamma_RW2[,3], lag.max = 50)


MSE_RW2_samples <- apply(gamma_RW2, 2, compute_MSE)
hist(MSE_RW2_samples)

####### Observation: The mean is much much smoother than the sample. The credible interval is quite wide. Sample ACF
####### for each sampled function only spikes at the first lag. 


## ARIMA:

# Mean:
gamma_ARIMA <- apply(samps3$samps,2, diff, differences = 2)
gamma_ARIMA_mean <- apply(gamma_ARIMA,1, mean)
gamma_ARIMA_upper <- apply(gamma_ARIMA,1,quantile, probs = 0.975)
gamma_ARIMA_lower <- apply(gamma_ARIMA,1,quantile, probs = 0.025)

plot(gamma_ARIMA_mean~x[-c(1,2)], type = 'l', col = 'red')
points(true_gamma ~ x[-c(1,2)], col = 'purple')
acf(gamma_RW2_mean, lag.max = 50)

# Samples
plot(gamma_ARIMA_mean~x[-c(1,2)], type = 'l', col = 'red', ylim = c(-0.1,0.1))
lines(true_gamma ~ x[-c(1,2)], col = 'purple')
lines(gamma_ARIMA_upper ~ x[-c(1,2)], col = 'orange')
lines(gamma_ARIMA_lower ~ x[-c(1,2)], col = 'orange')

for (i in sample.int(1000,1)) {
  lines(gamma_ARIMA[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

acf(gamma_ARIMA[,1], lag.max = 50)
acf(gamma_ARIMA[,2], lag.max = 50)
acf(gamma_ARIMA[,3], lag.max = 50)


MSE_ARIMA_samples <- apply(gamma_ARIMA, 2, compute_MSE)

hist_data <- data.frame(MSE = c(MSE_RW2_samples,MSE_ARIMA_samples), types = rep(c("RW2", "ARIMA"), each = 1000))

ggplot(data = hist_data) + geom_density(aes(x=MSE, fill = types), alpha = 0.3) + ylab("Distribution of MSE") + xlab("MSE of each sample")

credible_width_rw2 <- mean(gamma_RW2_upper - gamma_RW2_lower)
credible_width_rw2
credible_width_ARIMA <- mean(gamma_ARIMA_upper - gamma_ARIMA_lower)
credible_width_ARIMA






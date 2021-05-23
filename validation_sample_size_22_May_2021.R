### Sample size for validating a logistic regression based prediction model
### Riley RD, Debray TD, Collins GS, Archer L, Ensor J, van Smeden M, Snell KIE.
### Minimum sample size for external validation of a clinical prediction model
### with a binary outcome. Stat Med 2021 (to appear).
#
# written by Gary Collins (May 2021)
#
# needs PearsonDS package installed (for example in section 4)

### Sample size for O/E
N_OE <- function(phi, SE_OE = 0.051){
  # phi = prevalence
  # SE_OE = standard error of O/E
  round((1 - phi) / (phi * SE_OE^2), 0)
}

### Sample size for the calibration slope
N_slope <- function(LP, beta0 = 0, beta1 = 1, SE_slope = 0.051){
  # LP = linear predictor (vector)
  # beta = calibration intercept 
  # beta = calibration slope
  # SE_slope = standard error of the calibration slope
  I00 <- mean(exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I01 <- mean(LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I11 <- mean(LP * LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  round(I00 / (SE_slope^2 * (I00 * I11 - I01^2)), 0)
}

### Sample size for the c-statistic
N_C <- function(c.index = 0.7, phi = 0.1, target.se = 0.0255, min.opt = 100, max.opt = 10000){
  # needs c.index
  # phi = prevalence
  # target.se = se of the c.index
  se.c <- function(n, cc = c.index, phi2 = phi, target.se2 = target.se){
    zz <- sqrt((cc * (1 - cc)*(1 + (n/2 - 1)*((1 - cc)/(2 - cc)) + ((n/2 - 1) * cc / (1 + cc))))/(n^2 * phi2 * (1 - phi2)))
    abs(zz - target.se2)
  }
  round(optimize(se.c, c(min.opt, max.opt, tol = 0.0001))$minimum, 0)
}

### Sample size for net benefit
N_nb <- function(sens, spec, phi, pt, se.nb = 0.051){
  # sens = sensitivity
  # spec = specificity
  # phi = prevalence
  # pt = threshold
  # se.nb = standard error of net benefit
  w    <- (1 - phi) / phi * (pt / (1 - pt))
  round((sens * (1-sens)/phi + w^2*spec*(1-spec)/(1-phi) + w^2*(1-spec)^2/(phi*(1-phi))) / se.nb^2, 0)
}

# section 3.1.1
N_OE(phi = 0.5, SE_OE = 0.051) # prevalence = 0.5, CI width = 0.2 (i.e., se = 0.051)
N_OE(phi = 0.1, SE_OE = 0.051) # prevalence = 0.1, CI width = 0.2 (i.e., se = 0.051)
N_OE(phi = 0.1, SE_OE = 0.175) # prevalence = 0.1, CI width = 0.7 (i.e., se = 0.175)

# section 3.2.4
N <- 1000000
LP  <- data.frame(id = 1:N, value = rnorm(N, -1.75, 1.47))
p   <- exp(LP$value)/(1 + exp(LP$value))
outcome <- rbinom(N, 1, prob=p)
N_slope(LP$value, SE_slope = 0.051)

# figure 2
cc <- 0.6
phi2 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
n <- c(100, 200, 300, 400, 500, 1000, 1200, 1500)
se.c <- matrix(ncol = length(n), nrow = length(phi2))
for(i in 1:length(phi2)){
  for(j in 1:length(n)){
    se.c[i,j] <- sqrt((cc * (1 - cc)*(1 + (n[j]/2 - 1)*((1 - cc)/(2 - cc)) + ((n[j]/2 - 1) * cc / (1 + cc))))/(n[j]^2 * phi2[i] * (1 - phi2[i])))
  }
}
par(mfrow = c(2,2))
plot(n, se.c[1,], type="l", lty=1, ylim=c(0, 0.1), xlab = 'sample size', ylab = 'SE(C)', main = 'c-statistic = 0.6')
for(i in 2:nrow(se.c)) lines(n, se.c[i,], type='l', lty=i)
abline(h=0.0255)

cc <- 0.7
phi2 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
n <- c(100, 200, 300, 400, 500, 1000, 1200, 1500)
se.c <- matrix(ncol = length(n), nrow = length(phi2))
for(i in 1:length(phi2)){
  for(j in 1:length(n)){
    se.c[i,j] <- sqrt((cc * (1 - cc)*(1 + (n[j]/2 - 1)*((1 - cc)/(2 - cc)) + ((n[j]/2 - 1) * cc / (1 + cc))))/(n[j]^2 * phi2[i] * (1 - phi2[i])))
  }
}
plot(n, se.c[1,], type="l", lty=1, ylim=c(0, 0.1), xlab = 'sample size', ylab = 'SE(C)', main = 'c-statistic = 0.7')
for(i in 2:nrow(se.c)) lines(n, se.c[i,], type='l', lty=i)
abline(h=0.0255)

cc <- 0.8
phi2 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
n <- c(100, 200, 300, 400, 500, 1000, 1200, 1500)
se.c <- matrix(ncol = length(n), nrow = length(phi2))
for(i in 1:length(phi2)){
  for(j in 1:length(n)){
    se.c[i,j] <- sqrt((cc * (1 - cc)*(1 + (n[j]/2 - 1)*((1 - cc)/(2 - cc)) + ((n[j]/2 - 1) * cc / (1 + cc))))/(n[j]^2 * phi2[i] * (1 - phi2[i])))
  }
}
plot(n, se.c[1,], type="l", lty=1, ylim=c(0, 0.1), xlab = 'sample size', ylab = 'SE(C)', main = 'c-statistic = 0.8')
for(i in 2:nrow(se.c)) lines(n, se.c[i,], type='l', lty=i)
abline(h=0.0255)

cc <- 0.9
phi2 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
n <- c(100, 200, 300, 400, 500, 1000, 1200, 1500)
se.c <- matrix(ncol = length(n), nrow = length(phi2))
for(i in 1:length(phi2)){
  for(j in 1:length(n)){
    se.c[i,j] <- sqrt((cc * (1 - cc)*(1 + (n[j]/2 - 1)*((1 - cc)/(2 - cc)) + ((n[j]/2 - 1) * cc / (1 + cc))))/(n[j]^2 * phi2[i] * (1 - phi2[i])))
  }
}
plot(n, se.c[1,], type="l", lty=1, ylim=c(0, 0.1), xlab = 'sample size', ylab = 'SE(C)', main = 'c-statistic = 0.9')
for(i in 2:nrow(se.c)) lines(n, se.c[i,], type='l', lty=i)
abline(h=0.0255)

# section 3.3.1
N_C(c.index = 0.7, phi = 0.1, target.se = 0.0255) 
N_C(c.index = 0.8, phi = 0.5, target.se = 0.0255)

# section 3.4.1
N_nb(sens = 0.6, spec = 0.88, pt = 0.8, phi = 0.72) 

### section 4
moments <- c(mean = -5.8, variance = 5, skewness = -0.5, kurtosis = 4)
LP      <- data.frame(id = 1:N, value = PearsonDS::rpearson(N, moments = moments))
p       <- exp(LP$value)/(1 + exp(LP$value))
outcome <- rbinom(N, 1, prob=p)

N_OE(phi = mean(outcome), SE_OE = 0.245)
N_slope(LP = LP$value)
N_C(c.index = 0.80, phi = 0.018)

# section 4.4 (net benefit)
N_nb(sens = 0.2, spec = 0.99, pt = 0.08, phi = 0.018) 

# section 4.6.3 (only c reported)
c.index <- 0.8
var.eq <- 2 * qnorm(c.index)^2

f.mu2 <- function(x, c.index = 0.8, phi = 0.018, N.POP = 1000000){
  # function to get mu2 in eqn 10
  # phi = prevalence
  # N.POP = sample size of a 'large populaiton' to get mu2
  var.eq <- 2 * qnorm(c.index)^2
  y  <- rbinom(N, 1, phi)
  LP <- rnorm(N.POP, x, sqrt(var.eq))
  LP[y==1] <- rnorm(sum(y), x + sqrt(var.eq), sqrt(var.eq))
  p <- exp(LP) / (1 + exp(LP))
  outcome <- rbinom(N, 1, p)
  OUT <- as.numeric(prop.table(table(outcome))[2])
  abs(OUT - phi)
}

mu2 <- optimize(f.mu2, c(-10, 10, tol = 0.0001))$minimum

N <- 1000000
phi <- 0.018
y <- rbinom(N, 1, phi)
LP <- rnorm(N, mu2, sqrt(var.eq))
LP[y==1] <- rnorm(sum(y), mu2 + sqrt(var.eq), sqrt(var.eq))
p <- exp(LP) / (1 + exp(LP))
outcome2 <- rbinom(N, 1, p)
N_slope(LP = LP)




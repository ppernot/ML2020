# Estimate and test critical value for Shapiro-Wilk test
# Source: Chris C (https://stats.stackexchange.com/users/44952/chris-c),
# How to calculate critical values ($W_{\alpha}$) for the Shapiro-Wilk test?,
# URL (version: 2015-09-12): https://stats.stackexchange.com/q/172155

find.W <- function(alpha = 0.05, error = 0.000001, n = 100){
  not.done <- TRUE
  while(not.done){
    a <- shapiro.test(rnorm(n))
    if(a$p.value < alpha+error && a$p.value > alpha-error){
      not.done = FALSE
      W <- a$statistic
    }
  }
  return(W)
}

alpha = 0.05
n = 5000

# Estimate mean and sd of critical value
wc <- vector()
for(i in 1:100){
  wc[i] <- find.W(alpha = alpha, error = 0.0001, n = n)
}
mwc = mean(wc)
ErrViewLib::prettyUnc(mwc, sd(wc),numDig = 1)

# MC test consistency
nMC = 10000
p = W = c()
for(i in 1:nMC) {
  stat = shapiro.test(rnorm(n))
  p[i] = stat$p.value
  W[i] = stat$statistic
}
# Both rejection probas should be equal
mean(p < alpha)
mean(W < mwc)

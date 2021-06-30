
set.seed(514)
n <- 63
mciter <- 1e4
z <- 1
y <- 1
Time <- 0.25
delta <- Time/n

r <- .05

S0 <- 100
K = 105

sig_high <- 0.3
sig_low <- 0.3


uh <- exp(sig_high*sqrt(delta))
ul <- exp(sig_low*sqrt(delta))

dh <- exp(-sig_high*sqrt(delta))
dl <- exp(-sig_low*sqrt(delta))

qh <- (exp(r*delta) - dh)/(uh - dh)
ql <- (exp(r*delta) - dl)/(ul - dl)
###################################################################

tpm <- matrix(c(.5, .5, .2, .8), nrow = 2, byrow = TRUE)
G <- transition2Generator(P = tpm)
#P <- diag(c(1,1)) + delta*G
P <- expm(delta*G)

ph <- P[1,1]
pl <- P[2,2]
#df <- 1
df = exp(-r*delta)
###########################################################################

###################################################################

## For initial state Low
chain <- c()
chaing <- c()

for(i in 1:mciter) {
  x <- rep(0, n)
  x[1] <- 2
  for(j in 2:n) {
    x[j] <- sample(1:2,1, prob=P[x[j-1],])
  }
  
  s <- rep(0,n)
  sg <- rep(0,n)
  s[1] <- S0
  sg[1] <- S0
  
  for(j in 2:n) {
    if (x[j-1] ==1) s[j] <- s[j-1]*sample(c(uh,dh), 1, prob = c(qh, (1-qh)))
    if (x[j-1] ==2) s[j] <- s[j-1]*sample(c(ul,dl), 1, prob = c(ql, (1-ql)))
    if (x[j-1] ==1) sg[j] <- exp(log(sg[j-1]) + rnorm(1, 0, sig_high*sqrt(delta)))
    if (x[j-1] ==2) sg[j] <- exp(log(sg[j-1]) + rnorm(1, 0, sig_low*sqrt(delta)))
  }
  chain <- c(chain, s[n])
  chaing <- c(chaing, sg[n])
}

pay <- 0.5*(abs(chain - K) + (chain - K)) 
low_MC <- exp(-r)*mean(pay)
payg <- 0.5*(abs(chaing - K) + (chaing - K)) 
low_MCg <- exp(-r)*mean(payg)
#option_low_MC
#option_low

# #initial state 1 (high)
# chain <- c()
# 
# for(i in 1:mciter) {
#   x <- rep(0, n)
#   x[1] <- 1
#   for(j in 2:n) {
#     x[j] <- sample(1:2,1, prob=P[x[j-1],])
#   }
#   
#   s <- rep(0,n)
#   s[1] <- S0
#   
#   for(j in 2:n) {
#     if (x[j-1] ==1) s[j] <- s[j-1]*sample(c(uh,dh), 1, prob = c(qh, (1-qh)))
#     if (x[j-1] ==2) s[j] <- s[j-1]*sample(c(ul,dl), 1, prob = c(ql, (1-ql)))
#   }
#   chain <- c(chain, s[n])
# }
# 
# #chain
# 
# pay <- 0.5*(abs(chain - K) + (chain - K)) 
# high_MC <- exp(-r)*mean(pay)
# #option_high_MC
# #option_high
# 
# 
# #####################################################################
# # plotting the stock prices
# n <- 100
# mciter <- 1000
# z <- 1
# y <- 1
# Time <- 1.0
# delta <- Time/n
# 
# r <- .05
# 
# S0 <- 100
# K = 105
# 
# sig_high <-0.4
# sig_low <- 0.1
# #sig <- 0.3
# 
# uh <- exp(sig_high*sqrt(delta))
# ul <- exp(sig_low*sqrt(delta))
# 
# dh <- exp(-sig_high*sqrt(delta))
# dl <- exp(-sig_low*sqrt(delta))
# 
# qh <- (exp(r*delta) - dh)/(uh - dh)
# ql <- (exp(r*delta) - dl)/(ul - dl)
# 
# tpm <- matrix(c(.5, .5, .1, .9), nrow = 2, byrow = TRUE)
# G <- transition2Generator(P = tpm)
# #P <- diag(c(1,1)) + delta*G
# P <- expm(delta*G)
# 
# 
# ph <- P[1,1]
# pl <- P[2,2]
# #df <- 1
# df = exp(-r*delta)
# 
# x <- rep(0, n)
# x[1] <- 2
# for(i in 2:n) {
#   x[i] <- sample(1:2,1, prob=P[x[i-1],])
# }
# 
# s <- rep(0,n)
# s[1] <- S0
# 
# for(i in 2:n) {
#   if (x[i-1] ==1) s[i] <- s[i-1]*sample(c(uh,dh), 1, prob = c(qh, 1-qh))
#   if (x[i-1] ==2) s[i] <- s[i-1]*sample(c(ul,dl), 1, prob = c(ql, 1-ql))
# }
# x
# plot(s, type = "line")
# 
# ##############################################################################
# # How many times do we get regime change
# n <- 1e2
# mciter <- 1e3 
# z <- 1
# y <- 1
# Time <- 1.0
# delta <- Time/n
# 
# r <- .05
# 
# S0 <- 100
# K = 105
# 
# sig_high <-0.4
# sig_low <- 0.1
# #sig <- 0.3
# 
# uh <- exp(sig_high*sqrt(delta))
# ul <- exp(sig_low*sqrt(delta))
# 
# dh <- exp(-sig_high*sqrt(delta))
# dl <- exp(-sig_low*sqrt(delta))
# 
# qh <- (exp(r*delta) - dh)/(uh - dh)
# ql <- (exp(r*delta) - dl)/(ul - dl)
# 
# tpm <- matrix(c(.5, .5, .1, .9), nrow = 2, byrow = TRUE)
# G <- transition2Generator(P = tpm)
# #P <- diag(c(1,1)) + delta*G
# P <- expm(delta*G)
# 
# 
# ph <- P[1,1]
# pl <- P[2,2]
# #df <- 1
# df = exp(-r*delta)
# 
# 
# 
# tr <- vector("logical", mciter)
# for(i in 1:mciter) {
#   x <- rep(0, n)
#   x[1] <- 2
#   for(j in 2:n) {
#     x[j] <- sample(1:2,1, prob=P[x[j-1],])
#   }
# 
#   tr[i] <- sum(x!=2)!=0
# }
#   
# sum(tr)
low_MC
low_MCg

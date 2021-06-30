# # Pricing European Call option as done in Section 4.2 of ADM
# if (!require(markovchain)) install.packages('markovchain')
# if (!require(expm)) install.packages('expm')
# if (!require(fOptions)) install.packages('fOptions')
# if (!require(dplyr)) install.packages('dplyr')
# 
# if (!require(markovchain)) library(markovchain)
# if (!require(expm)) library(expm)
# if (!require(fOptions)) library(fOptions)
# if (!require(dplyr)) library(dplyr)


n <- 176
f <- 1 # +1 for European -1 for American
z <- 1 # +1 for Call option -1 for Put Option 
Time <- 0.25
delta <- Time/n

r <- .05

S0 <- 100
MF <- S0*exp(r*Time)
K <- 90


sig_high <- 0.4
sig_low <- 0.1
#sig <- 0.3

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

###################################################################

i <- n
high_leaf <- list()

for (j in 1:(i+1)) {
  high_leaf[[j]] <- 0.5*(z*(S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE) - K) + 
                           abs(z*(S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE) - K)))
}

low_leaf <- list()

for (j in 1:(i+1)) {
  low_leaf[[j]] <- 0.5*(z*(S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE) - K) + 
                          abs(z*(S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE) - K)))
}

stock_leaf <- list()

for (j in 1:(i+1)) {
  stock_leaf[[j]] <- S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE)
}


biglist_high <- list()
biglist_low <- list()
biglist_stock <- list()


biglist_high[[1]] <- high_leaf
biglist_low[[1]] <- low_leaf
biglist_stock[[1]] <- stock_leaf

####################################################################

for (k in 1:(n-1)) {
  
  #convert k to a more familiar variable i
  i <- n-k+1
  
  #create parent lists for this iteration
  high_parent <- biglist_high[[k]]
  low_parent <- biglist_low[[k]]
  stock_parent <- biglist_stock[[k]]
  
  #creating child lists for this iteration
  i <- i-1
  high_child <- list()
  for (j in 1:(i+1)) {
    #list_of_matrices[[j]] <- matrix(rep(99, j*(i+1-(j-1))), nrow = j)
    high_child[[j]] <- 0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE)
  }
  
  low_child <- list()
  for (j in 1:(i+1)) {
    #list_of_matrices[[j]] <- matrix(rep(99, j*(i+1-(j-1))), nrow = j)
    low_child[[j]] <- 0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE)
  }
  
  stock_child <- list()
  for (j in 1:(i+1)) {
    #list_of_matrices[[j]] <- matrix(rep(99, j*(i+1-(j-1))), nrow = j)
    stock_child[[j]] <- S0*diag(uh^(0:(j-1))*dh^((j-1):0))%*%matrix(rep(ul^(0:(i-j + 2 -1))*dl^((i-j + 2 -1):0) , j),  nrow =j, byrow=TRUE)
  }
  
  i <- i+1
  
  #populating the child lists from the parent lists
  
  for (j in 1:(i+1)){
    for (x in 1:j){
      for(y in 1:(2+i-j)){
        #high_parent[[j]][x,y]
        if(j!=1&x!=1) high_child[[j-1]][x-1,y] <- high_child[[j-1]][x-1,y] + df*ph*qh*high_parent[[j]][x,y]
        if(j!=1&x!=j) high_child[[j-1]][x,y] <- high_child[[j-1]][x,y] + df*ph*(1-qh)*high_parent[[j]][x,y]
        if(j!=(i+1)&y!=1) low_child[[j]][x,y-1] <- low_child[[j]][x,y-1] + df*(1-pl)*(ql)*high_parent[[j]][x,y]
        if(j!=(i+1)&y!=(2+i-j)) low_child[[j]][x,y] <- low_child[[j]][x,y] + df*(1-pl)*(1-ql)*high_parent[[j]][x,y]
        #low_parent[[j]][x,y]
        if(j!=(i+1)&y!=1) low_child[[j]][x,y-1] <- low_child[[j]][x,y-1] + df*(pl)*(ql)*low_parent[[j]][x,y]
        if(j!=(i+1)&y!=(2+i-j)) low_child[[j]][x,y] <- low_child[[j]][x,y] + df*(pl)*(1-ql)*low_parent[[j]][x,y]
        if(j!=1&x!=1) high_child[[j-1]][x-1,y] <- high_child[[j-1]][x-1,y] + df*(1-ph)*(qh)*low_parent[[j]][x,y]
        if(j!=1&x!=j) high_child[[j-1]][x,y] <- high_child[[j-1]][x,y] + df*(1-ph)*(1-qh)*low_parent[[j]][x,y]
      }
    }
  }
  
  #Finding American payoffs for the child nodes
  if(f == -1){
    i <- i - 1
    for (j in 1:(i+1)){
      for (x in 1:j){
        for(y in 1:(2+i-j)){
          high_child[[j]][x,y] <- max(high_child[[j]][x,y], z*(stock_child[[j]][x,y]-K))
          low_child[[j]][x,y] <- max(low_child[[j]][x,y], z*(stock_child[[j]][x,y]-K))
        }
      }
    }
    i <- i + 1
  }
  
  biglist_high[[k+1]] <- high_child
  biglist_low[[k+1]] <- low_child
  biglist_stock[[k+1]] <- stock_child
}


####################################################################

option_high <- df*ph*(sum(as.vector(biglist_high[[n]][[2]])*c(1-qh, qh))) + df*(1-ph)*(sum(as.vector(biglist_low[[n]][[2]])*c(1-qh, qh)))

option_low <- df*(1-pl)*(sum(as.vector(biglist_high[[n]][[1]])*c(1-ql, ql))) + df*(pl)*(sum(as.vector(biglist_low[[n]][[1]])*c(1-ql, ql)))

if (z == 1) type <- "c"
if (z== -1) type <- "p"
if (f == 1) flavour <- "e"
if (f== -1) flavour <- "a"
# type %>% cat("\n")
# flavour %>% cat("\n")
if(f==1) GBSOption(TypeFlag = type, S = S0, X = K, Time = Time, r = r, b = r, sigma = sig_low)@price %>% cat("GBS price is",., "\n")
# CRRBinomialTreeOption(TypeFlag = paste(type,flavour,sep=""), S = S0, X = K, Time = Time, r = r, b = r, sigma = sig_low, n = n)@price -> crr 
# crr %>% cat("CRR price is",., "\n")
# if (z == 1) cat("RHS =", option_low + K*exp(-r*Time), "\n")
# if (z== -1) cat("LHS =", S0 + option_low, "or", S0 + crr, "\n")
cat(type, flavour,"ADM price (high) for n =", n, "is :", option_high, "\n")
cat(type, flavour,"ADM price (low) for n =", n, "is :", option_low, "\n")


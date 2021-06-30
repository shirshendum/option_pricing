# # Pricing European Call option as done in Section 4.2 of ADM

# library(fOptions)
# install.packages("markovchain")
# library(markovchain)
# library(expm)

ivs_high <- vector()
ivs_low <- vector()

for(fac in c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)) {

n <- 25
y <- 1 # +1 for European -1 for American
z <- 1 # +1 for Call option -1 for Put Option 
Time <- 0.25


delta <- Time/n
r <- 0.05
S0 <- 100
K = fac*S0*exp(r*Time)

sig_high <-0.5
sig_low <- 0.2
#sig <- 0.3

uh <- exp(sig_high*sqrt(delta))
ul <- exp(sig_low*sqrt(delta))

dh <- exp(-sig_high*sqrt(delta))
dl <- exp(-sig_low*sqrt(delta))

qh <- (exp(r*delta) - dh)/(uh - dh)
ql <- (exp(r*delta) - dl)/(ul - dl)
###################################################################

tpm <- matrix(c(0.5, 0.5, 0.1, 0.9), nrow = 2, byrow = TRUE)
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
  if(y == -1){
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

#cat("price (high) for n =", n, "is : ", option_high, "\n")
#cat("price (low) for n =", n, "is : ", option_low, "\n")

iv_high <- GBSVolatility(price = option_high, TypeFlag = "c", S = 100, X = K, Time = Time, r = r, b = r)
iv_low <- GBSVolatility(price = option_low, TypeFlag = "c", S = 100, X = K, Time = Time, r = r, b = r)

ivs_high <- c(ivs_high, iv_high)
ivs_low <- c(ivs_low, iv_low)

}

df <- data.frame(high= ivs_high, low = ivs_low, 
                 x = rep(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3), 4), 
                 Maturity = as.character(c(rep(1,7), rep(0.75,7), rep(0.5, 7), rep(0.25, 7))))


ggplot(df, aes(x=x, y=low, group=Maturity)) +
  geom_line(aes(color=Maturity))+
  geom_point(aes(color=Maturity)) +
  ggtitle("Initial Volatility Low") + 
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Strike ATM factor", y = "Implied Volatility")


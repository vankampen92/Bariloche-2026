#Set the seed
set.seed(579423)

# Function that follows the stochastic evolution of our system
DifferentIslands <- function(n0, mu0, lambda, p, K, t){
  
  # gn <- rep(lambda, K) # It's constant so, I'll use lambda.
  mup <- mu0 * (p ^ (1:K - 1))
  rn <- mup * n0
  
  out <- matrix(NA, ncol = K + 1, nrow = length(t))
  out[1, ] <- c(t[1], n0)
  
  n <- n0 # We need to keep track of the populations before events happen.
  
  # Now, we assign "tickets". Each island can suffer only one event at a time. 
  # As the islands are independent, we can make a list of the events that are 
  # going to happen. For that, we make two vectors, one with the time at which 
  # the event is going to happen (according to the sum of the rates in the 
  # island) and another with the population after the event (event sampled with
  # the probabilities of that particular island). These are our tickets.
  
  dt <- rep(NA, K)
  for(i in 1:K){
    dt[i] <- t[1] + rexp(1, lambda + rn[i])
    n0[i] <- n0[i] + sample(c(1, -1), 1, prob = c(lambda, rn[i]))
  }
    
  for (j in 2:length(t)){
    #Now we select the ticket that happened first, make n change to n0, and draw
    #another ticket.
    while(min(dt) < t[j]){ 
      index <- which.min(dt)
      n[index] <- n0[index]
       
      rn[index] <- mup[index] * n0[index]
      dt[index] <- dt[index] + rexp(1, lambda + rn[index])
      n0[index] <- n0[index] + sample(c(1, -1), 1, prob = c(lambda, rn[index]))
    }
    
  out[j, ] <- c(t[j], n) # We sample.
  }
  
  colnames(out) <- c("Time", 1:K)
  out
}

# Now we set the  parameters for case i
t <- seq(0, 12, 0.2) # A vector indicating when we sample
lambda <- 10 # Inmigration
mu0 <- 1.0 
K <- 100 # Number of islands
p <- 0.99
n0 <- rep(0, K) # Initial population size


# Now we simulate one stochastic realization
sim <- DifferentIslands(n0, mu0, lambda, p, K, t)

plot(sim[, 1], rowSums(sim[, -1]), type = "l", ylab = "# Species", 
     xlab = "Time", main = "Species in 100 islands")

plot(sim[, c(1, 51)], type = "l", ylab = "# Species", xlab = "Time", 
     main = "Island 50")

# Parameters for case ii
t <- seq(0, 12, 0.2) 
lambda <- 10
mu0 <- 1.0
K <- 10000
p <- 0.9999
n0 <- rep(0, K) 


# Simulation
sim2 <- DifferentIslands(n0, mu0, lambda, p, K, t)

plot(sim2[, 1], rowSums(sim2[, -1]), type = "l", ylab = "# Species", 
     xlab = "Time", main = "Species in 10000 islands")

plot(sim2[, c(1, 5001)], type = "l", ylab = "# Species", xlab = "Time", 
     main = "Island 5000")

##### Let's see if the distribution of the number of species matches the
##### expectation.

TotalNumber <- matrix(NA, ncol = 11, nrow = 300) ### To store simulations
t <- c(0, seq(13, 18, 0.5))
lambda <- 10
mu0 <- 1.0
K <- 100
p <- 0.99
n0 <- rep(0, K) 

# Simulations
for(i in 1:300){
  TotalNumber[i,] <- rowSums(DifferentIslands(n0, mu0, lambda, p, K, t)[-1, -1])
}

hist(c(TotalNumber), breaks = 30, probability = T, xlab="Total # Species",
     main = "Stationary distribution", col = "lightgrey")

# Adding the poisson with the total rate.
TotalRate <- sum(lambda / mu0 / p^(1:K -1))
lines(1550:1850, dpois(1550:1850, TotalRate), col = "red", lty = 8)

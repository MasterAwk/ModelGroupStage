############## Data Explorations ##################
library(ggplot2)
### Import data
ECL <- read.csv("raw_ECL.csv", stringsAsFactors=FALSE)
ECL$Score_diff <- ECL$Home_score - ECL$Away_score #Add new column

### Some explorations
# Distribution of score difference
qplot(as.numeric(ECL$Score_diff), geom="histogram", binwidth = 0.5,  
      main = "Histogram of Score Difference", xlab = "Goal_diff",  
      fill=I("blue"),  alpha=I(.8))
mean(ECL$Score_diff)
# Mean = 0.385 >0, AKA "home advantage"

# Distribution of score difference by pairs
tapply(ECL$Score_diff, list(ECL$Home_seed, ECL$Away_seed), mean)

# Distribution of score difference by home seeds
ggplot(ECL, aes(factor(Home_seed), Score_diff)) + geom_boxplot() +
  labs(x="Home_seed", title="Score Difference by Home Seeds")

### Short version of data for optimization
dat <- matrix(unlist(ECL[ ,c(4,6,9)]), ncol=3, byrow=F)

############## Model Building & Inference #################
### Construct likelihood
lik.Skellam <- function(row,params=c(miu,H,A2,A3,A4,D2,D3,D4)){
  miu <- params[1]
  H <- params[2]
  A2 <- params[3]
  A3 <- params[4]
  A4 <- params[5]
  D2 <- params[6]
  D3 <- params[7]
  D4 <- params[8]
  
  A1 <- -(A2+A3+A4)
  D1 <- -(D2+D3+D4)
  z <- row[3]
  
  if (row[1] == 1 & row[2] == 2){
    lambda1 <- exp(miu+H+A1+D2)
    lambda2 <- exp(miu+A2+D1)
  }
  else if (row[1] == 1 & row[2] == 3){
    lambda1 <- exp(miu+H+A1+D3)
    lambda2 <- exp(miu+A3+D1)
  }
  else if (row[1] == 1 & row[2] == 4){
    lambda1 <- exp(miu+H+A1+D4)
    lambda2 <- exp(miu+A4+D1)
  }
  else if (row[1] == 2 & row[2] == 1){
    lambda1 <- exp(miu+H+A2+D1)
    lambda2 <- exp(miu+A1+D2)
  }
  else if (row[1] == 2 & row[2] == 3){
    lambda1 <- exp(miu+H+A2+D3)
    lambda2 <- exp(miu+A3+D2)
  }
  else if (row[1] == 2 & row[2] == 4){
    lambda1 <- exp(miu+H+A2+D4)
    lambda2 <- exp(miu+A4+D2)
  }
  else if (row[1] == 3 & row[2] == 1){
    lambda1 <- exp(miu+H+A3+D1)
    lambda2 <- exp(miu+A1+D3)
  }
  else if (row[1] == 3 & row[2] == 2){
    lambda1 <- exp(miu+H+A3+D2)
    lambda2 <- exp(miu+A2+D3)
  }
  else if (row[1] == 3 & row[2] == 4){
    lambda1 <- exp(miu+H+A3+D4)
    lambda2 <- exp(miu+A4+D3)
  }
  else if (row[1] == 4 & row[2] == 1){
    lambda1 <- exp(miu+H+A4+D1)
    lambda2 <- exp(miu+A1+D4)
  }
  else if (row[1] == 4 & row[2] == 2){
    lambda1 <- exp(miu+H+A4+D2)
    lambda2 <- exp(miu+A2+D4)
  }
  else if (row[1] == 4 & row[2] == 3){
    lambda1 <- exp(miu+H+A4+D3)
    lambda2 <- exp(miu+A3+D4)
  }
  
  bessel <- besselI(2*sqrt(lambda1*lambda2), abs(z))
  return(log(exp(-(lambda1+lambda2)) * (lambda1/lambda2)^(z/2) * bessel))
}

NLL <- function(data,params=c(miu,H,A2,A3,A4,D2,D3,D4)){
  return(-sum(apply(data, 1, lik.Skellam, params=params)))
}

### Optimization with regard to NLL
result <- optim(par=rnorm(8,0,1), fn=NLL, data=dat, method = "BFGS",
               hessian = TRUE , control = list(trace=0))

### Find se and confidence interval for parameters
fisher_info <- solve(result$hessian)
sigma_hat <- sqrt(diag(fisher_info))
upper_par <- result$par + 1.96 * sigma_hat
lower_par <- result$par - 1.96 * sigma_hat

### Visualize attack factors
A2 <- list(value=result$par[3], lower=lower_par[3], upper=upper_par[3])
A3 <- list(value=result$par[4], lower=lower_par[4], upper=upper_par[4])
A4 <- list(value=result$par[5], lower=lower_par[5], upper=upper_par[5])

a1_val <- -(A2$value+A3$value+A4$value)
a1_interval <- c(a1_val-1.96*mean(sigma_hat[3:5]), a1_val+1.96*mean(sigma_hat[3:5]))
A1 <- list(value=a1_val, lower=a1_interval[1], upper=a1_interval[2])

plot.new()
plot.window(xlim = c(-0.5, 0.5), ylim = c(4.2, 0.8), xaxs='i')
# Add x-axis
axis(side = 1, at = c(-0.5, 0, 0.5))
# Add y-axis
axis(side = 2, at = 4:1, las = 2)
# Add confidence intervals
lines(x=c(A1$lower,A1$upper), y=c(1,1), lwd=2)
points(A1$value, 1, pch=18, cex=1.5)
lines(x=c(A2$lower,A2$upper), y=c(2,2), lwd=2)
points(A2$value, 2, pch=18, cex=1.5)
lines(x=c(A3$lower,A3$upper), y=c(3,3), lwd=2)
points(A3$value, 3, pch=18, cex=1.5)
lines(x=c(A4$lower,A4$upper), y=c(4,4), lwd=2)
points(A4$value, 4, pch=18, cex=1.5)
# Add texts
avg_goals <- (tapply(ECL$Home_score, ECL$Home_seed, mean) + tapply(ECL$Away_score, ECL$Away_seed, mean))/2
text(x=c(A1$value, A2$value, A3$value, A4$value), y=1:4, pos=1,
     labels=paste0("seed_",1:4," (",round(avg_goals,2),")"))

################## Model Assessment/Application ###################
### Assess model by simulating 1-4 pairs, comparing to observed results.
### Estimating 3-2 pairs (closest match-up)
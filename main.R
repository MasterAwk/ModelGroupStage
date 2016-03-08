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

### Visualize defense factors
D2 <- list(value=result$par[6], lower=lower_par[6], upper=upper_par[6])
D3 <- list(value=result$par[7], lower=lower_par[7], upper=upper_par[7])
D4 <- list(value=result$par[8], lower=lower_par[8], upper=upper_par[8])

d1_val <- -(D2$value+D3$value+D4$value)
d1_interval <- c(d1_val-1.96*mean(sigma_hat[6:8]), d1_val+1.96*mean(sigma_hat[6:8]))
D1 <- list(value=d1_val, lower=d1_interval[1], upper=d1_interval[2])

plot.new()
plot.window(xlim = c(-0.5, 0.53), ylim = c(4.2, 0.8), xaxs='i')
# Add x-axis
axis(side = 1, at = c(-0.5, 0, 0.5))
# Add y-axis
axis(side = 2, at = 4:1, las = 2)
# Add confidence intervals
lines(x=-c(D1$lower,D1$upper), y=c(1,1), lwd=2)
points(-D1$value, 1, pch=18, cex=1.5)
lines(x=-c(D2$lower,D2$upper), y=c(2,2), lwd=2)
points(-D2$value, 2, pch=18, cex=1.5)
lines(x=-c(D3$lower,D3$upper), y=c(3,3), lwd=2)
points(-D3$value, 3, pch=18, cex=1.5)
lines(x=-c(D4$lower,D4$upper), y=c(4,4), lwd=2)
points(-D4$value, 4, pch=18, cex=1.5)
# Add texts
avg_concede <- (tapply(ECL$Away_score, ECL$Home_seed, mean) + tapply(ECL$Home_score, ECL$Away_seed, mean))/2
text(x=-c(D1$value, D2$value, D3$value, D4$value), y=1:4, pos=1,
     labels=paste0("seed_",1:4," (",round(avg_concede,2),")"))

################## Model Assessment/Application ###################
### Make prediction function
params <- list(miu=result$par[1], H=result$par[2], A1=A1$value, A2=A2$value, 
            A3=A3$value, A4=A4$value, D1=D1$value, D2=D2$value, D3=D3$value,
            D4=D4$value)

makePred <- function(home_seed, away_seed, params.=params){
  lambda1 <- exp(params$miu + params$H + params[[home_seed+2]] + params[[away_seed+6]])
  lambda2 <- exp(params$miu + params[[away_seed+2]] + params[[home_seed+6]])
  
  return(rpois(1, lambda1) - rpois(1, lambda2))
}

### Assess model by simulating 1-4 pairs, comparing to observed results.
score_diff_1_4 <- ECL[ECL$Home_seed == 1 & ECL$Away_seed == 4, "Score_diff"]
pred_1_4 <- replicate(n=length(score_diff_1_4), makePred(1,4))

# Visualize using side-by-side barplot
dat_1_4 <- rbind(data.frame(group=1, obs=score_diff_1_4),
            data.frame(group=2, obs=pred_1_4))

ggplot(dat_1_4, aes(x=obs, fill=factor(group))) +
  geom_histogram(binwidth=1, colour="black", position="dodge") +
  scale_fill_manual("Category\n",labels = c("Observed", "Predicted"), values=c("black","red")) +
  ggtitle(expression(atop("Histogram of Goal Difference (Projected v.s. Observed)", atop(italic("home seed = 1; away seed = 4"), "")))) +
  labs(x="Goal Difference")

### Simulating 3-2 pairs (closest match-up)
score_diff_3_2 <- ECL[ECL$Home_seed == 3 & ECL$Away_seed == 2, "Score_diff"]
pred_3_2 <- replicate(n=length(score_diff_3_2), makePred(3,2))

# Visualize using side-by-side barplot
dat_3_2 <- rbind(data.frame(group=1, obs=score_diff_3_2),
                 data.frame(group=2, obs=pred_3_2))

ggplot(dat_3_2, aes(x=obs, fill=factor(group))) +
  geom_histogram(binwidth=1, colour="black",position="dodge") +
  scale_fill_manual("Category\n",labels = c("Observed", "Predicted"), values=c("black","red")) +
  ggtitle(expression(atop("Histogram of Goal Difference (Projected v.s. Observed)", atop(italic("home seed = 3; away seed = 2"), "")))) +
  labs(x="Goal Difference")

### MSE of proposed model versus benchmark
compareMSE <- function(home_seed, away_seed, iter=1000){
  score_diff <- ECL[ECL$Home_seed == home_seed & ECL$Away_seed == away_seed, "Score_diff"]
  
  pool_sk <- c()
  pool_norm <- c()
  pool_unif <- c()
  for (i in 1:iter){
    pred_sk <- replicate(n=length(score_diff), makePred(home_seed,away_seed))
    pred_norm <- rnorm(n=length(score_diff), mean=mean(score_diff), sd=sd(score_diff))
    pred_unif <- runif(n=length(score_diff), min=min(score_diff), max=max(score_diff))
    pool_sk[i] <- sum((score_diff-pred_sk)^2)
    pool_norm[i] <- sum((score_diff-pred_norm)^2)
    pool_unif[i] <- sum((score_diff-pred_unif)^2)
  }
  return(c(Skellam=mean(pool_sk), norm=mean(pool_norm), unif=mean(pool_unif)))
}

### One sample K-S test
compareKS <- function(home_seed, away_seed, iter=1000){
  score_diff <- ECL[ECL$Home_seed == home_seed & ECL$Away_seed == away_seed, "Score_diff"]
  
  pool_sk <- c()
  pool_norm <- c()
  pool_unif <- c()
  for (i in 1:iter){
    pred_sk <- replicate(n=length(score_diff), makePred(home_seed,away_seed))
    pred_norm <- rnorm(n=length(score_diff), mean=mean(score_diff), sd=sd(score_diff))
    pred_unif <- runif(n=length(score_diff), min=min(score_diff), max=max(score_diff))
    pool_sk[i] <- ks.test(score_diff, pred_sk)$p
    pool_norm[i] <- ks.test(score_diff, pred_norm)$p
    pool_unif[i] <- ks.test(score_diff, pred_unif)$p
  }
  return(c(Skellam=sum(pool_sk<0.05), norm=sum(pool_norm<0.05), unif=sum(pool_unif<0.05)))
}

#################### Discuss rules ######################
### Make simulation functions
simulateGroup <- function(params.=params){
  # This function returns a matrix storing the results of each match-up
  matchMat <- matrix(rep(NA, 16), ncol=4)
  
  for (i in 1:4){
    for (j in 1:4){
      if (i == j){next}
      matchMat[i,j] <- makePred(i, j, params.=params.)
    }
  }
  return(matchMat)
}

makeTable <- function(matchMat, winPts=3){
  # This function takes in the matrix of results, outputs the final group table.
  goalDiff <- apply(matchMat, 1, sum, na.rm=T) - apply(matchMat, 2, sum, na.rm=T)
  wins <- apply(matchMat, 1, function(x){sum(x>0, na.rm=T)}) + apply(matchMat, 2, function(x){sum(x<0, na.rm=T)})
  draws <- apply(matchMat, 1, function(x){sum(x==0, na.rm=T)}) + apply(matchMat, 2, function(x){sum(x==0, na.rm=T)})
  
  table <- data.frame(Seed=1:4, Points=winPts*wins+draws, Goal_diff=goalDiff)
  return(table[order(table$Points, decreasing=T), ])
}

# Test
makeTable(simulateGroup())

############ A win = 2 pts? 3pts? (with tie-breaker being goal_diff)
pool_2pt <- list()
pool_3pt <- list()

for (i in 1:10000){ #Pick number of simulation cycles here
  pool_2pt[[i]] <- makeTable(simulateGroup(), winPts=2)
  pool_3pt[[i]] <- makeTable(simulateGroup(), winPts=3)
}

### Count number of ties
countTies <- function(pool){
  # This function takes a list of ranking tables, and outputs #ties in this pool.
  n <- length(pool)
  tie <- c()
  
  for (i in 1:n){
    tie[i] <- length(unique(pool[[i]]$Points)) < 4
  }
  return(c(ties=sum(tie), tieRatio=sum(tie)/n))
}

countTies(pool_2pt)
# Results in a tie 41.55% of the times.
countTies(pool_3pt)
# Results in a tie 32.1% of the times.

### Qualifying conditions
qualify <- function(table){
  # This function takes in the ranking table and outputs the qualified seeds,
  # given that the tie-breaker being goal_diff.
  if (length(unique(table$Points)) == 1){
    table <- table[order(table$Goal_diff, decreasing=T), ]
    return(c(table[1,1], table[2,1]))
  }
  else if(table[2,2]==table[3,2]){
    if(table[2,3]>table[3,3]){return(c(table[1,1],table[2,1]))}
    if(table[2,3]==table[3,3]){return(c(table[1,1],table[sample(2:3,1),1]))}
    else{return(c(table[1,1],table[3,1]))}
  }
  else {return(c(table[1,1], table[2,1]))}
}

# Test
qualify(makeTable(simulateGroup()))

countQualify <- function(pool){
  # This function takes in a list of tables and returns the qualifying ratios.
  qual <- list()
  for (i in 1:length(pool)){
    qual[[i]] <- qualify(pool[[i]])
  }
  qualifyRatio <- c(seed1=sum(unlist(qual)==1)/length(pool),seed2=sum(unlist(qual)==2)/length(pool),
                    seed3=sum(unlist(qual)==3)/length(pool),seed4=sum(unlist(qual)==4)/length(pool))
  return(qualifyRatio)
}

countQualify(pool_2pt)
countQualify(pool_3pt)
# Not much difference

################ A discussion of tie-breakers
### Generate simulations
filterTies <- function(tablePool, matPool){
  # Inputs: a list of simulated ranking tables; a list of simulated match result matrices;
  # Output: lists of results with ties
  tie_ind <- (1:length(tablePool))[sapply(tablePool, function(x){length(unique(x$Points)) != 4})]
  tieTables <- list()
  tieMats <- list()
  for (i in 1:length(tie_ind)){
    tieTables[[i]] <- tablePool[[tie_ind[i]]]
    tieMats[[i]] <- matPool[[tie_ind[i]]]
  }
  
  return(list(tieTables, tieMats))
}

matPool <- list()
tablePool <- list()
for (i in 1:1000){
  matPool[[i]] <- simulateGroup()
  tablePool[[i]] <- makeTable(matPool[[i]])
}

new.matPool <- filterTies(tablePool, matPool)[[2]]
new.tablePool <- filterTies(tablePool, matPool)[[1]]

### Summarize results
countTies <- read.csv("countTies.csv") #import manually summarized data

table(countTies$tie_seeds)
sum(countTies$goal_solve/nrow(countTies))
# Goal difference solves 89.97% of the ties
sum(countTies$history_solve/nrow(countTies))
# Match-up history solves 75.22% of the ties.

table(countTies$goal_favor)
table(countTies$history_favor)

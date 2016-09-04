################################################################################
# Loading required libraries
library(mets) # Should be version 1.1.1.1 or later
library(numDeriv)
library(mvtnorm)
library(reshape)

# Loading data
data <- read.table("data.csv", header=TRUE, sep="\t")
head(data)

# Loading composite loglikelihood function and composite score function
source("functions.R")

################################################################################
# Estimation of IPCWs using the variable "group"

# Formula
formula <- Surv(age2, status2==3)~as.factor(group)

# Design matrix
X <- model.matrix(formula, data)

# Fitting Aalen's additive model
fit <- aalen(formula, data, n.sim=0, robust=0)

# Predicting cumulative effects at observed event times
Gcxp <- Cpred(fit$cum, data$age2)[,-1]

# Calculating censoring probabilities at observed event times
Gcx <- exp(-apply(Gcxp*X,1,sum))
data$pc <- Gcx

# IPCW
data[,"ipcw"] <- as.numeric(data$status2!=3)/Gcx # if sibling is censored status2==3 then weight is zero
head(data)

################################################################################
# Fitting the model and estimating variance

# Initial values
par <- c(-1,1,1)

# Optimisation
op <- nlminb(par, loglik, gradient=score, outcome="out", data=data, rel="rel", ID="ID", weights="ipcw", control=list(trace=0))

# Checking whether model has converged or not
op$convergence # 0: succesfull convergence

# Variance-covariance matrix
H <- jacobian(score, op$par, data=data, outcome="out", ID="ID", rel="rel", weights="ipcw", grad=FALSE)
U <- score(op$par, data=data, outcome="out", ID="ID", rel="rel", weights="ipcw", grad=TRUE)
V <- t(U)%*%U
vcv <- solve(H)%*%V%*%solve(H)

################################################################################
# Results
parest(par=op$par,varcovar=vcv)

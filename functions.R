################################################################################
### Composite loglikelihood function ###
loglik <- function(par, data, outcome, rel, ID, weights, grad=FALSE){

# Parameters
beta <- par[1]
A <- par[2]
C <- par[3]

# New variable
data$kombo <- paste(data[,paste(outcome,1,sep="")], data[,paste(outcome,2,sep="")], sep="")
kombo <- "kombo"

# Dividing data
d1 <- data[which(data[,rel]=="f"),]
d2 <- data[which(data[,rel]=="h"),]

# Variance-covariance matrices
sigma1 <- matrix(0, nrow=2, ncol=2)
diag(sigma1) <- A + C + 1
sigma1[1,2] <- 0.5*A + C
sigma1[2,1] <- 0.5*A + C

sigma2 <- matrix(0, nrow=2, ncol=2)
diag(sigma2) <- A + C + 1
sigma2[1,2] <- 0.25*A + C
sigma2[2,1] <- 0.25*A + C

# Derivatives of variance-covariance matrices
dS01 <- rbind(c(1,0.5,0.5,1),c(1,1,1,1))
dS02 <- rbind(c(1,0.25,0.25,1),c(1,1,1,1))

# Joint loglikelihood of proband and family member
X <- as.matrix(rep(1,4))
est <- X%*%beta
Mu <- cbind(est,est)
Y <- cbind(c(0,0,1,1),c(0,1,0,1))
U1 <- .Call("biprobit2", Mu, cbind(X,X), sigma1, dS01, Y, NULL, FALSE, TRUE, TRUE, FALSE)
U2 <- .Call("biprobit2", Mu, cbind(X,X), sigma2, dS02, Y, NULL, FALSE, TRUE, TRUE, FALSE)
kombon1 <- 1*(d1[,kombo]=="00")+2*(d1[,kombo]=="01")+3*(d1[,kombo]=="10")+4*(d1[,kombo]=="11")
kombon2 <- 1*(d2[,kombo]=="00")+2*(d2[,kombo]=="01")+3*(d2[,kombo]=="10")+4*(d2[,kombo]=="11")
llj1 <- U1$loglik[kombon1,]
llj2 <- U2$loglik[kombon2,]

# Marginal loglikelihood of proband
o1 <- log(pnorm(est[1]/sqrt(A+C+1)))
o0 <- log(1-pnorm(est[1]/sqrt(A+C+1)))
llm1 <- (d1[,paste(outcome,1,sep="")]==1)*o1 + (d1[,paste(outcome,1,sep="")]==0)*o0
llm2 <- (d2[,paste(outcome,1,sep="")]==1)*o1 + (d2[,paste(outcome,1,sep="")]==0)*o0

# Combining
loglik1 <- cbind(llj1,llm1,d1[,weights])
loglik2 <- cbind(llj2,llm2,d2[,weights])

# Conditional composite loglikelihood
d1$cll <- loglik1[,3]*(loglik1[,1] - loglik1[,2])
d2$cll <- loglik2[,3]*(loglik2[,1] - loglik2[,2])

k <- -sum(d1$cll)-sum(d2$cll)

if(grad==TRUE){
    dg <- rbind(d1,d2)
    pc1 <- melt(tapply(dg$cll, dg[,ID], sum))[,2]
    return(pc1)
}
if(grad==FALSE){
    return(k)
}
}

################################################################################
### Marginal score function ###
scoremarg <- function(X, beta, A, C, Y){
    mu <- X%*%beta
    num <- (Y-pnorm(mu/sqrt(A+C+1)))*dnorm(mu/sqrt(A+C+1))
    denum <- pnorm(mu/sqrt(A+C+1))*(1-pnorm(mu/sqrt(A+C+1)))

    outB <- X/sqrt(A+C+1)
    outsigA <- -mu/(2*sqrt(A+C+1)^3)
    outsigC <- -mu/(2*sqrt(A+C+1)^3)

    uB <- (num/denum)*outB
    usigA2 <- (num/denum)*outsigA
    usigC2 <- (num/denum)*outsigC

    umarg <- cbind(uB, usigA2, usigC2)
    return(umarg)
}

################################################################################
### Composite score function ###
score <- function(par, data, outcome, rel, ID, weights, grad=FALSE){

# Parameters
beta <- par[1]
A <- par[2]
C <- par[3]

# New variable
data$kombo <- paste(data[,paste(outcome,1,sep="")], data[,paste(outcome,2,sep="")], sep="")
kombo <- "kombo"

# Dividing data
d1 <- data[which(data[,rel]=="f"),]
d2 <- data[which(data[,rel]=="h"),]

# Specifying variance-covariance matrices
sigma1 <- matrix(0, nrow=2, ncol=2)
diag(sigma1) <- A + C + 1
sigma1[1,2] <- 0.5*A + C
sigma1[2,1] <- 0.5*A + C

sigma2 <- matrix(0, nrow=2, ncol=2)
diag(sigma2) <- A + C + 1
sigma2[1,2] <- 0.25*A + C
sigma2[2,1] <- 0.25*A + C

# The derivatives of the variance-covariance matrices with respect to A and C
dS01 <- rbind(c(1,0.5,0.5,1),c(1,1,1,1))
dS02 <- rbind(c(1,0.25,0.25,1),c(1,1,1,1))

# Joint score of proband and family member
X <- as.matrix(rep(1,4))
est <- X%*%beta
Mu <- cbind(est,est)
Y <- cbind(c(0,0,1,1),c(0,1,0,1))
U1 <- .Call("biprobit2", Mu, cbind(X,X), sigma1, dS01, Y, NULL, FALSE, TRUE, TRUE, FALSE)
U2 <- .Call("biprobit2", Mu, cbind(X,X), sigma2, dS02, Y, NULL, FALSE, TRUE, TRUE, FALSE)
kombon1 <- 1*(d1[,kombo]=="00")+2*(d1[,kombo]=="01")+3*(d1[,kombo]=="10")+4*(d1[,kombo]=="11")
kombon2 <- 1*(d2[,kombo]=="00")+2*(d2[,kombo]=="01")+3*(d2[,kombo]=="10")+4*(d2[,kombo]=="11")
sj1 <- U1$score[kombon1,]
sj2 <- U2$score[kombon2,]

# Marginal score of proband
X <- as.matrix(1,1)
Y <- c(0,1)
sm <- scoremarg(X, beta, A, C, Y)
o1 <-  1*(d1[,paste(outcome,1,sep="")]=="0")+2*(d1[,paste(outcome,1,sep="")]=="1")
o2 <-  1*(d2[,paste(outcome,1,sep="")]=="0")+2*(d2[,paste(outcome,1,sep="")]=="1")
sm1 <- sm[o1,]
sm2 <- sm[o2,]

# Combining
scoreout1 <- cbind(sj1,sm1,d1[,weights])
scoreout2 <- cbind(sj2,sm2,d2[,weights])

cscore1 <- cbind(d1[,ID],scoreout1[,7]*(scoreout1[,1]-scoreout1[,4]), scoreout1[,7]*(scoreout1[,2]-scoreout1[,5]), scoreout1[,7]*(scoreout1[,3]-scoreout1[,6]))
cscore2 <- cbind(d2[,ID],scoreout2[,7]*(scoreout2[,1]-scoreout2[,4]), scoreout2[,7]*(scoreout2[,2]-scoreout2[,5]), scoreout2[,7]*(scoreout2[,3]-scoreout2[,6]))

# Conditional composite score
cscore <- rbind(cscore1, cscore2)
m <- -cbind(colSums(cscore[,c(2,3,4)]))

if(grad==TRUE){
    pc1 <- melt(tapply(cscore[,2], cscore[,1], sum))[,2]
    pc2 <- melt(tapply(cscore[,3], cscore[,1], sum))[,2]
    pc3 <- melt(tapply(cscore[,4], cscore[,1], sum))[,2]
    pc <- cbind(pc1,pc2,pc3)
    return(pc)
}
if(grad==FALSE){
    return(m)
}
}

################################################################################
### Output ###

parest <- function(par,varcovar){
est <- estimate(NULL, function(p)list(pnorm(p[1]/sqrt(p[2]+p[3]+1)),p[2]/(p[2]+p[3]+1), p[3]/(p[2]+p[3]+1),
pmvnorm(lower=rep(-Inf,2),upper=rep(p[1],2),mean=rep(0,2),sigma=rbind(c(p[2]+p[3]+1, 0.5*p[2]+p[3]),c(0.5*p[2]+p[3],p[2]+p[3]+1)))[1]/pnorm(p[1]/sqrt(p[2]+p[3]+1)),
pmvnorm(lower=rep(-Inf,2),upper=rep(p[1],2),mean=rep(0,2),sigma=rbind(c(p[2]+p[3]+1, 0.25*p[2]+p[3]),c(0.25*p[2]+p[3],p[2]+p[3]+1)))[1]/pnorm(p[1]/sqrt(p[2]+p[3]+1)),
(pmvnorm(lower=rep(-Inf,2),upper=rep(p[1],2),mean=rep(0,2),sigma=rbind(c(p[2]+p[3]+1, 0.5*p[2]+p[3]),c(0.5*p[2]+p[3],p[2]+p[3]+1)))[1]/pnorm(p[1]/sqrt(p[2]+p[3]+1)))/pnorm(p[1]/sqrt(p[2]+p[3]+1)),
(pmvnorm(lower=rep(-Inf,2),upper=rep(p[1],2),mean=rep(0,2),sigma=rbind(c(p[2]+p[3]+1, 0.25*p[2]+p[3]),c(0.25*p[2]+p[3],p[2]+p[3]+1)))[1]/pnorm(p[1]/sqrt(p[2]+p[3]+1)))/pnorm(p[1]/sqrt(p[2]+p[3]+1))),vcov=varcovar,coef=par)$coefmat
rownames(est) <- c("Cumulative risk at 75 years","sigA^2","sigC^2","Casewise concordance full siblings","Casewise concordance half siblings","Relative recurrence risk ratio full siblings","Relative recurrence risk ratio half siblings")
esto <- cbind(round(est[,c(1:4)],4),est[,5])
colnames(esto)[5] <- "P-value"
return(esto)
}

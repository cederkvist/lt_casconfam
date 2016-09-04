################################################################################
# Preparation

library(mets)

################################################################################
# Alpha functions
alpha_C <- function(t, Tau, q.1){
    (t <= Tau)*qnorm(pnorm(q.1)*(t/Tau))+(t>Tau)*(q.1)
    }
#alpha_D <- function(t) 0.1*(t-80)
alpha_D <- function(t) 0.06*(t-80)

################################################################################
# Cumulative incidence functions
Fc <- function(r.effect, t, Tau, q.1, a2, c2){
        alpha <- alpha_C(t, Tau, q.1)
        pc <- pnorm(alpha + r.effect, mean=0, sd=sqrt(1-(a2+c2)))
        return(pc)
}

Fd <- function(r.effect, t, c2){
        alpha <- alpha_D(t)
        pc <- pnorm(alpha + r.effect, mean=0, sd=sqrt(1-c2))
        return(pc)
}

################################################################################
# Functions: Marginal at time t divided by marginal at time Tau
fmarkC <- function(r.effect, t, Tau, q.1, a2, c2){
    fnum <- Fc(r.effect, t, Tau, q.1, a2, c2)
    fden <- Fc(r.effect, Tau, Tau, q.1, a2, c2)
    fm <- fnum/fden
    return(fm)
}

fmarkD <- function(r.effect, t, Tau, c2){
    fnum <- Fd(r.effect, t, c2)
    fden <- Fd(r.effect, Tau, c2)
    fm <- fnum/fden
    return(fm)
}


################################################################################
# Function for finding time points
time.new <- function(x, Tau1, Tau2, q.1, a2, c2){
    f1 <- Fc(x[,2], Tau1, Tau1, q.1, a2, c2)
    f2 <- Fd(x[,3], Tau2, c2)
    u <- runif(dim(x)[1],0.0001)
    time <- (x[,1]==1)*pnorm(qnorm(f1*u, mean=0, sd=sqrt(1-(a2+c2)))-x[,2])*Tau1/pnorm(q.1)+
        (x[,1]==0)*((qnorm(f2*u, mean=0, sd=sqrt(1-c2))-x[,3])/0.06+80)
    t <- ifelse(time<0,0,time)
    return(t)
    }

#time <- function(x, tt1, tt2, Tau1, Tau2, q.1, a2, c2){
#    if(x[1]==1){
#        fmark.tC <- fmarkC(x[2], tt1, Tau1, q.1, a2, c2)
#        time <- approx(fmark.tC, tt1, runif(1))$y
#        return(time)
#    }
#    if(x[1]==0){
#        fmark.tD <- fmarkD(x[3], tt2, Tau2, c2)
#        time <- approx(fmark.tD, tt2, runif(1))$y
#        return(time)
#    }
#}

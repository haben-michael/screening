
rnorm.trunc <- function(n,mean,sd,lower) {u <- runif(n); mean+sd*qnorm(u+pnorm((lower-mean)/sd)*(1-u))}
dnorm.trunc <- function(x,mean,sd,lower) (x>=lower) * 1/sd * dnorm((x-mean)/sd) / (1-pnorm((lower-mean)/sd))
pnorm.trunc <- function(q,mean,sd,lower) (q>=lower)*(pnorm((q-mean)/sd) - pnorm((lower-mean)/sd)) / (1-pnorm((lower-mean)/sd))


unif.Z <- local({
    dZ <- function(z)dunif(z,-1/2,1/2)
    pZ <- function(q)punif(q,-1/2,1/2)
    qZ <- function(u) (u-1/2)*(u>0)*(u<1)
    rZ.c <- function(n,cutoff) qZ(runif(n)*(1-pZ(cutoff))+pZ(cutoff))
    theta.to.cutoff <- function(theta)(2*theta-1/2)
    F.delta.c <- function(q,cutoff) (q>=0 & q<=1/2-cutoff)*(-q^2/2/(1/2-cutoff)^2+1/2+q/(1/2-cutoff)) + (q<0 & q>=-1/2+cutoff)*(1/2*(q/(1/2-cutoff)+1)^2) + 1*(q>1/2-cutoff)
    c.prime <- function(theta)2
    theta.max <- 1/3
    supp.Z <- c(-1/2,1/2)
    return(mget(ls()))
})
normal.Z <- local({
    dZ <- dnorm
    pZ <- pnorm
    qZ <- qnorm
    rZ.c <- function(n,c) rnorm.trunc(n,0,1,c)
    mu <- function(x)exp(log(dnorm(x))-log(1-pnorm(x)))
    theta.to.cutoff <- function(theta)uniroot(function(x)mu(x)-theta,c(-1,1),extendInt='yes')$root  
    kernel.mean <- integrate(function(x)dnorm(x)^2,-Inf,Inf)$val
    theta.max <- 1
    supp.Z <- c(-Inf,Inf)
    return(mget(ls()))
})


S.env <- function(supp.S=c(0,1)) {
    return(list(
    dS=function(s)dunif(s,supp.S[1],supp.S[2]),
    pS=function(q)punif(q,supp.S[1],supp.S[2]),
    rS=function(n)runif(n,supp.S[1],supp.S[2]),
    E.S2=diff(supp.S)^2/12 + mean(supp.S)^2,
    mean.S.pair=with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
    ))
}


tau <- function(z,s) {
    y <- z/s
    theta.fe <- sum(y*s^2)/sum(s^2)
    cor(z-theta.fe*s,s,method='kendall')
}



dminnorm <- function(x,mean=c(0,0),sd=c(1,1),cov=0) {
    rho <- cov/prod(sd)
    rho.det <- sqrt(1-rho^2)
    f1 <- 1/sd[1] * dnorm((x-mean[1])/sd[1]) * pnorm(rho*(x-mean[1])/sd[1]/rho.det - (x-mean[2])/sd[2]/rho.det)
    f2 <- 1/sd[2] * dnorm((x-mean[2])/sd[2]) * pnorm(rho*(x-mean[2])/sd[2]/rho.det - (x-mean[1])/sd[1]/rho.det)
    f1 + f2
}

pminnorm <- function(q,mean=c(0,0),sd=c(1,1),cov=0) {
    q.sorted <- c(-Inf,sort(q))
    pieces <- sapply(1:(length(q.sorted)-1), function(idx) integrate(function(x)dminnorm(x,mean,sd,cov),q.sorted[idx],q.sorted[idx+1])$val)
    return(cumsum(pieces)[order(q)])
}

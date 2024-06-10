args <- commandArgs(trailingOnly=TRUE)
if(length(args)==5) {
    B <- as.numeric(args[1])
    n <- as.numeric(args[2])
    rZ.idx <- as.numeric(args[3])
    Z.param <- as.numeric(args[4])
    theta <- as.numeric(args[5])
}

## print(B)
## print(n)
## print(rZ.idx)
## print(Z.param)
## print(theta)
## source('/mnt/c/Users/haben/OneDrive - University of Massachusetts/umass/research/begg/2/misc.R')
source('misc.R')
rZ <- switch(rZ.idx,
             function(n)rt(n,df=Z.param)/sqrt(Z.param/(Z.param-2)),
             with(power.Z(Z.param),rZ),
             with(beta.Z(c(Z.param,Z.param)),rZ)
                          )
## rZ <-     function(n)rt(n,df=df)/sqrt(df/(df-2))
rS <- function(n) runif(n,1,4)
## rS <- function(n) 1+rbeta(n,10,10)
## rS <- function(n) 1+rlnorm(n,0,1)

## start <- Sys.time()
## require(parallel)


## <- function(B,n,rZ,rS=runif) {
stats <-     replicate(B, {
    ## z <- with(power.Z(1), rZ(n))
    ## z <- rnorm(n)
    ## z <- rt(n,df=10)
    ## z <- with(beta.Z(c(.25,.25)),rZ(n))
    z <- rZ(n)+theta
    s <- rS(n)
    y <- z/s
    v <- 1/s^2
    theta.fe <- sum(y/v)/sum(1/v)
    var.theta.fe <- 1/sum(1/v)
    ## c(egger=egger.test(y,v),begg=begg.test(y,v,exact=TRUE),ma.stat=theta.fe/sqrt(1/sum(1/v)))
    ## begg.stat <- tau.hat.pi(z,s,0) + theta.fe*D
    begg.stat <- tau.hat(z,s,theta.fe)
    begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
    c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(var.theta.fe))        
})



## filename <- paste0('/scratch/users/habnice/save',as.integer(abs(rnorm(1))*1e8),'.RData')
filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(B,n,rZ.idx,Z.param,rS,theta,stats,file=filename)


## t call:
## for x in $(seq 1 10); do for param in $(seq 2.1 .1 4); do Rscript sim.R 100 75 1 $param; done; done


## power law call:
## for x in $(seq 1 10); do for param in $(seq -.9 .1 0); do Rscript sim.R 1000 75 2 $param; done; done

## beta call
## for x in $(seq 1 10); do for param in $(seq .2 .1 1); do Rscript sim.R 100 75 3 $param; done; done


## reps=10
## B=1000
## ns=(25 75)
## thetas=(0 0.2)
## ## echo $ns

## ## students t call
## for theta in ${thetas[@]}; do for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in 2.1 3 4; do Rscript sim.R $B $n 1 $param $theta; done; done; done; done


## ## power law call:
## for theta in ${thetas[@]}; do for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in -.1 -.5 -.9; do Rscript sim.R $B $n 2 $param $theta; done; done; done; done

## ## beta call
## for theta in ${thetas[@]}; do for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in .9 .5 .1; do Rscript sim.R $B $n 3 $param $theta; done; done; done; done


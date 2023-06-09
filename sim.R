args <- commandArgs(trailingOnly=TRUE)
B <- as.numeric(args[1])
n <- as.numeric(args[2])
rZ.idx <- as.numeric(args[3])
Z.param <- as.numeric(args[4])
theta <- as.numeric(args[5])

source('misc.R')
rZ <- switch(rZ.idx,
             function(n)rt(n,df=Z.param)/sqrt(Z.param/(Z.param-2)),
             with(power.Z(Z.param),rZ),
             with(beta.Z(c(Z.param,Z.param)),rZ)
                          )
rS <- function(n)runif(n)+1

stats <-     replicate(B, {
    z <- rZ(n)+theta
    s <- rS(n)
    y <- z/s
    v <- 1/s^2
    theta.fe <- sum(y/v)/sum(1/v)
    var.theta.fe <- 1/sum(1/v)
    begg.stat <- tau.hat(z,s,theta.fe)
    begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
    c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(var.theta.fe))        
})



filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(B,n,rZ.idx,Z.param,rS,theta,stats,file=filename)



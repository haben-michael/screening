filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
pb.alphas <- c(.05,.15)
ma.alpha <- .05
powers <- lapply(filelist, function(file) {
    load(file)
    rows <- lapply(pb.alphas, function(pb.alpha) {
        null.idx <- stats['begg.pval',] > pb.alpha
        begg.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
        null.idx <- stats['egger.pval',] > pb.alpha
        egger.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
        unconditional.power <- mean(1-pnorm(abs(stats['ma.stat',])) < ma.alpha/2)
        data.frame(B=B,n=n,pb.alpha=pb.alpha,rZ.idx=rZ.idx,Z.param=Z.param,theta=theta,begg=begg.power,egger=egger.power,unconditional=unconditional.power)
    })
    do.call(rbind,rows)
})
powers <- do.call(rbind,powers)
powers <- powers[order(powers$Z.param),]
powers <- by(powers,list(powers$n,powers$pb.alpha,powers$rZ.idx,powers$Z.param,powers$theta),FUN=function(df) {
    wts <- df$B/sum(df$B)
    c(B=sum(df$B),n=unique(df$n),pb.alpha=unique(df$pb.alpha),rZ.idx=unique(df$rZ.idx),Z.param=unique(df$Z.param),theta=unique(df$theta),colSums(wts*subset(df,select=c('begg','egger','unconditional'))))
})
powers <- as.data.frame(do.call(rbind,powers))
powers <- powers[order(powers$n,powers$rZ.idx,powers$Z.param),]


## table

rZ.idxs <- 1:3
Z.params <- list(c(low=4,med=3,high=2.1), c(low=-.1,med=-.5,high=-.9), c(low=.9,med=.5,high=.1))
theta0 <- 0

require(xtable)
require(abind)

for(theta0 in c(0,.2)) {
    out <- subset(powers, theta==theta0 & ((rZ.idx==rZ.idxs[[1]] & Z.param %in% Z.params[[1]]) | (rZ.idx==rZ.idxs[[2]] & Z.param %in% Z.params[[2]]) | (rZ.idx==rZ.idxs[[3]] & Z.param %in% Z.params[[3]])), select=c(n,pb.alpha,rZ.idx,Z.param,begg,egger,unconditional))
    out <- reshape(out,varying=c('begg','egger','unconditional'),v.names='power',timevar='condition',times=c('begg','egger','unconditional'),direction='long',new.row.names=NULL)
    out$id <- NULL
    out$zeta <- out$Z.param
    for(j in rZ.idxs) for(k in 1:length(Z.params[[j]])) out$zeta[out$rZ.idx==j & out$Z.param==Z.params[[j]][k] ]  <- names(Z.params[[j]])[k]
    out$zeta <- factor(out$zeta,levels=c('low','med','high'))
    out <- out[order(out$n,out$pb.alpha,out$rZ.idx,out$condition,out$zeta),]
    dims <- sapply(subset(out,select=rev(c(n,pb.alpha,rZ.idx,condition))),function(x)sort(unique(x)))
    dims <- c(list(zeta=levels(out$zeta)),dims)
    out <- out$power
    dim(out) <- sapply(dims,length)
    dimnames(out) <- dims
    out <- round(ftable(out,row.vars=c(2,4,5),col.vars=c(3,1)),3)
    out <- xtableFtable(out,method='compact')
    sink(paste0('power_',theta0,'.tex'))
    ## attr(out,'col.vars')$rZ.idx <- c('t','power','beta')
    attr(out,'col.vars')$rZ.idx <- c('t','Power','Beta')
    attr(out,'row.vars')$condition <- unname(sapply(attr(out,'row.vars')$condition ,function(str)paste0(toupper(substr(str,1,1)),substr(str,2,nchar(str)))))
    names(attr(out,'col.vars')) <- c('$f_Z$','$\\zeta$')
    names(attr(out,'row.vars'))[names(attr(out,'row.vars'))=='condition'] <- 'Condition'
    names(attr(out,'row.vars'))[names(attr(out,'row.vars'))=='pb.alpha'] <- '$\\alpha_0$'
    print.xtableFtable(out,floating=FALSE,latex.environment=NULL)
    sink()
}









reps <- 10
B <- 1000
ns <- c(25, 75)
thetas <- c(0, 0.2)

## students t call
with(list(rZ.idx=1), for(theta in thetas) for(n in ns) for(x in reps) for(Z.param in c(2.1,3,4)) source('sim.R',local=TRUE) )

## power law call:
with(list(rZ.idx=2), for(theta in thetas) for(n in ns) for(x in reps) for(Z.param in c(-.1,-.5,-.9)) source('sim.R',local=TRUE) )

## beta call
with(list(rZ.idx=3), for(theta in thetas) for(n in ns) for(x in reps) for(Z.param in c(.9,.5,.1)) source('sim.R',local=TRUE) )

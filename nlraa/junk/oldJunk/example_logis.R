## Simulating data from the beta growth function
## First load required packaged
setwd("~/Dropbox/NonLinearAgron/nlraa/R")
source("logis.R")
library(nlme)
##library(lme4)
##library(ggplot2)

x <- seq(1,100)

y <- logis(x, 20, 50, 5)

xyplot(y ~ x, type = 'l')

## Simulated data
x <- seq(1,100, length.out=8)

nreps <- 4
ef <- 20
set.seed(1234)
ntrt <- 2
trts <- letters[1:ntrt]

prms  <- c(20,50,6)
trt.df <- matrix(c(0,0,0,-5,10,-3), byrow=TRUE,nrow=2)
k <- 0
m <- 1

bdt <- data.frame(trt=NA, rep=NA, x=rep(x,ntrt*nreps), y=NA)  

for(i in 1:nreps){

  pv <- rnorm(3,0,2)
  for(j in trts){
    y <- logis(x, prms[1] + trt.df[m,1] + pv[1], prms[2] + trt.df[m,2] + pv[2], prms[3] + trt.df[m,3] + pv[3])
    err <- rnorm(length(x), 0, y/ef)
    sq <- (1 + length(x)*k):(length(x)*(k+1))
    bdt[sq,"trt"] <- j
    bdt[sq,"rep"] <- i
    bdt[sq,"y"] <- y + err
    k <- k + 1
    m <- m + 1
  }
  m <- 1
}

bdt$trt <- as.factor(bdt$trt)

write.csv(bdt, file="logGrowthSim.csv", row.names=FALSE)

png("./figs/logistic_trta.png")
xyplot(y ~ x, data = bdt, pch = 1, col="black", subset = trt == "a")
dev.off()

bdt.a <- subset(bdt, trt == "a")

fna <- nls(y ~ logis(x, Asym, xmid, scal), data = bdt.a, start = c(Asym=20, xmid=50, scal=4))

png("./figs/logistic_trtapb.png")
xyplot(y ~ x | trt, data = bdt, pch = 1, col="black")
dev.off()


xyplot(y ~ x | trt, type=c("p"),
       cex=1.5,
       pch=c(14:17),
       col="black",
       groups = rep,
       data = bdt,
       key=list(text=list(as.character(1:nreps)), pch=c(14:17),col="black",
                     points=TRUE))

## although we might label our reps 1, 2, 3 they are not the same experimental unit
bdt <- read.csv("../data/logGrowthSim.csv")
bdt$subject <- with(bdt, trt:factor(rep))

bdt.G <- groupedData(y ~ x | subject, data = bdt)
plot(bdt.G)

fit.nlist <- nlsList(y ~ SSlogis(x, Asym, xmid, scal), data = bdt.G)

plot(fit.nlist) ## some evidence of increasing variance with fitted values

plot(intervals(fit.nlist)) ## Evidence of variability in the three parameters

fm0 <- nlme(fit.nlist)
intervals(fm0)
pairs(ranef(fm0))
pairs(ranef(fm0))
## Let's add treatment effects
plot(augPred(fm0, level = 0:1))

fe <- fixef(fm0)
fm1 <- update(fm0, fixed = list(Asym + xmid + scal ~ trt),
                    start = c(fe[1],0,fe[2],0,fe[3],0))
fm1
anova(fm1)

plot(fm1) ## Still evidence of increasing variance
intervals(fm1) 

fm2 <- update(fm1, weights = varPower(form=~x))

plot(fm2)

## The model is now over parameterized
fm3 <- update(fm2, random = pdDiag(list(Asym + xmid + scal ~ 1)))

anova(fm3)




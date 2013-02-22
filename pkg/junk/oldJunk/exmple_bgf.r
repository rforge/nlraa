## Simulating data from the beta growth function
## First load required packaged
source("SSbgf.R")
library(nlme)
##library(lme4)
##library(ggplot2)

x <- seq(1,100)

y <- bgf(x, 20, 80, 45)

xyplot(y ~ x, type = 'l')

## Simulated data
x <- seq(1,70, length.out=7)

nreps <- 4
ef <- 12
set.seed(1234)
ntrt <- 2
trts <- letters[1:ntrt]

prms  <- c(20,60,30)
trt.df <- matrix(c(0,0,0,-5,10,-10), byrow=TRUE,nrow=2)
k <- 0
m <- 1

bdt <- data.frame(trt=NA, rep=NA, x=rep(x,ntrt*nreps), y=NA)  

for(i in 1:nreps){

  pv <- rnorm(3,0,2)
  for(j in trts){
    y <- bgf(x, prms[1] + trt.df[m,1] + pv[1], prms[2] + trt.df[m,2] + pv[2], prms[3] + trt.df[m,3] + pv[3])
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

png("./figs/betagrowth_trta.png")
xyplot(y ~ x, data = bdt, pch = 1, col="black", subset = trt == "a")
dev.off()

bdt.a <- subset(bdt, trt == "a")

fna <- nls(y ~ bgf(x, w.max, t.e, t.m), data = bdt, start = c(w.max=20, t.e=80, t.m=40))

png("./figs/betagrowth_trtapb.png")
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
bdt$subject <- with(bdt, trt:factor(rep))

bdt.G <- groupedData(y ~ x | subject, data = bdt)
plot(bdt.G)

fit.nlist <- nlsList(y ~ SSbgf(x, w.max, t.e, t.m), data = bdt.G)

fit.nlist <- nlsList(y ~ SSlogis(x, Asym, xmid, scal), data = bdt.G)

plot(fit.nlist) ## some evidence of increasing variance with fitted values

plot(intervals(fit.nlist)) ## Evidence of variability in the three parameters

fm0 <- nlme(fit.nlist)
intervals(fm0)
pairs(ranef(fm0))
## Let's simplify the model
fm2 <- update(fm1, random = pdDiag(w.max + t.e + t.m ~ 1))
anova(fm1, fm2)
## Even simpler
fm3 <- update(fm2, random = pdDiag(w.max + t.m ~ 1))
anova(fm1, fm2, fm3)

fe <- fixef(fm3)
fm4 <- update(fm3, fixed = list(w.max + t.e + t.m ~ trt),
                    start = c(fe[1],0,fe[2],0,fe[3],0))
fm4
## Variance for t.m is very small
fm5 <- update(fm4, random = list(w.max ~ 1))
plot(fm5) ## Still evidence of increasing variance
intervals(fm5) 

fm6 <- update(fm5, weights = varPower())


fm6 <- update(fm5, weights = varPower(0.1, form = ~x), verbose=TRUE,
              control=list(minScale=1e-316))


fm6 <- update(fm5, weights = varPower(value=0.15, form=~x),
              control = list(pnlsTol = 0.01), verbose=TRUE)

fm6 <- update(fm5, weights = varPower(form=~x),
              control = list(pnlsTol = 0.1), verbose=TRUE)

fm6 <- update(fm5, weights = varPower(form=~x),
              control = list(pnlsTol = 4), verbose=TRUE)

fm6 <- update(fm5, weights = varPower(),
              control = list(pnlsTol = 4), verbose=TRUE)

fm6 <- update(fm5, weights = varExp(form=~x), verbose=TRUE, control = list(pnlsTol = 0.01))
fm6 <- update(fm5, weights = varExp(), verbose=TRUE, control = list(pnlsTol = 0.01))
fm6 <- update(fm5, weights = varConstPower(power=0.15), verbose=TRUE)
fm6 <- update(fm5, weights = varPower(form=~x), verbose=TRUE)
fm6 <- update(fm5, weights = varFunc(form=~x), verbose=TRUE)
fm6 <- update(fm5, weights = varPower(form=~x), control = list(pnlsTol=0.1))


anova(fit4.nlme, fit5.nlme)

plot(augPred(fit4.nlme, level = 0:1))
fit.nlme <- nlme(fit.nlist, weights = varPower())
fit.nlme <- nlme(fit.nlist, random = pdDiag(w.max + t.e + t.m ~ 1))


## Trying the nlmer function first
fit.nlmer <- nlmer(y ~ SSbgf(x,w.max, t.e,t.m) ~ w.max + t.e + t.m |subject, data = yc,
                   start = c(w.max=20, t.e=80, t.m=40))

fit.nlmer <- nlmer(y ~ SSbgf(x,w.max, t.e,t.m) ~ w.max |subject, data = yc,
                   start = c(w.max=20, t.e=80, t.m=40))

yc.b3 <- subset(yc , subject=="b:3")
yc.b2 <- subset(yc , subject=="b:2")
yc.b1 <- subset(yc , subject=="b:1")

xyplot(y + y2 ~ x, data = yc.b2 , type =c("p","l"), distribute.type=TRUE)

fit.nls.b2 <- nls(y ~ SSbgf(x, w.max, t.e, t.m), data = yc.b2, trace = TRUE)
##start = c(w.max = 20, t.e=80, t.m=45))

## JUNK

qplot(data = yc, x = x, y = y, color = factor(rep), 
      facets = . ~ trt, cex=2)

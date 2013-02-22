## Read in data from Sotiris Archontoulis
library(nlme)
source("SSbgf.R")
sm <- read.csv("../data/MaizeSorghumGreece2009.csv", skip = 8, na.strings=".")
head(sm)
doy <- c(135, 165, 173, 180, 185, 200, 222, 240, 263, 282)

datc <- NULL

for(i in 1:nrow(sm)){
  
  block <- sm[i,1]
  input <- sm[i,2]
  crop <- sm[i,3]
  yield <- c(t(sm[i,4:ncol(sm)]))
  
  dat <- data.frame(block=block, input=input, crop=crop, doy=doy, yield=yield)
  datc <- rbind(datc, dat)
}

head(datc)

xyplot(yield ~ doy | factor(block)*factor(input) , 
       groups = crop, data = datc, type = 'o',
       auto.key=TRUE)

datc$eu <- with(datc, factor(block):factor(input):factor(crop))

datc.G <- groupedData(yield ~ doy | eu, data = datc)
plot(datc.G)

datc.G.nmd <- na.omit(datc.G)
##datc.G.nmd[datc.G.nmd$doy == 135,5] <- 0.1
datc.G.nmd <- subset(datc.G.nmd, doy != 135)

plot(datc.G.nmd)

fm.lisB <- nlsList(yield ~ SSbgf(doy, w.max, t.e, t.m), data = datc.G.nmd)
fm.lisB2 <- nlsList(yield ~ bgf2(doy, w.max, w.b=0, t.e, t.m, t.b),
                    data = datc.G.nmd,
                    start=c(w.max=20,t.e=260,t.m=225,t.b=0))
fm.lisG <- nlsList(yield ~ SSgompertz(doy, Asym, b1, b2), data = datc.G.nmd)
fm.lisL <- nlsList(yield ~ SSlogis(doy, Asym, xmid, scal), data = datc.G.nmd)
fm.lisW <- nlsList(yield ~ SSweibull(doy, Asym, Drop, lrc, pwr), data = datc.G.nmd)

fmB <- nlme(fm.lisB, control = list(minScale=1e-50, pnlsTol = 0.01))
#fmB2 <- nlme(fm.lisB2, control = list(minScale=1e-50, pnlsTol = 0.01))
fmG <- nlme(fm.lisG, control = list(minScale=1e-50, pnlsTol = 0.01), verbose=TRUE)
fmL <- nlme(fm.lisL)
fmW <- nlme(fm.lisW)

## At this stage the gompertz and weibull do not converge

## The logistic seems more promising
fmL

fm0
plot(augPred(fm0, level = 0:1))
intervals(fm0)
plot(ranef(fm0))
pairs(ranef(fm0))
## Let's choose a simpler covariance structure since the model is overparameterized now
fm1 <- update(fm0, random = pdDiag(w.max + t.e + t.m ~ 1), 
              control = list(minScale=1e-300, pnlsTol = 0.1))
## Let's introduce the fixed effect of crop
fe <- fixef(fm1)
fm2 <- update(fm1, fixed = list(w.max + t.e + t.m ~ crop), 
              start = c(fe[1],0,fe[2],0,fe[3],0))

## The standard deviations for t.m and t.e are now very small we can exclude them
fm3 <- update(fm2, random = list(w.max ~ 1))

anova(fm2, fm3) ## The models are the same, simpler is better

## What about the effect of the input
fe2 <- fixef(fm3)
fm4 <- update(fm3, fixed = list(w.max + t.m + t.e ~ crop + input),
              start = c(fe2[1:2], 0, fe2[3:4], 0, fe2[5:6], 0))

plot(augPred(fm3, level = 0:1))

## Ok let's try to include the interaction first
fe3 <- fixef(fm4)
fm5 <- update(fm4, fixed = list(w.max + t.m + t.e ~ crop + input + crop:input),
              start = c(fe3[1:3], 0, fe3[4:6], 0, fe3[7:9], 0),
              control = list(minScale = 1e-310, pnlsTol=0.9), verbose=TRUE)

## The random effect of w.max is virtually zero. 
## However, the residuals show unequal variance
## We need to fit a gnls model

fm5 <- gnls(yield ~ SSbgf(doy, w.max, t.e, t.m),
            params = list(w.max + t.e + t.m ~ crop),
            start = c(w.max=31.4,0,t.e=266.1,0,t.m=236.8,0),
            data = datc.G.nmd, weights = varPower())

cfs <- fm5$coefficients
fm6 <- gnls(yield ~ SSbgf(doy, w.max, t.e, t.m),
            params = list(w.max + t.e + t.m ~ crop + input),
            start = c(cfs[1:2],0,cfs[3:4],0,cfs[5:6],0),
            data = datc.G.nmd, weights = varPower(),
            control = list(nlsTol=0.02, minScale = 1e-70))




## Adding the effect of the interaction
fe3 <- fixef(fm4)
fm5 <- update(fm4, fixed = list(w.max + t.m + t.e ~ crop + input + crop:input),
              start = c(fe3[1:3], 0, fe3[4:6], 0, fe3[7:9], 0),
              control = list(pnlsTol = 0.1, minScale = 1e-300))






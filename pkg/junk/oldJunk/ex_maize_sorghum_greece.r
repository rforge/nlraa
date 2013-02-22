## Read in data from Sotiris Archontoulis
library(nlme)
library(ggplot2)
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
plot(augPred(fmL, level = 0:1))
plot(augPred(fmB, level = 0:1))
anova(fmL, fmB)


## The logistic seems more promising
fmL ## Not very strong evidence to simplify the model yet

fe <- fixef(fmL)
fmL2 <- update(fmL, fixed = list(Asym + xmid + scal ~ crop),
               start = c(fe[1],0,fe[2],0,fe[3],0))

anova(fmL2)
## Strongest evidence of crop effect on Asymptote
## Let's add the effect of input
fe2 <- fixef(fmL2)
fmL3 <- update(fmL2, fixed = list(Asym + xmid + scal ~ crop + input),
               start = c(fe2[1:2],0,fe2[3:4],0,fe2[5:6],0))

## Can we also try to include the interaction
fe3 <- fixef(fmL3)
fmL4 <- update(fmL3, fixed = list(Asym + xmid + scal ~ crop + input + crop:input),
               start = c(fe3[1:3],0,fe3[4:6],0,fe3[7:9],0))

## A simpler model
fe4 <- fixef(fmL4)
fmL5 <- update(fmL4, fixed = list(Asym ~ crop + input + crop:input, xmid + scal ~ 1),
               start = c(fe4[1:5],fe4[9]))

## address unequal variances
plot(fmL5)
## Do we really need this?
fmL6 <- update(fmL5, weights = varPower())
anova(fmL5, fmL6)

plot(fmL6)

## What about autocorrelation?
plot(ACF(fmL6), alpha = 0.05)
## Not a very strong autocorrelation function
fmL7 <- update(fmL6, correlation = corAR1())

## Model fmL7 does not seem to be needed
anova(fmL6, fmL7)

## Can probably eliminate random effect for scal
fmL8 <- update(fmL6, random = pdDiag(Asym + xmid ~ 1))

plot(augPred(fmL8, level = 0:1))

## Very small random effects
## Try the gnls model
fmL9 <- gnls(yield ~ SSlogis(doy, Asym, xmid, scal),
             data = datc.G.nmd,
             params = list(Asym ~ crop + input + crop:input, xmid + scal ~ 1),
             weights = varPower(),
             start = fixef(fmL8))

fmB1 <- gnls(yield ~ SSbgf(doy, w.max, t.e, t.m),
             data = datc.G.nmd,
             weights = varPower(),
             control = gnlsControl(minScale=1e-10))

## Examine the model
plot(fmL9) ## Residuals look good

## plot observed and predicted
## preds <- predict(fmL9,
##                  newdata = expand.grid(doy=unique(datc.G.nmd$doy),
##                    crop = unique(datc.G.nmd$crop),
##                    input = c(1,2)))

newdat <- datc.G.nmd
newdat$preds <- fitted(fmL9)

xyplot(yield ~ preds, data = newdat)

xyplot(yield + preds ~ doy | factor(input) , groups = crop,
       data = newdat,
       type = c("p","a"), distribute.type=TRUE)

xyplot(yield ~ doy | factor(input) , groups = crop,
       data = datc.G.nmd,
       distribute.type=TRUE)

f1 <- ggplot(data = newdat, aes(y = yield, x = doy, color = crop)) +
      facet_grid(. ~ input) +
      geom_line(aes(x=doy, y = preds, color=crop)) +
      geom_point(aes(shape = crop, fill=crop)) +
      scale_shape_manual(values=c(24,21)) +
      scale_fill_manual(values = c("white","black")) +
      theme_bw()
print(f1)




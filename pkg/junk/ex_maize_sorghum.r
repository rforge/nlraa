## Load libraries
library(nlraa)
library(grid)

## Read in the data
## sm <- read.csv("../data/CornSorghumAllValues2008.csv", skip=6)
data(sm)

sm$inputHL <- as.factor(ifelse(sm$Input < 1.5, "Low","High"))
sm$input <- NULL
## Graph for plot  
p1 <- ggplot(data = sm, aes(y = Yield, x = DOY)) +
  facet_grid(. ~ inputHL) +
  geom_point(aes(fill=Crop, shape=Crop), size=2) +
  scale_shape_manual(values=c(24,21,1)) +
  scale_fill_manual(values = c("grey","black","black")) +
  scale_x_continuous("Day of the Year") +
  scale_y_continuous("Biomass (Mg/ha)") +
  theme_bw()
print(p1)
ggsave(p1, file="./figs/MaizeSorghum.pdf", height=5)


sm$eu <- with(sm, factor(Block):factor(inputHL):factor(Crop))

sm2 <- subset(sm, DOY != 141)
smG <- groupedData(Yield ~ DOY | eu, data = sm2)

fit.lis <- nlsList(Yield ~ SSbgf(DOY, w.max, t.e, t.m), data = smG)

pdf("./figs/nlslisResid.pdf", height=5, width=4)
plot(fit.lis, col="black", main="Eq. 2.5", ylim=c(-4,4))
dev.off()

fit.me <- nlme(fit.lis, control = list(minScale =1e-50, pnlsTol=0.01))

plot(augPred(fit.me, level = 0:1))

## Trying the bgf2

fit.lis2 <- nlsList(Yield ~ bgf2(DOY, w.max, w.b = 0, t.e, t.m, t.b = 141),
                    data = smG,
                    start = c(w.max = 30, t.e=280, t.m=240))

pdf("./figs/nlslisResid_revised.pdf", height=5, width=4)
plot(fit.lis2, main="Eq. 2.11", ylim=c(-4,4))
dev.off()
# to merge


fit.me2 <- nlme(fit.lis2)
fit2.me2 <- update(fit.me2, random = pdDiag(w.max + t.e + t.m ~ 1))

anova(fit.me, fit.me2)

anova(fit.me2, fit2.me2)

fe <- fixef(fit2.me2) ## Some starting values with visual help
fit3.me2 <- update(fit2.me2, fixed = list(w.max + t.e + t.m ~ Crop),
                  start = c(fe[1], -10, 20, fe[2], -40, 0, fe[3], -40, 0))

fe2 <- fixef(fit3.me2)
fit4.me2 <- update(fit3.me2, fixed = list(w.max + t.e + t.m ~ Crop + inputHL),
                  start = c(fe2[1:3], 0, fe2[4:6], 0, fe2[7:9], 0))

fe3 <- fixef(fit4.me2)
fit5.me2 <- update(fit4.me2, fixed = list(w.max + t.e + t.m ~ Crop + inputHL + Crop:inputHL),
                  start = c(fe3[1:4], 0, 0,
                            fe3[5:8], 0, 0,
                            fe3[9:12], 0, 0))

fit6.me2 <- update(fit5.me2, weights = varPower(form = ~ fitted(.) | Crop))

fit7.me2 <- update(fit6.me2, weights = varPower(form = ~ fitted(.)))

anova(fit6.me2, fit7.me2)

##plot(augPred(fit7.me2, level =0:1), ylim = c(0,60))

## Random effects are almost zero
fit8.me2 <- gnls(Yield ~ bgf2(DOY, w.max, t.e, t.m, w.b=0, t.b=141),
                 data = smG,
                 params = list(w.max + t.e + t.m ~ Crop + inputHL + Crop:inputHL),
                 weights = varPower(form = ~ fitted(.) | Crop),
                 start = fixef(fit7.me2))

anova(fit6.me2, fit8.me2)

fit9.me2 <- gnls(Yield ~ bgf2(DOY, w.max, t.e, t.m, w.b=0, t.b=141),
                 data = smG,
                 params = list(w.max + t.e + t.m ~ Crop + inputHL + Crop:inputHL),
                 weights = varPower(form = ~ fitted(.)),
                 start = fixef(fit7.me2))

anova(fit8.me2, fit9.me2)
## The more complex model is favored
anova(fit8.me2)

## Simplifying the model
fit10.me2 <- gnls(Yield ~ bgf2(DOY, w.max, t.e, t.m, w.b=0, t.b=141),
                 data = smG,
                 params = list(w.max + t.e ~ Crop + inputHL + Crop:inputHL, t.m ~ Crop),
                 weights = varPower(form = ~ fitted(.) | Crop),
                 start = fixef(fit7.me2)[1:15])

anova(fit8.me2, fit10.me2)

datfs <- data.frame(obs = smG$Yield, fit = fitted(fit8.me2))
datfs.low <- datfs[datfs$obs < 10,]
datfs.med <- datfs[datfs$obs > 10 & datfs$obs < 20,]
datfs.high <- datfs[datfs$obs > 20,]

(mse.low <- rmse(datfs.low$obs, datfs.low$fit))
(mse.med <- rmse(datfs.med$obs, datfs.med$fit))
(mse.high <- rmse(datfs.high$obs, datfs.high$fit))

## Some statistics
mse <- rmse(smG$Yield, fitted(fit10.me2))
cc <- concor(smG$Yield, fitted(fit10.me2))
mef <- ef(smG$Yield, fitted(fit10.me2))


plot(fit10.me2)
## plot(ACF(fit8.me2), alpha = 0.05)

smG$prds <- fitted(fit8.me2)

## Graph
smGd <- as.data.frame(smG)


p2 <- qplot(data = smGd, x = Yield, y = prds, xlab="Observed biomass (Mg/ha)",
            ylab = "Simulated biomass (Mg/ha)", xlim=c(0,55), ylim=c(0,55)) +
  geom_abline(intercept=0,slope=1) +
  theme_bw() +
  geom_text(aes(label = paste("low. rMSE = ",round(mse.low$rmse,1)," (Mg/ha)", sep=""), x = 11, y = 47)) +
  geom_text(aes(label = paste("med. rMSE = ",round(mse.med$rmse,1)," (Mg/ha)", sep=""), x = 11, y = 49.5)) +
  geom_text(aes(label = paste("high. rMSE = ",round(mse.high$rmse,1)," (Mg/ha)", sep=""), x = 11, y = 52))
print(p2)
ggsave(p2, file="./figs/ModAgg.pdf")


doys <- 168:303
ndat <- expand.grid(DOY=doys, Crop= unique(smG$Crop), inputHL=unique(smG$inputHL))
ndat$preds <- predict(fit8.me2, newdata = ndat)

ndat2 <- ndat
ndat2[ndat2$Crop == "M" & ndat2$DOY > 270,"preds"] <- NA
ndat2 <- na.omit(ndat2)

p3 <- ggplot(data = smG, aes(y = Yield, x = DOY)) +
  facet_grid(. ~ inputHL) +
##  geom_point(aes(fill=Crop, shape=Crop), size=2) +
  stat_summary(fun.data="mean_cl_boot", aes(y = Yield, x = DOY, shape=Crop, fill=Crop)) +
  geom_line(aes(x = DOY, y = preds, linetype = Crop), data=ndat2) +
  scale_shape_manual(values=c(24,21,1)) +
  scale_fill_manual(values = c("grey","black","black")) +
  scale_x_continuous("Day of the Year") +
  scale_y_continuous("Biomass (Mg/ha)") +
  theme_bw() + 
  theme(legend.position = c(0.1,0.8),
       legend.key.size = unit(1,"cm"),
       legend.text = element_text(size=12))
print(p3)

ggsave(p3, file="./figs/SorghumMaizePredsCI.pdf", height=5)

p4 <- ggplot(data = smG, aes(y = Yield, x = DOY)) +
  facet_grid(. ~ inputHL) +
  geom_point(aes(fill=Crop, shape=Crop), size=2) +
  geom_line(aes(x = DOY, y = preds, linetype = Crop), data=ndat2) +
  scale_shape_manual(values=c(24,21,1)) +
  scale_fill_manual(values = c("grey","black","black")) +
  scale_x_continuous("Day of the Year") +
  scale_y_continuous("Biomass (Mg/ha)") +
  theme_bw() +
  theme(legend.position = c(0.1,0.8),
       legend.key.size = unit(1, "cm"),
       legend.text = element_text(size=12))
print(p4)

ggsave(p4, file="./figs/SorghumMaizePredsRaw.png", dpi=300)

## Maximum biomass and the point at which it was reached.
ndat[which.max(ndat$preds),]
ndat[which.max(ndat[ndat$Crop == "M" & ndat$input == "High",]$preds),]


## Testing the bgf2
library(nlraa)
x <- 1:300
y <- bgf2(x, 20, 2, 260, 220, 0)
plot(x, y)

library(coenocliner)

# Gradients
xy <- cbind(x = seq(from = 4, to = 7, length.out = 100),
            y = seq(from = 1, to = 100, length.out = 100))

# Species responses
opt <- c(4,5,6)
tol <- rep(0.25, 3)
h <- c(10,20,30)
parm <- cbind(opt = opt, tol = tol, h = h)

opty <- c(25, 50, 75)
tol <- c(5, 10, 20)
pars <- list(px = parm,
             py = cbind(opt = opty, tol = tol))




M <- 20                                    # number of species
ming <- 3.5                                # gradient minimum...
maxg <- 7                                  # ...and maximum
locs <- seq(ming, maxg, length = 100)      # gradient locations
opt  <- runif(M, min = ming, max = maxg)   # species optima
tol  <- rep(0.25, M)                       # species tolerances
h    <- ceiling(rlnorm(M, meanlog = 3))    # max abundances
pars <- cbind(opt = opt, tol = tol, h = h) # put in a matrix


# coenocline implied by these parameters
mu <- coenocline(locs, responseModel = "gaussian", params = pars,
                 expectation = TRUE)
matplot(locs, mu, lty = "solid", type = "l", xlab = "pH", ylab = "Abundance")

# negative binomial count data
simnb <- coenocline(locs, responseModel = "gaussian", params = pars,
                    countModel = "negbin", countParams = list(alpha = 0.5))

# Realistic species data
matplot(locs, simnb, lty = "solid", type = "p", pch = 1:10, cex = 0.8,
        xlab = "pH", ylab = "Abundance")

head(simnb)

# Species counts
species = data.frame(simnb[,])

# mvabund input
species.mva = mvabund(species)

m1 = manyglm(species.mva ~ locs, family = "negative.binomial")
plot(m1, which = c(1:4))

m2 = manyglm(species.mva ~ locs + I(locs^2), family = "negative.binomial")
anova(m1,m2, nBoot = 100)

m2.summary = summary(m2, nBoot = 100, test="LR", p.uni = "adjusted")
m2.summary$uni.p

# differential observation depth will be random poisson noise
ObsDiff = rpois(100, lambda = 2)*1000+500
hist(ObsDiff)

# A part of the bias is added only if the species is observed in a sample.
species2 = as.data.frame(matrix(0, nrow=100,ncol=20))
colnames(species2) = colnames(species)

# The replacement solution from here:
# http://stackoverflow.com/questions/8214303/conditional-replacement-of-values-in-a-data-frame
species2[,1][species[,1] != 0] <- species[,1][species[,1] != 0] + ObsDiff[species[,1] != 0]/20

cbind(species$X1, species2$X1, ObsDiff/20)
plot(species$X1, species2$X1)

# replace all
for (i in 1:20) {
  species2[,i][species[,i] != 0] <- species[,i][species[,i] != 0] + ObsDiff[species[,i] != 0]/20
}

ObsCounts = apply(species2,1,sum)
plot(ObsCounts, ObsDiff)

# some more controls
cbind(species$X11, species2$X11, ObsDiff/20)
plot(species$X11, species2$X11)
  
# Models
sp2.mva = mvabund(species2)

mc1 = manyglm(sp2.mva ~ locs + I(locs^2), family = "negative.binomial")
mc1.summary = summary(mc1, nBoot = 50, test="LR", p.uni = "adjusted")

# Coefficients 
plot(coef(m2)["locs",], coef(mc1)["locs",])
plot(coef(m2)["I(locs^2)",], coef(mc1)["I(locs^2)",])
cbind(coef(m2)["locs",],coef(mc1)["locs",])

mc2 = manyglm(sp2.mva ~ ObsCounts + locs + I(locs^2), family = "negative.binomial")
anova(mc1, mc2, nBoot=50)
mc2.summary = summary(mc2, nBoot = 50, test="LR", p.uni = "adjusted")

plot(coef(m2)["locs",], coef(mc2)["locs",])
plot(coef(m2)["I(locs^2)",], coef(mc2)["I(locs^2)",])

                      
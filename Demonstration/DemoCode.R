####################################
# Dm: manipulation of the mediator #
# M: measurement of the mediator   #
# X: independent variable          #
# Y: dependent variable            #
####################################


# IV Mediation Illustration
library('haven')
library('janitor')
library('lavaan')
library('mediation')

# reading data
rawdata = read_sav('/Users/luzhiming/000我/sysu/自己/主要/EM/importantApplied/illustration/Illu2.sav')
rawdata = data.frame(rawdata)

# effect coding and standardization
rawdata$Gender = rawdata$Gender-1.5; rawdata$Perf = rawdata$Perf-1.5
data = data.frame(X = rawdata$Gender, Dm = rawdata$Perf, M = rawdata$CEOperf, Y = rawdata$Jobfit)
data$M = scale(data$M); data$Y = scale(data$Y)



## ImaiMed
m1.med.std = lm(M~X*Dm, data = data)
m2.med.std = lm(Y~X*M, data = data)
fit.med.std = mediate(model.m = m1.med.std, model.y = m2.med.std, treat = "X", mediator = "M",
                      control.value = -0.5, treat.value = 0.5)
summary(fit.med.std)



## NON-IV 
model0 = 'M~a*X+d*Dm
        Y~cp*X+b*M+g*X:M'

fit0 = sem(model=model0, data=data)
summary(fit0)  # b.est = 0.18, Z-value = 2.24, p = 0.025

PNIE.est0 = coef(fit0)['a']*(coef(fit0)['b']-0.5*coef(fit0)['g'])
TNIE.est0 = coef(fit0)['a']*(coef(fit0)['b']+0.5*coef(fit0)['g'])

# estimating the standard error of the indirect effect using delta method

covs0 = vcov(fit0)
IE.vcov0 = matrix(c(covs0['a','a'], covs0['a','b'], covs0['a','g'],
                    covs0['a','b'], covs0['b','b'], covs0['b','g'],
                    covs0['a','g'], covs0['b','g'], covs0['g','g']),3,3)

delta.f.PNIE0 = matrix(c(coef(fit0)['b']-0.5*coef(fit0)['g'], coef(fit0)['a'], -0.5*coef(fit0)['a']), 3, 1)
delta.f.TNIE0 = matrix(c(coef(fit0)['b']+0.5*coef(fit0)['g'], coef(fit0)['a'], 0.5*coef(fit0)['a']), 3, 1)

PNIE.SEest0 = sqrt(t(delta.f.PNIE0) %*% IE.vcov0 %*% delta.f.PNIE0)
TNIE.SEest0 = sqrt(t(delta.f.TNIE0) %*% IE.vcov0 %*% delta.f.TNIE0)

PNIE.Z.est0 = abs(PNIE.est0)/PNIE.SEest0; PNIE.p0 = pnorm(abs(PNIE.Z.est0))
TNIE.Z.est0 = abs(TNIE.est0)/TNIE.SEest0; TNIE.p0 = pnorm(abs(TNIE.Z.est0))


## IV mediation
model1 = 'M~a*X+d*Dm
        Y~cp*X+b*M+g*X:M
        M~~Y'

fit1 = sem(model=model1, data=data)
summary(fit1)  # b.est = 0.25, Z-value = 2.89, p = 0.004.

PNIE.est = coef(fit1)['a']*(coef(fit1)['b']-0.5*coef(fit1)['g'])
TNIE.est = coef(fit1)['a']*(coef(fit1)['b']+0.5*coef(fit1)['g'])

# estimating the standard error of the indirect effect using delta method

covs = vcov(fit1)
IE.vcov = matrix(c(covs['a','a'], covs['a','b'], covs['a','g'],
                   covs['a','b'], covs['b','b'], covs['b','g'],
                   covs['a','g'], covs['b','g'], covs['g','g']),3,3)

delta.f.PNIE = matrix(c(coef(fit1)['b']-0.5*coef(fit1)['g'], coef(fit1)['a'], -0.5*coef(fit1)['a']), 3, 1)
delta.f.TNIE = matrix(c(coef(fit1)['b']+0.5*coef(fit1)['g'], coef(fit1)['a'], 0.5*coef(fit1)['a']), 3, 1)

PNIE.SEest = sqrt(t(delta.f.PNIE) %*% IE.vcov %*% delta.f.PNIE)
TNIE.SEest = sqrt(t(delta.f.TNIE) %*% IE.vcov %*% delta.f.TNIE)

PNIE.Z.est = abs(PNIE.est)/PNIE.SEest; PNIE.p = pnorm(abs(PNIE.Z.est))
TNIE.Z.est = abs(TNIE.est)/TNIE.SEest; TNIE.p = pnorm(abs(TNIE.Z.est))


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


## IV mediation
model1 = 'M~a*X+d*Dm
         Y~cp*X+b*M
         M~~Y'

fit1 = sem(model=model1, data=data)
summary(fit1)  # b.est = 0.24, Z-value = 2.74

IE.est = coef(fit1)['a']*coef(fit1)['b']

# estimating the standard error of the indirect effect using delta method
covs = vcov(fit1)
IE.se.est = sqrt(coef(fit1)['b']^2*covs['a','a']+
                   coef(fit1)['a']^2*covs['b','b']+
                   2*coef(fit1)['a']*coef(fit1)['b']*covs['a','b'])

Z.est = abs(IE.est)/IE.se.est

SIG.IVmed = Z.est > qnorm(1-.025)#not significant
pnorm(Z.est)

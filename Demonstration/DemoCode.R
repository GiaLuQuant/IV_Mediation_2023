#################################
# Dm: manipulation of the mediator 
# M: measurement of the mediator
# X: independent variable 
# Y: dependent variable 



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



## non-IV SEM
model0 = 'M~a*X+d*Dm
         Y~cp*X+b*M'

fit0 = sem(model=model0, data=data)
summary(fit0)  # b.est = 0.171, Z-value = 2.086 
IE.est0 = coef(fit0)['a']*coef(fit0)['b']

covs0 = vcov(fit0)
IE.se.est0 = sqrt(coef(fit0)['b']^2*covs0['a','a']+
                   coef(fit0)['a']^2*covs0['b','b']+
                   2*coef(fit0)['a']*coef(fit0)['b']*covs0['a','b'])

Z.est0 = abs(IE.est0)/IE.se.est0
SIG.nonIV = Z.est0 > qnorm(1-.025) #not significant



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

SIG.IVmed = Z.est0 > qnorm(1-.025)#not significant

###############################################################################
#  X: the independent variable
#  M: the measurement of the mediator
#  Dm: the manipulation of the mediator
#  Y: the dependent variable
#  V: the confounder representing the potentially omitted variables

gen_data <- function(N,a,b,cp,d,g,v,NT,p=.5){
  # preset variances
  v.X = v.Dm = p*(1-p); v.M = v.Y = v.V = 1
  
  # computing variances of error terms and interaction terms
  v.eM = v.M - a^2*v.X - d^2*v.Dm - v^2*v.V
  
  v.XM = d^2*v.X*v.Dm + v^2*v.X*v.V + v.X*v.eM
  
  v.eY = v.Y - cp^2*v.X - b^2*v.M - g^2*v.XM - v^2*v.V 
  
  # covariance matrix of the error terms
  Ve = matrix(c(v.V, 0, 0,
                0, v.eM, 0,
                0, 0, v.eY), 3,3,byrow = T)
  colnames(Ve)=rownames(Ve) = c('eV','eM','eY')
  
  errors = mvrnorm(N, c(0,0,0), Ve)
  
  if(!NT){
    errors = apply(errors, 2, tras_ADF, skew=-2, kurtosis = 6)
  }
  
  X = rbinom(N, 1, 0.5)-0.5; Dm = rbinom(N, 1, 0.5)-0.5
  
  V = errors[,'eV']
  M = a*X+d*Dm+v*V+errors[,'eM']
  Y = cp*X+b*M+g*X*M+v*V+errors[,'eY']
  
  dat = data.frame(X=X, Dm=Dm, M=M, Y=Y, V=V)  
  
  return(dat)
}

fn = function(pars, skew, kurtosis){
  b = pars[1]; c = pars[2]; d = pars[3]
  eq1 = b^2 + 6*b*d + 2*c^2 + 15*d^2 -1
  eq2 = 2*c* (b^2 + 24*b*d + 105*d^2 + 2) - skew
  eq3 = 24*(b*d + c^2*(1+b^2+28*b*d) + d^2*(12 + 48*b*d + 141*c^2 + 225*d^2)) - kurtosis
  v = c(eq1, eq2, eq3)
  Q = sum(v^2)
  return(Q)
}

tras_ADF <- function(error, skew, kurtosis){
  bcd = nlm(fn, c(1,1,1), skew, kurtosis)$estimate
  b = bcd[1]; c=bcd[2]; d=bcd[3]; a=-c
  e.ADF  = a + b*error+ c*error^2 + d*error^3
  return(e.ADF)
}

# proposed method
IVSEM_main <- function(data){
  model = 'M~a*X+d*Dm
           Y~cp*X+b*M+g*X:M
           M~~Y'
  fit = try(sem(model=model, data=data), silent=T)
  if(inherits(fit,"try-error")){
    res = rep(NA,16)
  }else{
    if(fit@optim[["converged"]]==F){
      res = rep(NA,16)
    }else{
      para.est = c(coef(fit)['a'], coef(fit)['b'], coef(fit)['cp'], 
                   coef(fit)['d'],coef(fit)['g'])
      sum = summary(fit)
      pValues = sum$pe$pvalue[1:5]
      ## need further revision
      resV.M = sum$pe$est[7]; resV.Y = sum$pe$est[8]
      covs = vcov(fit)
      if(det(covs)<0){
        sigma.est=rep(NA,5); IE.Xn.SEest = IE.Xp.SEest = NA; infoDef=-1
      }else{
        sigma.est = c(covs['a','a'],covs['b','b'],covs['cp','cp'],
                      covs['d','d'],covs['g','g'])
        
        IE.vcov = matrix(c(covs['a','a'], covs['a','b'], covs['a','g'],
                           covs['a','b'], covs['b','b'], covs['b','g'],
                           covs['a','g'], covs['b','g'], covs['g','g']),3,3)
        
        delta.f.Xn = matrix(c(coef(fit)['b']-0.5*coef(fit)['g'], coef(fit)['a'], -0.5*coef(fit)['a']), 3, 1)
        delta.f.Xp = matrix(c(coef(fit)['b']+0.5*coef(fit)['g'], coef(fit)['a'], 0.5*coef(fit)['a']), 3, 1)
        
        IE.Xn.SEest = sqrt(t(delta.f.Xn) %*% IE.vcov %*% delta.f.Xn)
        IE.Xp.SEest = sqrt(t(delta.f.Xp) %*% IE.vcov %*% delta.f.Xp)
        infoDef=1
      }
      res = c(para.est, sigma.est, pValues, IE.Xn.SEest, IE.Xp.SEest, resV.M, resV.Y, infoDef)
    }
  }
  return(res)
}

TradMed <- function(data){
  m1 = lm(M~X+Dm,data)
  m2 = lm(Y~X*M,data)
  fit = try(mediate(m1,m2,treat='X',mediator='M', sims=1000, treat.value = 0.5, control.value = -0.5), silent=T)
  if (inherits(fit, 'try-error')){
    res = rep(NA,8)
  }else{
    sum=summary(fit)
    IE = c(sum$d0, sum$d1)
    IE.CI = c(sum$d0.ci, sum$d1.ci)
    IE.p = c(sum$d0.p, sum$d1.p)
    res = c(IE, IE.CI, IE.p)
  }
  return(res)
}



AnaRES <- function(nrep, N,a,b,cp,d,g,v,NT,p=.5){

  res.IVSEM = matrix(NA,nrep,20)
  res.ImaiMed = matrix(NA,nrep,8); 
  
  colnames(res.IVSEM) <- c('a.est','b.est','cp.est','d.est','g.est',
                           'a.sigma','b.sigma','cp.sigma','d.sigma','g.sigma',
                           'a.p','d.p','cp.p','b.p','g.p','PNIE.SEest','TNIE.SEest','resV.M','resV.Y','InfoDef')
  colnames(res.ImaiMed) <- c('PNIE','TNIE','PNIE.lower','PNIE.upper',
                             'TNIE.lower','TNIE.upper','PNIE.p','TNIE.p')
  
  sd.X = sd.Dm = sqrt(p*(1-p))
  a.u = a/sd.X; cp.u = cp/sd.X; d.u = d/sd.Dm; g.u = g/sd.X
  
  for(i in 1:nrep){
    data = try(gen_data(N,a.u,b,cp.u,d.u,g.u,v,NT,p=.5))
    
    if(inherits(data, 'try-error')){
      res.IVSEM[i,] = rep(NA,20); 
      res.ImaiMed[i,] = rep(NA,8);
    }else{
      # the proposed method
      res.IVSEM[i,] = IVSEM_main(data)
      
      # mediation package (measurement of mediation)  -- Imai, 2010
      res.ImaiMed[i,] = TradMed(data)
      
    }   
  }
  res = list(IVSEM=res.IVSEM, ImaiMed=res.ImaiMed)
  return(res)
}

## performance measures
# estimation bias
cal.eBias <- function(theta.hat, theta){
  if(theta==0){
    Ebias <- mean(theta.hat, na.rm=T)-theta
  }else{Ebias <- (mean(theta.hat, na.rm=T)-theta)/theta}
  return(Ebias)
}

# coverage rate
cal.CR <- function(SE, theta, theta.hat){
  cv.error = qnorm(1-0.05/2) * SE
  l = theta.hat-cv.error  # lower bound
  u = theta.hat+cv.error  # upper bound
  lcriterion = (theta >= l)
  ucriterion = (theta <= u)
  CR = mean(lcriterion*ucriterion, na.rm=T)
  return(CR)
}

# rejection rate
cal.RR <- function(Est, seEst, para){
  Z.est = abs(Est)/seEst
  Z.cri = qnorm(1-.025)
  RR = mean(Z.est>=Z.cri, na.rm=T)
  
  if(para == 0){
    power = NA
    typeIerror = RR
  }else{
    power = RR
    typeIerror = NA
  }
  return(c(power,typeIerror))
}

sep.RR <- function(RR, para){
  if(para == 0){
    power = NA
    typeIerror = RR
  }else{
    power = RR
    typeIerror = NA
  }
  return(c(power,typeIerror))
}

# deleting results with erroneous variance estimation 
del_vcov <- function(resi){
  range = 5:8
  for(i in 1:nrow(resi)){
    if(is.na(resi[i,'resV.M'])){
      resi[i,] = rep(NA, ncol(resi))
    }else if(resi[i,'resV.M']<0|resi[i,'resV.Y']<0){
      resi[i,] = rep(NA, ncol(resi))
    }
  }
  resi = remove_empty(resi, which='rows')
  for(i in 1:nrow(resi)){
    for(j in range){
      if(is.na(resi[i,j])){
        resi[i,]=rep(NA, ncol(resi))
      }
    }
  }
  resi = remove_empty(resi, which='rows')
  return(resi)
}

IVSEM.main.reg <- function(res, a, b, cp, d,g){
  
  ConvRate = 1-mean(is.na(res[,'a.est']))
  res = remove_empty(res, which='rows')
  
  res = del_vcov(res)
  infoDef = mean(res[,'InfoDef'], na.rm = T)
  
  ## performance measures for the proposed method
  # point Estimate Bias
  a.ebias = cal.eBias(res[,'a.est'], a)
  b.ebias = cal.eBias(res[,'b.est'], b)
  cp.ebias = cal.eBias(res[,'cp.est'], cp)
  d.ebias = cal.eBias(res[,'d.est'], d)
  g.ebias = cal.eBias(res[,'g.est'], g)
  
  # se para
  a.se = sd(res[,'a.est'])
  b.se = sd(res[,'b.est'])
  cp.se = sd(res[,'cp.est'])
  d.se = sd(res[,'d.est'])
  g.se = sd(res[,'g.est'])
  
  # se bias
  a.sebias = cal.eBias(sqrt(res[,'a.sigma']), a.se)
  b.sebias = cal.eBias(sqrt(res[,'b.sigma']), b.se)
  cp.sebias = cal.eBias(sqrt(res[,'cp.sigma']), cp.se)
  d.sebias = cal.eBias(sqrt(res[,'d.sigma']), d.se)
  g.sebias = cal.eBias(sqrt(res[,'g.sigma']), g.se)
  
  # coverage rate
  a.CR = cal.CR(sqrt(res[,'a.sigma']), a, res[,'a.est'])
  b.CR = cal.CR(sqrt(res[,'b.sigma']), b, res[,'b.est'])
  cp.CR = cal.CR(sqrt(res[,'cp.sigma']), cp, res[,'cp.est'])
  d.CR = cal.CR(sqrt(res[,'d.sigma']), d, res[,'d.est'])
  g.CR = cal.CR(sqrt(res[,'g.sigma']), g, res[,'g.est'])
  
  # rejection rate
  a.RR = sep.RR(mean(res[,'a.p']<.05), a)
  b.RR = sep.RR(mean(res[,'b.p']<.05), b)
  cp.RR = sep.RR(mean(res[,'cp.p']<.05), cp)
  d.RR = sep.RR(mean(res[,'d.p']<.05), d)
  g.RR = sep.RR(mean(res[,'g.p']<.05), g)
  
  reg.res = matrix(c(a.ebias, b.ebias, cp.ebias, d.ebias, g.ebias,
                     a.se, b.se, cp.se, d.se, g.se,
                     a.sebias,b.sebias,cp.sebias,d.sebias, g.sebias,
                     a.CR, b.CR, cp.CR, d.CR, g.CR,
                     a.RR, b.RR, cp.RR, d.RR, g.RR, ConvRate, infoDef), 1,32)
  colnames(reg.res) = c('a.ebias','b.ebias','cp.ebias','d.ebias', 'g.ebias',
                        'a.empSE','b.empSE','cp.empSE','d.empSE','g.empSE',
                        'a.sebias','b.sebias','cp.sebias','d.sebias','g.sebias',
                        'a.CR','b.CR','cp.CR','d.CR','g.CR',
                        'a.power','a.typeIerror','b.power','b.typeIerror',
                        'cp.power','cp.typeIerror','d.power','d.typeIerror',
                        'g.power','g.typeIerror','ConvRate', 'NegDef')
  
  return(reg.res)
}

IVSEM.main.IE <- function(res, a, b,g){
  ConvRate = 1-mean(is.na(res[,'a.est']))
  res = remove_empty(res, which='rows')
  
  res = del_vcov(res)
  infoDef = mean(res[,'InfoDef'], na.rm = T)

  PNIE.tv=a*(b-0.5*g)
  TNIE.tv=a*(b+0.5*g)
  
  PNIE.est=res[,'a.est']*(res[,'b.est']-0.5*res[,'g.est'])
  TNIE.est=res[,'a.est']*(res[,'b.est']+0.5*res[,'g.est'])
  
  PNIE.ebias=cal.eBias(PNIE.est, PNIE.tv)
  TNIE.ebias=cal.eBias(TNIE.est, TNIE.tv)
  
  PNIE.empSE = sd(PNIE.est); TNIE.empSE = sd(TNIE.est)
  
  PNIE.sebias = cal.eBias(res[,'PNIE.SEest'], PNIE.empSE)
  TNIE.sebias = cal.eBias(res[,'TNIE.SEest'], TNIE.empSE)
  
  PNIE.CR = cal.CR(res[,'PNIE.SEest'], PNIE.tv, PNIE.est)
  TNIE.CR = cal.CR(res[,'TNIE.SEest'], TNIE.tv, TNIE.est)
  PNIE.RR = cal.RR(PNIE.est, res[,'PNIE.SEest'], PNIE.tv)
  TNIE.RR = cal.RR(TNIE.est, res[,'TNIE.SEest'], TNIE.tv)
  
  LMIE.res = matrix(c(PNIE.tv,TNIE.tv, PNIE.ebias, TNIE.ebias, 
                      PNIE.empSE, TNIE.empSE, PNIE.sebias, TNIE.sebias,
                      PNIE.CR, TNIE.CR, PNIE.RR, TNIE.RR),1,14)
  
  colnames(LMIE.res) = c('PNIE.para','TNIE.para','PNIE.eBias','TNIE.eBias',
                         'PNIE.empSE','TNIE.empSE','PNIE.seBias','TNIE.seBias',
                         'PNIE.CR','TNIE.CR','PNIE.power','PNIE.typeIerror', 'TNIE.power','TNIE.typeIerror')

  return(LMIE.res)
}

TradMed.res.IE <- function(res, a,b,cp,d,g){
  ## performace measures for Imai's mediation package
  ConvRate = 1-mean(is.na(res[,1]))
  res = remove_empty(res, which='rows')
  
  PNIE.tv=a*(b-0.5*g)
  TNIE.tv=a*(b+0.5*g)
  
  PNIE.ebias = cal.eBias(res[,'PNIE'], PNIE.tv)
  TNIE.ebias = cal.eBias(res[,'TNIE'], TNIE.tv)
  
  PNIE.CR = mean((res[,'PNIE.lower']<=PNIE.tv) * (res[,'PNIE.upper']>=PNIE.tv), na.rm=T)
  TNIE.CR = mean((res[,'TNIE.lower']<=TNIE.tv) * (res[,'TNIE.upper']>=TNIE.tv), na.rm=T)
  
  PNIE.RR = sep.RR(mean(res[,'PNIE.p']<=0.05, na.rm=T), PNIE.tv)
  TNIE.RR = sep.RR(mean(res[,'TNIE.p']<=0.05, na.rm=T), TNIE.tv)
  
  TradMed.res = matrix(c(PNIE.tv, TNIE.tv, PNIE.ebias, TNIE.ebias, 
                         PNIE.CR, TNIE.CR, PNIE.RR, TNIE.RR, ConvRate),1, 11)
  colnames(TradMed.res) = c('PNIE.para', 'TNIE.para', 'PNIE.eBias', 'TNIE.eBias','PNIE.CR', 'TNIE.CR', 
                            'PNIE.power','PNIE.typeIerror', 'TNIE.power','TNIE.typeIerror','ConvRate')
  
  return(TradMed.res)
}


simRES <- function(res, a, b, cp, d,g){
  IVSEM.res = IVSEM.main.reg(res[[1]], a, b, cp, d,g)
  IVSEM.IE.res = IVSEM.main.IE(res[[1]], a, b,g)
  TradMed.res = TradMed.res.IE(res[[2]],a,b,cp,d,g)
  
  final_res = list(IVSEM.res = IVSEM.res, IVSEM.IE.res = IVSEM.IE.res, 
                   TradMed.res=TradMed.res)
  return(final_res)
}



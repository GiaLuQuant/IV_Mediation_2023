###############################################################################
#  X: the independent variable
#  M: the measurement of the mediator
#  Dm: the manipulation of the mediator
#  Y: the dependent variable
#  V: the confounder representing the potentially omitted variables

gen_data <- function(N,a,b,cp,d,v,NT,p=.5){
  # preset variances
  v.X = v.Dm = p*(1-p); v.M = v.Y = v.V = 1
  
  # computing variances of error terms and interaction terms
  v.eM = v.M - a^2*v.X - d^2*v.Dm - v^2*v.V
  v.eY = v.Y - cp^2*v.X - b^2*v.M - v^2*v.V 
  
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
  Y = cp*X+b*M+v*V+errors[,'eY']
  
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
           Y~cp*X+b*M
           M~~Y'
  fit = try(sem(model=model, data=data), silent=T)
  if(inherits(fit,"try-error")){
    res = rep(NA,16)
  }else{
    if(fit@optim[["converged"]]==F){
      res = rep(NA,16)
    }else{
      para.est = c(coef(fit)['a'], coef(fit)['b'], coef(fit)['cp'], coef(fit)['d'])
      pValues = summary(fit)$pe$pvalue[1:4]
      resV.M = var(data$M) - coef(fit)['a']^2*var(data$X) - coef(fit)['d']^2*var(data$Dm)
      resV.Y = var(data$Y) - coef(fit)['cp']^2*var(data$X) - coef(fit)['b']^2*var(data$M)
      covs = try(vcov(fit), silent=T)
      if (inherits(covs, 'try-error')){sigma.est=rep(NA,4); IE.se.est=NA; infoDef = -1
      }else{if(det(covs)<0){
        sigma.est=rep(NA,4); IE.se.est=NA; infoDef = -1}
        else{
          sigma.est = c(covs['a','a'],covs['b','b'],covs['cp','cp'],covs['d','d'])
          IE.se.est = sqrt(coef(fit)['b']^2*covs['a','a']+
                             coef(fit)['a']^2*covs['b','b']+
                             2*coef(fit)['a']*coef(fit)['b']*covs['a','b'])
          infoDef = 1
        }
      }
      res = c(para.est, sigma.est, pValues, IE.se.est, resV.M, resV.Y, infoDef)
    }
  }
  return(res)
}

TradMed <- function(data){
  m1 = lm(M~X+Dm,data)
  m2 = lm(Y~X+M,data)
  fit = try(mediate(m1,m2,treat='X',mediator='M', sims=1000, treat.value = 0.5, control.value = -0.5), silent=T)
  if (inherits(fit, 'try-error')){
    res = rep(NA,4)
  }else{
    sum=summary(fit)
    IE = sum$d0 
    IE.CI = sum$d0.ci
    IE.p = sum$d0.p
    res = c(IE, IE.CI, IE.p)
  }
  return(res)
}



AnaRES <- function(nrep, N,a,b,cp,d,v,NT=TRUE,skew=0,kurtosis=0,p=.5){

  res.IVSEM = matrix(NA,nrep,16)
  res.ImaiMed = matrix(NA,nrep,4); 
  
  colnames(res.IVSEM) <- c('a.est','b.est','cp.est','d.est',
                           'a.sigma','b.sigma','cp.sigma','d.sigma',
                           'a.p','d.p','cp.p','b.p','IE.SEest','resV.M','resV.Y','InfoDef')
  colnames(res.ImaiMed) <- c('IE.est','IE.lower','IE.upper','IE.p')
  
  sd.X = sd.Dm = sqrt(p*(1-p))
  a.u = a/sd.X; cp.u = cp/sd.X; d.u = d/sd.Dm
  
  for(i in 1:nrep){
    data = try(gen_data(N,a.u,b,cp.u,d.u,v,NT,p=.5))
    
    if(inherits(data, 'try-error')){
      res.IVSEM[i,] = rep(NA,16); 
      res.ImaiMed[i,] = rep(NA,4);
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

IVSEM.main.reg <- function(res, a, b, cp, d){
  
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
  
  # se para
  a.se = sd(res[,'a.est'])
  b.se = sd(res[,'b.est'])
  cp.se = sd(res[,'cp.est'])
  d.se = sd(res[,'d.est'])
  
  # se bias
  a.sebias = cal.eBias(sqrt(res[,'a.sigma']), a.se)
  b.sebias = cal.eBias(sqrt(res[,'b.sigma']), b.se)
  cp.sebias = cal.eBias(sqrt(res[,'cp.sigma']), cp.se)
  d.sebias = cal.eBias(sqrt(res[,'d.sigma']), d.se)
  
  # coverage rate
  a.CR = cal.CR(sqrt(res[,'a.sigma']), a, res[,'a.est'])
  b.CR = cal.CR(sqrt(res[,'b.sigma']), b, res[,'b.est'])
  cp.CR = cal.CR(sqrt(res[,'cp.sigma']), cp, res[,'cp.est'])
  d.CR = cal.CR(sqrt(res[,'d.sigma']), d, res[,'d.est'])
  
  # rejection rate
  a.RR = sep.RR(mean(res[,'a.p']<.05), a)
  b.RR = sep.RR(mean(res[,'b.p']<.05), b)
  cp.RR = sep.RR(mean(res[,'cp.p']<.05), cp)
  d.RR = sep.RR(mean(res[,'d.p']<.05), d)
  
  reg.res = matrix(c(a.ebias, b.ebias, cp.ebias, d.ebias, 
                     a.se, b.se, cp.se, d.se, 
                     a.sebias,b.sebias,cp.sebias,d.sebias,
                     a.CR, b.CR, cp.CR, d.CR, 
                     a.RR, b.RR, cp.RR, d.RR, ConvRate, infoDef), 1,26)
  colnames(reg.res) = c('a.ebias','b.ebias','cp.ebias','d.ebias',
                        'a.empSE','b.empSE','cp.empSE','d.empSE',
                        'a.sebias','b.sebias','cp.sebias','d.sebias',
                        'a.CR','b.CR','cp.CR','d.CR',
                        'a.power','a.typeIerror','b.power','b.typeIerror',
                        'cp.power','cp.typeIerror','d.power','d.typeIerror','ConvRate', 'NegDef')
  
  return(reg.res)
}

IVSEM.main.IE <- function(res, a, b){
  ConvRate = 1-mean(is.na(res[,'a.est']))
  res = remove_empty(res, which='rows')
  
  res = del_vcov(res)
  infoDef = mean(res[,'InfoDef'], na.rm = T)

  IE.tv=a*b
  IE.est=res[,'a.est']*res[,'b.est']
  IE.ebias=cal.eBias(IE.est, IE.tv)
  IE.empSE = sd(IE.est)
  IE.sebias = cal.eBias(res[,'IE.SEest'], IE.empSE)
  IE.CR = cal.CR(res[,'IE.SEest'], IE.tv, IE.est)
  IE.RR = cal.RR(IE.est, res[,'IE.SEest'], IE.tv)
  
  LMIE.res = matrix(c(IE.tv, IE.ebias, IE.empSE, IE.sebias, IE.CR, IE.RR),1,7)
  
  colnames(LMIE.res) = c('IE.para','IE.eBias','IE.empSE','IE.seBias','IE.CR','IE.power','IE.typeIerror')

  return(LMIE.res)
}

TradMed.res.IE <- function(res, a,b,cp,d,g,h){
  ## performace measures for Imai's mediation package
  ConvRate = 1-mean(is.na(res[,1]))
  res = remove_empty(res, which='rows')
  
  IE.tv = a*b
  TM.IE.ebias = cal.eBias(res[,'IE.est'], IE.tv)
  TM.IE.CR = mean((res[,'IE.lower']<=IE.tv) * (res[,'IE.upper']>=IE.tv), na.rm=T)
  TM.IE.RR = sep.RR(mean(res[,'IE.p']<=0.05, na.rm=T), IE.tv)
  TradMed.res = matrix(c(IE.tv, TM.IE.ebias, TM.IE.CR, TM.IE.RR, ConvRate),1, 6)
  colnames(TradMed.res) = c('IE.para', 'IE.eBias','IE.CR','IE.power','IE.typeIerror','ConvRate')
  
  return(TradMed.res)
}


simRES <- function(res, a, b, cp, d){
  IVSEM.res = IVSEM.main.reg(res[[1]], a, b, cp, d)
  IVSEM.IE.res = IVSEM.main.IE(res[[1]], a, b)
  TradMed.res = TradMed.res.IE(res[[2]],a,b,cp,d)
  
  final_res = list(IVSEM.res = IVSEM.res, IVSEM.IE.res = IVSEM.IE.res, 
                   TradMed.res=TradMed.res)
  return(final_res)
}



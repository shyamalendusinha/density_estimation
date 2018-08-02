
rm(list = ls())

library(foreach)
library(doParallel)

#********Data Storing**************************

outputnumber = '2'

setwd("~/Desktop/finalRdir/DPMM/v12")
outputfilecsv = paste("diffq",outputnumber,".csv",sep="")
outputfile = paste("diffq",outputnumber,".Rdata",sep="")

#************Data Simulation parameters*********************
Bq = 100
q.upper.lim = 1000
q.lower.lim = 100
q.interval = 100
B = Bq*((q.upper.lim-q.lower.lim)/q.interval+1)


#****Data simulation parameters******************
#mu.i, sigma.i^2 are independently generated
meand = 0
musd = sqrt(3)
sigma2shape = 9
sigma2rate = 3

r.biv.density = function(q){
  sigma2d = rgamma(q, shape=sigma2shape, rate=sigma2rate)
  mud = rnorm(q, mean=meand, sd=musd)
  return(list(mud=mud, sigma2d=sigma2d))
}

d.biv.density = function(x, y) dnorm(x, mean=meand, sd=musd)*
  dgamma(y, shape=sigma2shape, rate=sigma2rate)

#marginal densities
d.mu.marginal.density = function(x) dnorm(x, mean=meand, sd=musd)
d.sigma2.marginal.density = function(y) dgamma(y, shape=sigma2shape, rate=sigma2rate)

# d.mu.marginal.density = function(x,mu.grid) density(x, n=length(mu.grid), from=mu.grid[1], to=mu.grid[length(mu.grid)])$y
# d.sigma2.marginal.density = function(y,sigma2.grid) density(y, n=length(sigma2.grid), from=sigma2.grid[1], to=sigma2.grid[n2])$y




source("DPMMfunction.R")
source("SUREmethods.R")
source("grouplinearfunction_all.R")
source("heteroscedastic_abhra.R")

cl = min(40,B)
registerDoParallel(cl)

alloutput_diffq = foreach(b = 1:B,.combine='cbind',.inorder=FALSE) %dopar% {
  
  set.seed(1234+47*b)
  
  n = 4
  q = ceiling(b/Bq)*q.interval+(q.lower.lim-q.interval)
  ni.vec = rep(n,q)
  n1 = 100
  n2 = 100
  
  random.sample = r.biv.density(q)
  mud = random.sample$mud
  sigma2d = random.sample$sigma2d
  errord = rnorm(n*q)
  x.vec = rep(mud,ni.vec)+sqrt(rep(sigma2d,ni.vec))*errord
  ids = rep(1:q,ni.vec)
  xbar.vec = tapply(x.vec,ids,mean)
  S2.vec = tapply(x.vec,ids,var)
  sigma2d.xbar = sigma2d/ni.vec
  S2.vec.xbar = S2.vec/ni.vec
  mu.grid = seq(min(xbar.vec)-IQR(xbar.vec), max(xbar.vec)+IQR(xbar.vec), 
                length.out=2*n1+1)[seq(2,2*n1,2)]
  sigma2.grid = seq(max(min(S2.vec)-IQR(S2.vec),0.001), max(S2.vec)+IQR(S2.vec), 
                    length.out=2*n2+1)[seq(2,2*n2,2)]
  
  biv.grid = array(0,c(n1*n2,2))
  biv.grid[,1] = rep(mu.grid,n2)
  biv.grid[,2] = rep(sigma2.grid,each=n1)
  
  # true densities on the grid
  d.biv.true = d.biv.density(biv.grid[,1], biv.grid[,2])
  d.mu.true = d.mu.marginal.density(mu.grid)
  d.sigma2.true = d.sigma2.marginal.density(sigma2.grid)
  
  d.biv.true.kde =c(kde2d(c(mud,mud),c(sigma2d,-sigma2d),n=c(n1,n2),lims=c(range(mu.grid),
                                                                           range(sigma2.grid)))$z)
  d.mu.true.kde = density(mud, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  d.sigma2.true.kde = density(c(sigma2d,-sigma2d), n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.true = sum((d.biv.true.kde-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.true = sum((d.mu.true.kde-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.true = sum((d.sigma2.true.kde-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  #naive density estimation
  avg.loss.mu.est.sample.mean = mean((mud-xbar.vec)^2)
  avg.loss.sigma2d.est.sample.var = mean((sigma2d-S2.vec)^2) 
  
  
  d.biv.sample.est.kde.est = c(kde2d(c(xbar.vec,xbar.vec),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                              range(sigma2.grid)))$z)
  d.mu.sample.est.kde.est = density(xbar.vec, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  d.sigma2.sample.est.kde.est = density(c(S2.vec,-S2.vec), n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.sample.est = sum((d.biv.sample.est.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.sample.est = sum((d.mu.sample.est.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.sample.est = sum((d.sigma2.sample.est.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  #********EBMLE***********
  
  mu.est.EBMLE = thetahat.EBMLE.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.EBMLE = mean((mud-mu.est.EBMLE)^2)
  
  d.biv.EBMLE.kde.est = c(kde2d(c(mu.est.EBMLE,mu.est.EBMLE),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                 range(sigma2.grid)))$z)
  d.mu.EBMLE.kde.est = density(mu.est.EBMLE, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.EBMLE= sum((d.biv.EBMLE.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.EBMLE= sum((d.mu.EBMLE.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #**********************
  
  
  #*******EBMOM**********
  
  mu.est.EBMOM = thetahat.EBMOM.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.EBMOM = mean((mud-mu.est.EBMOM)^2)
  
  d.biv.EBMOM.kde.est = c(kde2d(c(mu.est.EBMOM,mu.est.EBMOM),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                 range(sigma2.grid)))$z)
  d.mu.EBMOM.kde.est = density(mu.est.EBMOM, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.EBMOM= sum((d.biv.EBMOM.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.EBMOM= sum((d.mu.EBMOM.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  #*****************************
  
  #*************JS estimate*****
  
  mu.est.JS = thetahat.JS.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.JS = mean((mud-mu.est.JS)^2)
  
  d.biv.JS.kde.est = c(kde2d(c(mu.est.JS,mu.est.JS),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                        range(sigma2.grid)))$z)
  d.mu.JS.kde.est = density(mu.est.JS, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.JS= sum((d.biv.JS.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.JS= sum((d.mu.JS.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #****************************
  
  #********Oracle**************
  
  mu.est.oracle = thetahat.oracle.XKB(xbar.vec, S2.vec.xbar, mud)
  
  avg.loss.mu.est.oracle = mean((mud-mu.est.oracle)^2)
  
  d.biv.oracle.kde.est = c(kde2d(c(mu.est.oracle,mu.est.oracle),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                    range(sigma2.grid)))$z)
  d.mu.oracle.kde.est = density(mu.est.oracle, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.oracle= sum((d.biv.oracle.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.oracle= sum((d.mu.oracle.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  
  #************************************
  
  #********large class oracle**************
  
  m.or = mean(mud)
  or.loss = vector()
  
  for (i in 1:q){
    if(mud[i] < xbar.vec[i] & mud[i] > m.or){
      or.loss[i] = 0
    }
    else if(mud[i] > xbar.vec[i] & mud[i] < m.or){
      or.loss[i] = 0
    } else{
      or.loss[i] = min((xbar.vec-mud)^2,(mud-m.or)^2)
    }
  }
  avg.loss.mu.est.true.oracle = mean(or.loss)
  
  #*****SURE.G method*****************
  
  mu.est.SURE.G = thetahat.SURE.G.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.G = mean((mud-mu.est.SURE.G)^2)
  
  d.biv.SURE.G.kde.est = c(kde2d(c(mu.est.SURE.G,mu.est.SURE.G),c(S2.vec,-S2.vec),n=c(n1,n2),
                                 lims=c(range(mu.grid),range(sigma2.grid)))$z)
  d.mu.SURE.G.kde.est = density(mu.est.SURE.G, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.SURE.G= sum((d.biv.SURE.G.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.G= sum((d.mu.SURE.G.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  #************************************
  
  #*****SURE.M method*****************
  
  mu.est.SURE.M = thetahat.SURE.M.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.M = mean((mud-mu.est.SURE.M)^2)
  
  d.biv.SURE.M.kde.est = c(kde2d(c(mu.est.SURE.M,mu.est.SURE.M),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                    range(sigma2.grid)))$z)
  d.mu.SURE.M.kde.est = density(mu.est.SURE.M, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.SURE.M= sum((d.biv.SURE.M.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.M= sum((d.mu.SURE.M.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  
  #************************************
  
  #*****SURE.SG  method*****************
  
  mu.est.SURE.SG = thetahat.SURE.SG.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.SG = mean((mud-mu.est.SURE.SG)^2)
  
  d.biv.SURE.SG.kde.est = c(kde2d(c(mu.est.SURE.SG,mu.est.SURE.SG),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                       range(sigma2.grid)))$z)
  d.mu.SURE.SG.kde.est = density(mu.est.SURE.SG, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  MISE.biv.kde.SURE.SG= sum((d.biv.SURE.SG.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.SG= sum((d.mu.SURE.SG.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #************************************
  
  #*****SURE.SM  method*****************
  
  mu.est.SURE.SM = thetahat.SURE.SM.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.SM = mean((mud-mu.est.SURE.SM)^2)
  
  d.biv.SURE.SM.kde.est = c(kde2d(c(mu.est.SURE.SM,mu.est.SURE.SM),c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                       range(sigma2.grid)))$z)
  d.mu.SURE.SM.kde.est = density(mu.est.SURE.SM, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.SURE.SM= sum((d.biv.SURE.SM.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.SM= sum((d.mu.SURE.SM.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  #************************************
  
  #************************************************************  
  
  # group_linear: num bins = n^1/3
  mu.est.gl = grouplinear(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.gl = mean((mud-mu.est.gl)^2)
  
  d.biv.gl.kde.est = c(kde2d(c(mu.est.gl,mu.est.gl),c(S2.vec,-S2.vec),
                             n=c(n1,n2),lims=c(range(mu.grid),
                                               range(sigma2.grid)))$z)
  d.mu.gl.kde.est = density(mu.est.gl, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl= sum((d.biv.gl.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl= sum((d.mu.gl.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #************************************************************
  
  # group_linear: sure(equal-bins)
  mu.est.gl.SURE = grouplinear.sure(xbar.vec, S2.vec.xbar, kmax=min(ceiling(q^(1/3)/0.8),q))
  
  avg.loss.mu.est.gl.SURE = mean((mud-mu.est.gl.SURE)^2)
  
  d.biv.gl.SURE.kde.est = c(kde2d(c(mu.est.gl.SURE,mu.est.gl.SURE),c(S2.vec,-S2.vec),
                                  n=c(n1,n2),lims=c(range(mu.grid),
                                                    range(sigma2.grid)))$z)
  d.mu.gl.SURE.kde.est = density(mu.est.gl.SURE, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl.SURE= sum((d.biv.gl.SURE.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl.SURE= sum((d.mu.gl.SURE.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #************************************************************
  
  # dynamic_group_linear_all_division
  mu.est.gl.dynamic = GroupSure(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.gl.dynamic = mean((mud-mu.est.gl.dynamic)^2)
  
  d.biv.gl.dynamic.kde.est = c(kde2d(c(mu.est.gl.dynamic,mu.est.gl.dynamic),
                                     c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                         range(sigma2.grid)))$z)
  d.mu.gl.dynamic.kde.est = density(mu.est.gl.dynamic, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl.dynamic = sum((d.biv.gl.dynamic.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl.dynamic = sum((d.mu.gl.dynamic.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  mu.est.gl.dynamicMin = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*0.8)
  
  avg.loss.mu.est.gl.dynamicMin = mean((mud-mu.est.gl.dynamicMin)^2)
  
  d.biv.gl.dynamicMin.kde.est = c(kde2d(c(mu.est.gl.dynamicMin,mu.est.gl.dynamicMin),
                                        c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                            range(sigma2.grid)))$z)
  d.mu.gl.dynamicMin.kde.est = density(mu.est.gl.dynamicMin, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl.dynamicMin = sum((d.biv.gl.dynamicMin.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl.dynamicMin = sum((d.mu.gl.dynamicMin.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  
  mu.est.gl.dynamicMin2 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3))
  
  avg.loss.mu.est.gl.dynamicMin2 = mean((mud-mu.est.gl.dynamicMin2)^2)
  
  d.biv.gl.dynamicMin2.kde.est = c(kde2d(c(mu.est.gl.dynamicMin2,mu.est.gl.dynamicMin2),
                                         c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                             range(sigma2.grid)))$z)
  d.mu.gl.dynamicMin2.kde.est = density(mu.est.gl.dynamicMin2, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl.dynamicMin2 = sum((d.biv.gl.dynamicMin2.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl.dynamicMin2 = sum((d.mu.gl.dynamicMin2.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  mu.est.gl.dynamicMin3 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*1.2)
  
  avg.loss.mu.est.gl.dynamicMin3 = mean((mud-mu.est.gl.dynamicMin3)^2)
  
  d.biv.gl.dynamicMin3.kde.est = c(kde2d(c(mu.est.gl.dynamicMin3,mu.est.gl.dynamicMin3),
                                         c(S2.vec,-S2.vec),n=c(n1,n2),lims=c(range(mu.grid),
                                                                             range(sigma2.grid)))$z)
  d.mu.gl.dynamicMin3.kde.est = density(mu.est.gl.dynamicMin3, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  
  MISE.biv.kde.gl.dynamicMin3 = sum((d.biv.gl.dynamicMin3.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.gl.dynamicMin3 = sum((d.mu.gl.dynamicMin3.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  #************************************************************
  
  
  #**********SURE.M Double**************************
  
  mu.est.SURE.M.Double.object = estimates.SURE.M.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.M.Double = mu.est.SURE.M.Double.object$mu.est
  
  sigma2d.est.SURE.M.Double = mu.est.SURE.M.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.M.Double = mean((mud-mu.est.SURE.M.Double)^2)
  avg.loss.sigma2d.est.SURE.M.Double = mean((sigma2d-sigma2d.est.SURE.M.Double)^2)
  
  d.biv.SURE.M.Double.kde.est = c(kde2d(c(mu.est.SURE.M.Double,mu.est.SURE.M.Double),
                                        c(sigma2d.est.SURE.M.Double,-sigma2d.est.SURE.M.Double),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                                  range(sigma2.grid)))$z)
  d.mu.SURE.M.Double.kde.est = density(mu.est.SURE.M.Double, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.SURE.M.Double.kde.est = density(c(sigma2d.est.SURE.M.Double,-sigma2d.est.SURE.M.Double),
                                           n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.SURE.M.Double= sum((d.biv.SURE.M.Double.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.M.Double= sum((d.mu.SURE.M.Double.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.SURE.M.Double= sum((d.sigma2.SURE.M.Double.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  #*************************************************************
  
  #**********SURE.G Double**************************************
  
  
  mu.est.SURE.G.Double.object = estimates.SURE.G.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.G.Double = mu.est.SURE.G.Double.object$mu.est
  
  sigma2d.est.SURE.G.Double = mu.est.SURE.G.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.G.Double = mean((mud-mu.est.SURE.G.Double)^2)
  avg.loss.sigma2d.est.SURE.G.Double = mean((sigma2d-sigma2d.est.SURE.G.Double)^2)
  
  d.biv.SURE.G.Double.kde.est = c(kde2d(c(mu.est.SURE.G.Double,mu.est.SURE.G.Double),
                                        c(sigma2d.est.SURE.G.Double,-sigma2d.est.SURE.G.Double),
                                        n=c(n1,n2),lims=c(range(mu.grid),range(sigma2.grid)))$z)
  
  d.mu.SURE.G.Double.kde.est = density(mu.est.SURE.G.Double, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.SURE.G.Double.kde.est = density(c(sigma2d.est.SURE.G.Double,-sigma2d.est.SURE.G.Double),
                                           n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  
  MISE.biv.kde.SURE.G.Double= sum((d.biv.SURE.G.Double.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.G.Double= sum((d.mu.SURE.G.Double.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.SURE.G.Double= sum((d.sigma2.SURE.G.Double.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  #*************************************************************
  
  #**********SURE.SM Double**************************************
  
  
  mu.est.SURE.SM.Double.object = estimates.SURE.SM.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$mu.est
  
  sigma2d.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.SM.Double = mean((mud-mu.est.SURE.SM.Double)^2)
  avg.loss.sigma2d.est.SURE.SM.Double = mean((sigma2d-sigma2d.est.SURE.SM.Double)^2)
  
  d.biv.SURE.SM.Double.kde.est = c(kde2d(mu.est.SURE.SM.Double,sigma2d.est.SURE.SM.Double,n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                            range(sigma2.grid)))$z)
  d.mu.SURE.SM.Double.kde.est = density(mu.est.SURE.SM.Double, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.SURE.SM.Double.kde.est = density(sigma2d.est.SURE.SM.Double, n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  
  MISE.biv.kde.SURE.SM.Double= sum((d.biv.SURE.SM.Double.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.SM.Double= sum((d.mu.SURE.SM.Double.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.SURE.SM.Double= sum((d.sigma2.SURE.SM.Double.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  #*************************************************************
  
  #**********SURE.SG Double**************************************
  
  mu.est.SURE.SG.Double.object = estimates.SURE.SG.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$mu.est
  
  sigma2d.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.SG.Double = mean((mud-mu.est.SURE.SG.Double)^2)
  avg.loss.sigma2d.est.SURE.SG.Double = mean((sigma2d-sigma2d.est.SURE.SG.Double)^2)
  
  d.biv.SURE.SG.Double.kde.est = c(kde2d(c(mu.est.SURE.SG.Double,mu.est.SURE.SG.Double),
                                         c(sigma2d.est.SURE.SG.Double,-sigma2d.est.SURE.SG.Double),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                                     range(sigma2.grid)))$z)
  d.mu.SURE.SG.Double.kde.est = density(mu.est.SURE.SG.Double, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.SURE.SG.Double.kde.est = density(c(sigma2d.est.SURE.SG.Double,-sigma2d.est.SURE.SG.Double),
                                            n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.SURE.SG.Double= sum((d.biv.SURE.SG.Double.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.SURE.SG.Double= sum((d.mu.SURE.SG.Double.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.SURE.SG.Double= sum((d.sigma2.SURE.SG.Double.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  #*************************************************************
  
  #*********NIG mixture method**********************************
  
  MCMCoutput = DPMM.nig(x.vec, ids, gamma.pi=0.1, 
                        k=10, Burnin=5000, Simsize=10000, mu.grid=mu.grid, sigma2.grid=sigma2.grid, 
                        hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM = mean((MCMCoutput$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM = mean((MCMCoutput$sigma2.est-sigma2d)^2)
  
  d.biv.DPMM.est = MCMCoutput$d.biv.grid.est
  d.mu.DPMM.est = MCMCoutput$d.mu.grid.est
  d.sigma2.DPMM.est = MCMCoutput$d.sigma2.grid.est
  
  MISE.biv.DPMM= sum((d.biv.DPMM.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.DPMM= sum((d.mu.DPMM.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.DPMM= sum((d.sigma2.DPMM.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  d.biv.DPMM.kde.est = c(kde2d(c(MCMCoutput$mu.est,MCMCoutput$mu.est),
                               c(MCMCoutput$sigma2.est,-MCMCoutput$sigma2.est),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                 range(sigma2.grid)))$z)
  d.mu.DPMM.kde.est = density(MCMCoutput$mu.est, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.DPMM.kde.est = density(c(MCMCoutput$sigma2.est,-MCMCoutput$sigma2.est), n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.DPMM = sum((d.biv.DPMM.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.DPMM = sum((d.mu.DPMM.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.DPMM= sum((d.sigma2.DPMM.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  #*********NG mixture method**********************************
  
  prec.grid = 1/sigma2.grid
  
  MCMCoutput_ng = DPMM.ng(x.vec, ids, gamma.pi=0.1, 
                          k=10, Burnin=5000, Simsize=10000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM_ng = mean((MCMCoutput_ng$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng = mean((1/MCMCoutput_ng$prec.est-sigma2d)^2)
  
  d.biv.DPMM_ng.est = MCMCoutput_ng$d.biv.grid.est*rep(prec.grid^2,each=n1)
  d.mu.DPMM_ng.est = MCMCoutput_ng$d.mu.grid.est
  d.sigma2.DPMM_ng.est = MCMCoutput_ng$d.prec.grid.est*prec.grid^2
  
  MISE.biv.DPMM_ng = sum((d.biv.DPMM_ng.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.DPMM_ng = sum((d.mu.DPMM_ng.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.DPMM_ng= sum((d.sigma2.DPMM_ng.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  d.biv.DPMM_ng.kde.est = c(kde2d(c(MCMCoutput_ng$mu.est,MCMCoutput_ng$mu.est),
                                  c(1/MCMCoutput_ng$prec.est,-1/MCMCoutput_ng$prec.est),n=c(n1,n2),lims=c(range(mu.grid),
                                                                                                          range(sigma2.grid)))$z)
  d.mu.DPMM_ng.kde.est = density(MCMCoutput_ng$mu.est, n=n1, from=mu.grid[1], to=mu.grid[n1])$y
  
  d.sigma2.DPMM_ng.kde.est = density(c(1/MCMCoutput_ng$prec.est,-1/MCMCoutput_ng$prec.est), n=n2, from=sigma2.grid[1], to=sigma2.grid[n2])$y
  
  MISE.biv.kde.DPMM_ng = sum((d.biv.DPMM_ng.kde.est-d.biv.true)^2)*
    (mu.grid[2]-mu.grid[1])*(sigma2.grid[2]-sigma2.grid[1])
  
  MISE.mu.kde.DPMM_ng = sum((d.mu.DPMM_ng.kde.est-d.mu.true)^2)*
    (mu.grid[2]-mu.grid[1])
  
  MISE.sigma2.kde.DPMM_ng= sum((d.sigma2.DPMM_ng.kde.est-d.sigma2.true)^2)*
    (sigma2.grid[2]-sigma2.grid[1])
  
  
  #**************abhra density*********************************
  
  MCMCoutput_abhra = UNIVARIATE_DECON_HETEROSCEDASTIC(ws=x.vec, mis=ni.vec, error.fit=c("normal"), x.grid=mu.grid,
                                                      simsize=10000, burnin=5000, plot_results=FALSE)
  
  MISE.mu.abhra = sum((MCMCoutput_abhra$density.x.est-d.mu.true)^2)* (mu.grid[2]-mu.grid[1])
  
  #************************************************************  
  
  
  if(b%%50==0){
    print(b)
  }
  
  #************************************************************
  
  
  #*************************************************************
  
  obj1 = 
    c(q,
      avg.loss.mu.est.true.oracle,
      avg.loss.mu.est.sample.mean,
      avg.loss.mu.est.EBMLE,
      avg.loss.mu.est.EBMOM,
      avg.loss.mu.est.JS,
      avg.loss.mu.est.oracle,
      avg.loss.mu.est.SURE.G,
      avg.loss.mu.est.SURE.M,
      avg.loss.mu.est.SURE.SG,
      avg.loss.mu.est.SURE.SM,
      avg.loss.mu.est.gl,
      avg.loss.mu.est.gl.SURE,
      avg.loss.mu.est.gl.dynamic,
      avg.loss.mu.est.gl.dynamicMin,
      avg.loss.mu.est.gl.dynamicMin2,
      avg.loss.mu.est.gl.dynamicMin3,
      avg.loss.mu.est.SURE.M.Double,
      avg.loss.mu.est.SURE.G.Double,
      avg.loss.mu.est.SURE.SG.Double,
      avg.loss.mu.est.SURE.SM.Double,
      avg.loss.mu.est.DPMM,
      avg.loss.mu.est.DPMM_ng,
      
      avg.loss.sigma2d.est.sample.var,
      avg.loss.sigma2d.est.SURE.M.Double,
      avg.loss.sigma2d.est.SURE.G.Double,
      avg.loss.sigma2d.est.SURE.SG.Double,
      avg.loss.sigma2d.est.SURE.SM.Double,
      avg.loss.sigma2d.est.DPMM,
      avg.loss.sigma2d.est.DPMM_ng,
      
      MISE.biv.kde.true,
      MISE.biv.kde.sample.est,
      MISE.biv.kde.EBMLE,
      MISE.biv.kde.EBMOM,
      MISE.biv.kde.JS,
      MISE.biv.kde.oracle,
      MISE.biv.kde.SURE.G,
      MISE.biv.kde.SURE.M,
      MISE.biv.kde.SURE.SG,
      MISE.biv.kde.SURE.SM,
      MISE.biv.kde.gl,
      MISE.biv.kde.gl.SURE,
      MISE.biv.kde.gl.dynamic,
      MISE.biv.kde.gl.dynamicMin,
      MISE.biv.kde.gl.dynamicMin2,
      MISE.biv.kde.gl.dynamicMin3,
      MISE.biv.kde.SURE.G.Double,
      MISE.biv.kde.SURE.M.Double,
      MISE.biv.kde.SURE.SG.Double,
      MISE.biv.kde.SURE.SM.Double,
      MISE.biv.kde.DPMM,
      MISE.biv.DPMM,
      MISE.biv.kde.DPMM_ng,
      MISE.biv.DPMM_ng,
      
      
      MISE.mu.kde.true,
      MISE.mu.kde.sample.est,
      MISE.mu.kde.EBMLE,
      MISE.mu.kde.EBMOM,
      MISE.mu.kde.JS,
      MISE.mu.kde.oracle,
      MISE.mu.kde.SURE.G,
      MISE.mu.kde.SURE.M,
      MISE.mu.kde.SURE.SG,
      MISE.mu.kde.SURE.SM,
      MISE.mu.kde.gl,
      MISE.mu.kde.gl.SURE,
      MISE.mu.kde.gl.dynamic,
      MISE.mu.kde.gl.dynamicMin,
      MISE.mu.kde.gl.dynamicMin2,
      MISE.mu.kde.gl.dynamicMin3,
      MISE.mu.kde.SURE.G.Double,
      MISE.mu.kde.SURE.M.Double,
      MISE.mu.kde.SURE.SG.Double,
      MISE.mu.kde.SURE.SM.Double,
      MISE.mu.abhra,
      MISE.mu.kde.DPMM,
      MISE.mu.DPMM,
      MISE.mu.kde.DPMM_ng,
      MISE.mu.DPMM_ng,
      
      MISE.sigma2.kde.true,
      MISE.sigma2.kde.sample.est,
      MISE.sigma2.kde.SURE.G.Double,
      MISE.sigma2.kde.SURE.M.Double,
      MISE.sigma2.kde.SURE.SG.Double,
      MISE.sigma2.kde.SURE.SM.Double,
      MISE.sigma2.kde.DPMM,
      MISE.sigma2.DPMM,
      MISE.sigma2.kde.DPMM_ng,
      MISE.sigma2.DPMM_ng
      
    )
  
  
  obj1
  
}

distinctq = unique(as.vector(alloutput_diffq[1,]))
paperoutput = array(0,c(dim(alloutput_diffq)[1],length(distinctq)))

for(kk in 1:length(distinctq)){
  paperoutput[,kk] = apply(alloutput_diffq[,which(as.vector(alloutput_diffq[1,])==distinctq[kk])],1,mean)
}

finalans = data.frame(c("q",
                        "avg.loss.mu.est.true.oracle",
                        "avg.loss.mu.est.sample.mean",
                        "avg.loss.mu.est.EBMLE",
                        "avg.loss.mu.est.EBMOM",
                        "avg.loss.mu.est.JS",
                        "avg.loss.mu.est.oracle",
                        "avg.loss.mu.est.SURE.G",
                        "avg.loss.mu.est.SURE.M",
                        "avg.loss.mu.est.SURE.SG",
                        "avg.loss.mu.est.SURE.SM",
                        "avg.loss.mu.est.gl",
                        "avg.loss.mu.est.gl.SURE",
                        "avg.loss.mu.est.gl.dynamic",
                        "avg.loss.mu.est.gl.dynamicMin",
                        "avg.loss.mu.est.gl.dynamicMin2",
                        "avg.loss.mu.est.gl.dynamicMin3",
                        "avg.loss.mu.est.SURE.M.Double",
                        "avg.loss.mu.est.SURE.G.Double",
                        "avg.loss.mu.est.SURE.SG.Double",
                        "avg.loss.mu.est.SURE.SM.Double",
                        "avg.loss.mu.est.DPMM",
                        "avg.loss.mu.est.DPMM_ng",
                        
                        "avg.loss.sigma2d.est.sample.var",
                        "avg.loss.sigma2d.est.SURE.M.Double",
                        "avg.loss.sigma2d.est.SURE.G.Double",
                        "avg.loss.sigma2d.est.SURE.SG.Double",
                        "avg.loss.sigma2d.est.SURE.SM.Double",
                        "avg.loss.sigma2d.est.DPMM",
                        "avg.loss.sigma2d.est.DPMM_ng",
                        
                        "MISE.biv.kde.true",
                        "MISE.biv.kde.sample.est",
                        "MISE.biv.kde.EBMLE",
                        "MISE.biv.kde.EBMOM",
                        "MISE.biv.kde.JS",
                        "MISE.biv.kde.oracle",
                        "MISE.biv.kde.SURE.G",
                        "MISE.biv.kde.SURE.M",
                        "MISE.biv.kde.SURE.SG",
                        "MISE.biv.kde.SURE.SM",
                        "MISE.biv.kde.gl",
                        "MISE.biv.kde.gl.SURE",
                        "MISE.biv.kde.gl.dynamic",
                        "MISE.biv.kde.gl.dynamicMin",
                        "MISE.biv.kde.gl.dynamicMin2",
                        "MISE.biv.kde.gl.dynamicMin3",
                        "MISE.biv.kde.SURE.G.Double",
                        "MISE.biv.kde.SURE.M.Double",
                        "MISE.biv.kde.SURE.SG.Double",
                        "MISE.biv.kde.SURE.SM.Double",
                        "MISE.biv.kde.DPMM",
                        "MISE.biv.DPMM",
                        "MISE.biv.kde.DPMM_ng",
                        "MISE.biv.DPMM_ng",
                        
                        "MISE.mu.kde.true",
                        "MISE.mu.kde.sample.est",
                        "MISE.mu.kde.EBMLE",
                        "MISE.mu.kde.EBMOM",
                        "MISE.mu.kde.JS",
                        "MISE.mu.kde.oracle",
                        "MISE.mu.kde.SURE.G",
                        "MISE.mu.kde.SURE.M",
                        "MISE.mu.kde.SURE.SG",
                        "MISE.mu.kde.SURE.SM",
                        "MISE.mu.kde.gl",
                        "MISE.mu.kde.gl.SURE",
                        "MISE.mu.kde.gl.dynamic",
                        "MISE.mu.kde.gl.dynamicMin",
                        "MISE.mu.kde.gl.dynamicMin2",
                        "MISE.mu.kde.gl.dynamicMin3",
                        "MISE.mu.kde.SURE.G.Double",
                        "MISE.mu.kde.SURE.M.Double",
                        "MISE.mu.kde.SURE.SG.Double",
                        "MISE.mu.kde.SURE.SM.Double",
                        "MISE.mu.abhra",
                        "MISE.mu.kde.DPMM",
                        "MISE.mu.DPMM",
                        "MISE.mu.kde.DPMM_ng",
                        "MISE.mu.DPMM_ng",
                        
                        "MISE.sigma2.kde.true",
                        "MISE.sigma2.kde.sample.est",
                        "MISE.sigma2.kde.SURE.G.Double",
                        "MISE.sigma2.kde.SURE.M.Double",
                        "MISE.sigma2.kde.SURE.SG.Double",
                        "MISE.sigma2.kde.SURE.SM.Double",
                        "MISE.sigma2.kde.DPMM",
                        "MISE.sigma2.DPMM",
                        "MISE.sigma2.kde.DPMM_ng",
                        "MISE.sigma2.DPMM_ng"),
                      paperoutput)

names(finalans) = c("measure", rep("value",((q.upper.lim-q.lower.lim)/q.interval+1)))

print(finalans)

write.csv(finalans, file=outputfilecsv)

save(list = ls(), file = outputfile)


  #************************************************************************************
  # R program for estimating mean and variance of heteroscedastic normal error
  # Paper: SURE estimates for a heteroscedastic hierarchical model(Xie, Xianchao and Kou, SC and Brown, Lawrence D)
  # Paper: On SURE-Type Double Shrinkage Estimation(Jing, Bing-Yi and Li, Zhouping and Pan, Guangming and Zhou, Wang)
  #************************************************************************************  

  #************************************************************************************
  # We did not the find the codes online so we wrote the code after reading the paper. Later, we found a piece of code
  # at https://github.com/MaZhuang/grouplinear/functions_XKB.R, from there we took the functions thetahat.SURE.M.XKB,
  # thetahat.SURE.G.XKB, thetahat.SURE.SG.XKB. Rest all functions was written by us. 
  # Author of the code: Shyamalendu Sinha
  #************************************************************************************


# if (!require("rootSolve")) install.packages("rootSolve")
# if (!require("quadprog")) install.packages("quadprog")
# if (!require("isotone")) install.packages("isotone")
 library(rootSolve)
 library(quadprog)
 library(isotone)

#********EBMLE
thetahat.EBMLE.XKB = function(X,A){
  lambda0 = max(mean(X^2-A)-mean(X)^2,0.001)
  m0 = mean(X)
  param0 = c(lambda0,m0)
  
  marginal.log.likelihood = function(param){
    lambda = param[1]
    m = param[2]
    ans = sum(0.5*log(lambda+A)+(X-m)^2/(2*(lambda+A)))
    return(ans)
  }
  
  grad.marginal.log.likelihood = function(param){
    lambda = param[1]
    m = param[2]
    ans = c(sum(0.5/(lambda+A)-(X-m)^2/(2*(lambda+A)^2)),sum((m-X)/(2*(lambda+A))))
    return(ans)
  }
  
  param.EBMLE = constrOptim(param0, marginal.log.likelihood, 
                            grad.marginal.log.likelihood, ui=t(c(1,0)), ci=as.vector(0))$par
  mu.est.EBMLE = (param.EBMLE[1]*X+A*param.EBMLE[2])/(param.EBMLE[1]+A)
  return(as.vector(mu.est.EBMLE))
  
}


#***********************

#*******EBMOM**********

thetahat.EBMOM.XKB = function(X,A){
  q = length(X)
  lambda0 = max(mean(X^2-A)-mean(X)^2,0.001)
  m0 = mean(X)
  param0 = c(lambda0,m0)
  
  MOM.equation = function(param){
    lambda = param[1]
    m = param[2]
    ans1 = m*sum(1/(A+lambda))-sum(X/(A+lambda))
    ans2 = (q-1)*lambda-sum((X-m)^2)+(q-1)/q*sum(A)
    ans3 = c(ans1,ans2)
    return(ans3)
  }
  
  param.EBMOM = multiroot(MOM.equation,start=param0)$root
  param.EBMOM[1] = ifelse(param.EBMOM[1]>0,param.EBMOM[1],0)
  mu.est.EBMOM = (param.EBMOM[1]*X+A*param.EBMOM[2])/(param.EBMOM[1]+A)
  return(as.vector(mu.est.EBMOM))
}


#*****************************

#*************JS estimate*****

thetahat.JS.XKB = function(X,A){
  q = length(X)
  m.JS = sum(X/A)/sum(1/A)
  mu.est.JS = m.JS+ifelse((1-(q-3)/sum((X-m.JS)^2/A))>0,
                          (1-(q-3)/sum((X-m.JS)^2/A)),0)*(X-m.JS)
  return(mu.est.JS)
}


#****************************

#********Oracle**************

thetahat.oracle.XKB = function(X,A,mud){
  oracle.loss = function(param){
    lambda = param[1]
    m = param[2]
    ans = mean((lambda^2*(mud^2+A)+A^2*m^2+2*lambda*m*A*mud)/(lambda+A)^2+
                 -2*mud*(lambda*mud+A*m)/(lambda+A)+mud^2)
    return(ans)
  }
  
  lambda0 = max(mean(X^2-A)-mean(X)^2,0.001)
  m0 = mean(X)
  param0 = c(lambda0,m0)
  param.oracle = nlm(oracle.loss,param0)$estimate
  mu.est.oracle = (param.oracle[1]*X+A*param.oracle[2])/(param.oracle[1]+A)
  return(as.vector(mu.est.oracle))
}


#************************************
# 
# #*****SURE.G method*****************
# 
# thetahat.SURE.G.XKB = function(X,A){
#   q = length(X)
#   SURE.G.loss = function(lambda){
#     mean(A^2*(A+lambda)^(-2)*(X-mean(X))^2)+
#       mean(A*(A+lambda)^(-1)*(lambda-A+2*q^(-1)*A))
#   }
#   
#   lambda0 = max(mean(X^2-A)-mean(X)^2,0.001)
#   lambda.SURE.G = nlm(SURE.G.loss,lambda0)$estimate
#   mu.est.SURE.G = (lambda.SURE.G*X+A*mean(X))/(lambda.SURE.G+A)
#   return(as.vector(mu.est.SURE.G))
# }
# 
# 
# #************************************
# 
# #*****SURE.M method*****************
# 
# thetahat.SURE.M.XKB = function(X,A){
#   
#   SURE.M.loss = function(param){
#     lambda = param[1]
#     m = param[2]
#     ans = mean(A*(A+lambda)^(-2)*(A*(X-m)^2+lambda^2-A^2))
#     return(ans)
#   }
#   
#   lambda0 = max(mean(X^2-A)-mean(X)^2,0.001)
#   m0 = mean(X)
#   param0 = c(lambda0,m0)
#   param.SURE.M = nlm(SURE.M.loss,param0)$estimate
#   mu.est.SURE.M = (param.SURE.M[1]*X+A*param.SURE.M[2])/(param.SURE.M[1]+A)
#   return(as.vector(mu.est.SURE.M))
# }
# 
# 
# #************************************
# 
# #*****SURE.SG  method*****************
# 
# thetahat.SURE.SG.XKB = function(X,A){
#   X.sort = X[order(A)]
#   A.sort = A[order(A)]
#   
#   SURE.SG.loss = function(b.vec){
#     ans = mean(b.vec^2*(X.sort-mean(X.sort))^2+(1-2*(1-q^(-1)*b.vec))*A.sort)
#     return(ans)
#   }
#   
#   D.mat = diag(2*(X.sort-mean(X.sort))^2)
#   d.vec = 2*(1-q^(-1))*A.sort
#   A0.mat = diag(q)
#   A1.mat = -A0.mat
#   A2.mat = matrix(0,nrow = q-1, ncol = q)
#   for(i in 1:(q-1)){
#     A2.mat[i,i] = -1
#     A2.mat[i,i+1] = 1
#   }
#   A.mat = as.matrix(rbind(A0.mat,A1.mat,A2.mat))
#   b.vec.SURE.SG = solve.QP(Dmat = D.mat, dvec = d.vec, 
#                            Amat = t(A.mat),bvec = c(rep(0,q),rep(-1,q),rep(0,q-1)), meq = 0)$solution
#   mu.est.SURE.SG = (1-b.vec.SURE.SG)*X.sort+mean(X.sort)*b.vec.SURE.SG
#   return(as.vector(mu.est.SURE.SG[order(order(A))]))
# }



#************************************

#*****SURE.SM  method*****************

thetahat.SURE.SM.XKB = function(X,A){
  q = length(X)
  X.sort = X[order(A)]
  A.sort = A[order(A)]
  
  SURE.SM.loss = function(b.vec,m){
    ans = mean(b.vec^2*(X.sort-m)^2+(1-2*b.vec)*A.sort)
    return(ans)
  }
  
  m.old = mean(X)-1
  m = mean(X)
  count = 0
  while(abs(m.old-m)>0.00001){
    m.old = m
    D.mat = diag(2*(X.sort-m)^2)
    d.vec = 2*A.sort
    A0.mat = diag(q)
    A1.mat = -A0.mat
    A2.mat = matrix(0,nrow = q-1, ncol = q)
    for(i in 1:(q-1)){
      A2.mat[i,i] = -1
      A2.mat[i,i+1] = 1
    }
    A.mat = as.matrix(rbind(A0.mat,A1.mat,A2.mat))
    b.vec.SURE.SM = solve.QP(Dmat = D.mat, dvec = d.vec, 
                             Amat = t(A.mat),
                             bvec = c(rep(0,q),rep(-1,q),
                                      rep(0,q-1)), meq = 0)$solution
    m= sum(b.vec.SURE.SM^2*X.sort)/sum(b.vec.SURE.SM^2)
    #print(m)
    count = count + 1
  }
  mu.est.SURE.SM = (1-b.vec.SURE.SM)*X.sort+m*b.vec.SURE.SM
  return(as.vector(mu.est.SURE.SM[order(order(A))]))
}


#************************************


#**********SURE.M Double**************************

estimates.SURE.M.JLPZ = function(X,S2,n){
  
  SURE.M.Double.loss = function(param){
    m = param[1]
    lambda = param[2]
    nu = param[3]
    A = param[4]
    ans = mean((1+lambda)^(-2)*S2+lambda^2/(1+lambda)^2*
                 (m^2+X^2-S2-2*m*X))+
      mean((n+nu-3)^(-2)*2/(n+1)*((n-1)*S2)^2+(nu-2)^2/(n+nu-3)^2*
             ((nu*A/(nu-2))^2+((n-1)*S2)^2/(n^2-1)-2*nu*A*S2/(nu-2)))
    return(ans)
  }
  
  m0 = mean(X)
  lambda0 = mean(S2)/var(X)
  nu0 = 1
  A0 = 1/(mean(1/S2))
  param0 = c(m0,lambda0,nu0,A0)
  param.SURE.M.Double = nlm(SURE.M.Double.loss,param0)$estimate
  mu.est.SURE.M.Double = (X+param.SURE.M.Double[2]*param.SURE.M.Double[1])/
    (1+param.SURE.M.Double[2])
  sigma2d.est.SURE.M.Double = n*(((n-1)*S2+param.SURE.M.Double[3]*param.SURE.M.Double[4])/
                                 (n+param.SURE.M.Double[3]-3))
  return(list(mu.est=mu.est.SURE.M.Double, sigma2.est=sigma2d.est.SURE.M.Double))
}


#*************************************************************

#**********SURE.G Double**************************************

estimates.SURE.G.JLPZ = function(X,S2,n){
  q = length(X)
  SURE.G.Double.loss = function(param){
    lambda = param[1]
    nu = param[2]
    A = param[3]
    ans = mean((1+lambda)^(-2)*(1+2*lambda/q)*S2+lambda^2/(1+lambda)^2*
                 (mean(X)^2+X^2-S2-2*mean(X)*X+2/q*S2))+
      mean((n+nu-3)^(-2)*2/(n+1)*((n-1)*S2)^2+(nu-2)^2/(n+nu-3)^2*
             ((nu*A/(nu-2))^2+((n-1)*S2)^2/(n^2-1)-2*nu*A*S2/(nu-2)))
    return(ans)
  }
  
  lambda0 = mean(S2)/var(X)
  nu0 = 1
  A0 = 1/(mean(1/S2))
  param0 = c(lambda0,nu0,A0)
  param.SURE.G.Double = nlm(SURE.G.Double.loss,param0)$estimate
  mu.est.SURE.G.Double = (X+param.SURE.G.Double[1]*mean(X))/
    (1+param.SURE.G.Double[1])
  sigma2d.est.SURE.G.Double = n*(((n-1)*S2+param.SURE.G.Double[2]*param.SURE.G.Double[3])/
                                 (n+param.SURE.G.Double[2]-3))
  return(list(mu.est=mu.est.SURE.G.Double, sigma2.est=sigma2d.est.SURE.G.Double))
}

#*************************************************************

#**********SURE.SM Double**************************************

estimates.SURE.SM.JLPZ = function(X,S2,n){
  
  SURE.SM.Double.loss = function(param){
    b1 = param[1]
    b2 = param[2]
    m = param[3]
    sigma2.center = param[4]
    
    ans = mean((1-2*b1)*S2+b1^2*(X-m)^2)+
      mean((1-2*b2)*2/(n-1)*((n-1)*S2)^2/(n^2-1)+b2^2*(sigma2.center-S2)^2)
    return(ans)
  }
  
  grad.SURE.SM.Double.loss = function(param){
    b1 = param[1]
    b2 = param[2]
    m = param[3]
    sigma2.center = param[4]
    ans = c(2*b1*mean((X-m)^2)-2*mean(S2),
            2*b2*mean((sigma2.center-S2)^2)-2*2/(n-1)*mean(((n-1)*S2)^2)/(n^2-1),
            b1^2*2*mean(m-X),
            b2^2*2*mean(sigma2.center-S2)
    )
    return(ans)
  }
  
  b1.0 = 0.01
  b2.0 = 0.01
  m0 = mean(X)
  sigma2.center0 = mean(S2)
  param0 = c(b1.0,b2.0,m0,sigma2.center0)
  
  param.SURE.SM.Double = constrOptim(param0, SURE.SM.Double.loss, grad.SURE.SM.Double.loss,
                                     ui=rbind(c(-1,0,0,0),c(0,-1,0,0)), ci=c(-1,-1))$par
  
  
  mu.est.SURE.SM.Double = X*(1-param.SURE.SM.Double[1])+
    param.SURE.SM.Double[1]*param.SURE.SM.Double[3]
  
  sigma2d.est.SURE.SM.Double = n*(S2*(1-param.SURE.SM.Double[2])+
                                  param.SURE.SM.Double[2]*param.SURE.SM.Double[4])
  
  return(list(mu.est=mu.est.SURE.SM.Double, sigma2.est=sigma2d.est.SURE.SM.Double))
}

#*************************************************************

#**********SURE.SG Double**************************************

estimates.SURE.SG.JLPZ = function(X,S2,n){
  q = length(X)
  SURE.SG.Double.loss = function(param){
    b1 = param[1]
    b2 = param[2]
    sigma2.center = param[3]
    
    ans = mean((1-2*b1*(1-1/q))*S2+b1^2*(X-mean(X))^2)+
      mean((1-2*b2)*2/(n-1)*((n-1)*S2)^2/(n^2-1)+b2^2*(sigma2.center-S2)^2)
    return(ans)
  }
  
  grad.SURE.SG.Double.loss = function(param){
    b1 = param[1]
    b2 = param[2]
    sigma2.center = param[3]
    ans = c(2*b1*mean((X-mean(X))^2)-2*(1-1/q)*mean(S2),
            2*b2*mean((sigma2.center-S2)^2)-2*2/(n-1)*mean(((n-1)*S2)^2)/(n^2-1),
            b2^2*2*mean(sigma2.center-S2)
    )
    return(ans)
  }
  
  b1.0 = 0.01
  b2.0 = 0.01
  sigma2.center0 = mean(S2)
  param0 = c(b1.0,b2.0,sigma2.center0)
  
  param.SURE.SG.Double = constrOptim(param0, SURE.SG.Double.loss, grad.SURE.SG.Double.loss,
                                     ui=rbind(c(-1,0,0),c(0,-1,0)), ci=c(-1,-1))$par
  
  
  mu.est.SURE.SG.Double = X*(1-param.SURE.SG.Double[1])+
    param.SURE.SG.Double[1]*mean(X)
  
  sigma2d.est.SURE.SG.Double = n*(S2*(1-param.SURE.SG.Double[2])+
                                    param.SURE.SG.Double[2]*param.SURE.SG.Double[3])
  
  return(list(mu.est=mu.est.SURE.SG.Double, sigma2.est=sigma2d.est.SURE.SG.Double))
}


#*************************************************************

# Rest of the code was copied from https://github.com/MaZhuang/grouplinear/functions_XKB.R

thetahat.SURE.M.XKB = function(X,A){
  
  # SURE(lambda,mu)
  f = function(par,X,A){
    lambda = par[1]
    mu = par[2]
    sum(A/(lambda+A)^2 * (A * (X-mu)^2 + lambda^2 - A^2))
  }
  
  # d/dlambda{SURE(lambda,mu=mu.hat.SURE(lambda))} (proportional to)
  g = function(lambda,X,A){  
    sum(A^2/(lambda+A)^3 * (X-(sum(A^2/(lambda+A)^2 *X) / sum(A^2/(lambda+A)^2)))^2 - A^2/(lambda+A)^2)
    #equivalent to the following(which is just easier to read):
    #mu = sum(A^2/(lambda+A)^2*X) / sum( A^2/(lambda+A)^2)
    #sum(A^2/(lambda+A)^3 * (X-mu)^2 - A^2/(lambda+A)^2)
  }
  
  # bayes rule for fixed lambda,mu
  thetahat = function(X,A,lambda,mu){
    lambda/(lambda+A) * X + A/(lambda+A) * mu
  }
  
  lambda.sure = ifelse(g(0,X=X,A=A)*g(max(A)*1000,X=X,A=A) < 0, 
                        uniroot(g,c(0,max(A)*1000),X=X,A=A, tol=1e-9)$root, 
                        optim(c(mean(pmax((X-mean(X))^2 - A,0)), mean(X)),f,X=X,A=A, 
                              method = "L-BFGS-B",lower=c(0,-Inf))$par[1])
  mu.sure = sum(A^2/(lambda.sure+A)^2 * X ) / sum(A^2/(lambda.sure+A)^2)
  return(thetahat(X,A,lambda.sure,mu.sure))
}



thetahat.SURE.G.XKB = function(X,A){
  q = length(X)
  # SURE.G(lambda)
  f.G = function(lambda,X,A){  
    sum((A/(A+lambda))^2 * (X - mean(X))^2 + A/(A+lambda)*(lambda - A + 2/q*A))
  }
  
  # bayes rule for fixed lambda,mu
  thetahat = function(X,A,lambda,mu){
    lambda/(lambda+A) * X + A/(lambda+A) * mu
  } 
  
  lambda = optimize(f.G,lower=0,upper=1000,X=X,A=A)$minimum
  return(thetahat(X,A,lambda,mean(X)))
}

thetahat.SURE.SG.XKB = function(X,A){
  q = length(X)
  fit = gpava(z = A, y = A * (1-1/q) /(X-mean(X))^2, weights = (X-mean(X))^2, 
                solver = weighted.mean, ties="primary")
  bhat = pmin(pmax(fit$x,0),1)
  return((1-bhat) * X + bhat * mean(X))
}


#************************************************************************************
# R program for fitting mixture of Normal-invgamma or Normal-gamma in 
# with heteroscedastic measurement error with \sigma_i^2 known or unknown
#************************************************************************************  

#************************************************************************************
# Author(s)
#
# Author of the code: Shyamalendu Sinha
#************************************************************************************

# if (!require("MCMCpack")) install.packages("MCMCpack")
# if (!require("invgamma")) install.packages("invgamma")
# if (!require("truncnorm")) install.packages("truncnorm")
library(MCMCpack)
library(invgamma)
library(truncnorm)
rinvgamma = invgamma::rinvgamma
dinvgamma = invgamma::dinvgamma

#**********************************************************************************

# X_{ij} = \mu_i+sqrt(\sigma_i^2) e_{ij}. (\mu_i,\sigma_i^2) follows mixture of normal-inverse-gamma. 
# (\mu_i,\sigma_i^2) not observed. i-1,...,q, j=1,..,n_i, and n_i>=2 for all i.

#**********************************************************************************

# input: x.vec = all X_{ij} in one vector of length \sum_{i=1}^q n_i.
#        ids = vector of indicators which denote the the $i$ of the corresponding X_{ij}. Range is between 1,..,q. 
#             The number i each appears n_i times.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.nig = function(x.vec, ids, gamma.pi = 0.1, k = 10,
                    Burnin = 5000, Simsize = 10000, mu.grid = c('Default'), sigma2.grid = c('Default'),
                    hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("invgamma")) install.packages("invgamma")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(invgamma)
  # library(truncnorm)
  
  # rinvgamma = invgamma::rinvgamma
  # dinvgamma = invgamma::dinvgamma
  
  orderid = order(ids)
  x.vec = x.vec[orderid]
  ids = ids[orderid]
  
  #*********************************Data Management****************************************#
  #number of replictions of mu.i. nivec[i]=5 means 5 replicates for mu.i. For balanced data all elements nivec is same.
  #For unbalanced data we need each elements to greater than equal 2 for identifiability of the joint distribution.
  ni.vec = as.data.frame(table(ids))$Freq
  q = length(ni.vec)
  
  if(sum(as.vector(ni.vec)<2)>0){
    print("Not enough data to estimate mean variance simultaneously")
  } else{
    
    #centering and scaling the whole dataset
    MM = mean(x.vec)
    SS = sd(x.vec)
    # MM = 0
    # SS = 1
    x1.vec = (x.vec-MM)/SS
    
    #******************************End of Data Management************************************#
    
    #******************************Initialize MCMC parameters************************************#
    
    #intial value of mu.i vector which is muvec. Similarly, sigma2vec is the vector sigma2.i's of length q. 
    xbar.vec = mu.vec = tapply(x1.vec,ids,mean)
    S2.vec = sigma2.vec = tapply(x1.vec,ids,var)
    
    #k-means clusters to determine initial values of cluster parameters.
    mscluster = kmeans(cbind(xbar.vec*SS+MM,S2.vec*SS^2),k)
    
    #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
    #initial value of mvec
    m.vec = data.frame(mscluster$center)[,1]
    
    m.vec = (m.vec-MM)/SS
    
    #latent variable which indicates cluster number. zi=r means (mu.i,sigma.i^2|z.i=r) \sim NIG(m.r, lambda.r, alpha.r,
    # beta.r)
    z.vec = mscluster$cluster
    
    #vector of M-H ratios of length k.
    mh.ratio1.vec = mh.ratio2.vec = rep(0,k)
    
    #matrix to store mu.i, sigma2.i and Z.i after each iteration.
    mu.vec.tot = sigma2.vec.tot = rep(0,q)
    mu.vec2.tot = sigma2.vec2.tot = rep(0,q)
    zvec.mat = array(0,c(q,Simsize))
    
    #vector of length of unique z.i's in each iteration.
    k.vec = vector()
    
    #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters. 
    #Empty matrix to sum of each iteration.
    theta.all = array(0,c(k,5))
    theta2.all = array(0,c(k,5))
    
    #initial values alphavec, betavec, lambdavec.
    lambda.vec = rep(1,k)
    alpha_0 = (mean(S2.vec))^2/var(S2.vec)+2
    beta_0 = ((mean(S2.vec))^2/var(S2.vec)+1)*mean(S2.vec)
    alpha.vec = rep(alpha_0,k)
    beta.vec = rep(beta_0,k)
    pi.vec = rep(1/k,k)
    
    
    if(sum(hyperparameters==c('Default'))>0){
      #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). 
      #So, setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
      m.0 = mean(x1.vec)
      tau2 = var(xbar.vec)
      
      #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
      #lambda is generally less than 1.
      a.lambda = 1
      b.lambda = 1
      
      #alphavec follows Gamma(a.alpha,b.alpha).
      a.alpha = 1
      b.alpha = 1
      
      #betavec follows Gamma(a.beta,b.beta).
      a.beta = 1
      b.beta = 1/(sd(x.vec)^2)*SS^2
    }
    
    if (sum(hyperparameters!=c('Default'))>0){
      
      m.0 = (hyperparameters[1]-MM)/SS
      tau2 = hyperparameters[2]/SS^2
      a.lambda = hyperparameters[3]
      b.lambda = hyperparameters[4]
      a.alpha = hyperparameters[5]
      b.alpha = hyperparameters[6]
      a.beta = hyperparameters[7]
      b.beta = hyperparameters[8]*SS^2
    }
    
    
    #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
    #Dirichlet(gamma.pi/k,...,gamma.pi/k). 
    #Low value of gamma.pi is non-informative, allowing higher variation in pi.r's. gamma.pi = 0.1
    #[mu.grid,sigma2.grid] is a bivariate grid on which we are estimating the density.
    if (sum(mu.grid==c('Default'))>0){
      mu.grid = seq(min(xbar.vec)-IQR(xbar.vec), max(xbar.vec)+IQR(xbar.vec),
                    length.out=2*100+1)[seq(2,2*100,2)]
    }
    if (sum(sigma2.grid==c('Default'))>0){
      sigma2.grid = seq(max(min(S2.vec)-IQR(S2.vec),0.001), max(S2.vec)+IQR(S2.vec),
                        length.out=2*100+1)[seq(2,2*100,2)]
    }
    
    n1 = length(mu.grid)
    n2 = length(sigma2.grid)
    
    biv.grid = array(0,c(n1*n2,2))
    biv.grid[,1] = rep(mu.grid,n2)
    biv.grid[,2] = rep(sigma2.grid,each=n1)
    
    #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate, mu-marginal, 
    #and sigma2-marginal distribution.
    d.biv.grid.tot = rep(0,n1*n2)
    d.mu.grid.tot = rep(0,n1)
    d.sigma2.grid.tot = rep(0,n2)
    total.log.like.vec = data.log.like.vec = rep(0,Simsize)
    
    #function to calculate bivariate density on each grid points
    d.biv.mixture = function(mu.grid, sigma2.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dnorm(mu.grid, mean=m.vec[r], sd=sqrt(sigma2.grid/lambda.vec[r]))*
          dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    #function to calculate mu-marginal density on each grid points
    d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
          dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
      }
      return(ans)
    }
    
    #function to calculate sigma2-marginal density on each grid points
    d.sigma2.mixture = function(sigma2.grid, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    #function to calculate log-likelihood of involving alpha
    loglike.alpha = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.alpha, b.alpha, r){
      ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
      return(ans)
    }
    
    #function to calculate log-likelihood of involving beta
    loglike.beta = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.beta, b.beta, r){
      ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
      return(ans)
    }
    
    #parameters should be rescaled, it is on x1.vec
    total.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      mu.vec.large = rep(mu.vec,ni.vec)
      sigma2.vec.large = rep(sigma2.vec,ni.vec)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x1.vec, mean=mu.vec.large, sd=sqrt(sigma2.vec.large), log=TRUE))+
        sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
        sum(count.vec*log(pi.vec))+
        sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
        sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
        sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
        sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
        lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
      
      return(ans)
    }
    
    #parameters should be in original scale, it is on x.vec
    data.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      mu.vec.large = rep(mu.vec,ni.vec)
      sigma2.vec.large = rep(sigma2.vec,ni.vec)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x.vec, mean=mu.vec.large, sd=sqrt(sigma2.vec.large), log=TRUE))+
        sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
        sum(count.vec*log(pi.vec))
      return(ans)
    }
    
    #*******************************MCMC algorithm*******************************************#
    for (l in 1:(Burnin+Simsize)){
      
      #updating mu.is
      mu.vec = rnorm(q, mean=(ni.vec*xbar.vec+m.vec[z.vec]*lambda.vec[z.vec])/(ni.vec+lambda.vec[z.vec]), 
                     sd=sqrt(sigma2.vec/(ni.vec+lambda.vec[z.vec])))
      
      #updating sigma2.is
      mu.vec.large = rep(mu.vec,ni.vec)
      
      sigma2.vec.shape = 0.5*(ni.vec+1)+alpha.vec[z.vec]
      
      sigma2.vec.rate = 0.5*tapply((x1.vec-mu.vec.large)^2,ids,sum)+
        0.5*lambda.vec[z.vec]*(mu.vec-m.vec[z.vec])^2+beta.vec[z.vec]
      
      sigma2.vec = rinvgamma(q, shape=sigma2.vec.shape, rate=sigma2.vec.rate)
      
      #updating z.is
      lognig.mat = array(0,dim=c(q,k))
      for (r in 1:k){
        lognig.mat[,r] = dnorm(mu.vec,mean=m.vec[r],sd=sqrt(sigma2.vec/lambda.vec[r]),log=TRUE) + 
          dinvgamma(sigma2.vec,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
      }
      
      lognig1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + lognig.mat
      
      #need it for computational purpose. Otherwise exp(small log values) going to be 0.
      lognig2.mat = lognig1.mat - apply(lognig1.mat,1,max)
      
      nig.mat = exp(lognig2.mat)
      
      cum.nig.mat = nig.mat %*% upper.tri(diag(k), diag = TRUE)
      
      random.q.vec = runif(q)*rowSums(nig.mat)
      
      z.vec = rowSums(random.q.vec > cum.nig.mat)+1
      
      #countvec is a vector of length k where each element is number of z.is in each cluster
      count.vec = tabulate(z.vec,nbins=k)
      
      #updating pi.rs
      pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
      
      #updating m.rs
      for (r in 1:k){
        if(count.vec[r] > 0){
          m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mu.vec/sigma2.vec)[which(z.vec==r)])+m.0/tau2)/
                             (lambda.vec[r]*sum((1/sigma2.vec)[which(z.vec==r)])+1/tau2),
                           sd=sqrt(1/(lambda.vec[r]*sum((1/sigma2.vec)[which(z.vec==r)])+1/tau2)))
        }
        if(count.vec[r] == 0){
          m.vec[r] = rnorm(1, mean=m.0, sd=sqrt(tau2))
        }
      }
      
      #updating lambda.rs
      for (r in 1:k){
        if(count.vec[r] > 0){
          lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                                 rate=sum((mu.vec[which(z.vec==r)]-m.vec[r])^2/
                                            (2*sigma2.vec[which(z.vec==r)]))+b.lambda)
        }
        if(count.vec[r] == 0){
          lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
        }
      }
      
      #updating alpha.rs
      alphacand.vec = rnorm(k, mean=alpha.vec, sd=1)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = loglike.alpha(alphacand.vec,beta.vec,sigma2.vec,z.vec,a.alpha,b.alpha,r)-
            loglike.alpha(alpha.vec,beta.vec,sigma2.vec,z.vec,a.alpha,b.alpha,r)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
            dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
      }
      
      #updating beta.rs
      betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS^2)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = loglike.beta(alpha.vec,betacand.vec,sigma2.vec,z.vec,a.beta,b.beta,r)-
            loglike.beta(alpha.vec,beta.vec,sigma2.vec,z.vec,a.beta,b.beta,r)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
            dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
      }
      
      d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS+MM, 
                                  lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS+MM, lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.sigma2.iter1 = d.sigma2.mixture(sigma2.grid, alpha.vec, beta.vec*SS^2, pi.vec)
      
      theta = data.frame(cbind(m.vec*SS+MM,lambda.vec,alpha.vec,beta.vec*SS^2,pi.vec))
      names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
      
      # if(l%%1000==0){
      #   theta1 = theta
      #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
      #   theta1 = theta1[order(c(-theta1$pi.vec)),]
      #   print(l)
      #   print(theta1)
      # }
      
      if (l > Burnin){
        #1st stage parameters
        mu.vec.tot = mu.vec.tot+mu.vec
        mu.vec2.tot = mu.vec2.tot+mu.vec^2
        sigma2.vec.tot = sigma2.vec.tot+sigma2.vec
        sigma2.vec2.tot = sigma2.vec2.tot+sigma2.vec^2
        zvec.mat[,(l-Burnin)] = z.vec
        
        #nig parameters
        theta.all = theta.all+as.matrix(theta)
        theta2.all = theta2.all+as.matrix(theta^2)
        
        #estimated density
        d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
        d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
        d.sigma2.grid.tot = d.sigma2.grid.tot + d.sigma2.iter1
        
        #total.log.like is computed based on x1.vec
        total.log.like.vec[(l-Burnin)] = total.log.like(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, 
                                                        beta.vec, z.vec, pi.vec)
        #data.log.like is compared on x.vec
        data.log.like.vec[(l-Burnin)] = data.log.like(mu.vec*SS+MM, sigma2.vec*SS^2, m.vec*SS+MM, lambda.vec, 
                                                      alpha.vec, beta.vec*SS^2, z.vec, pi.vec)
        k.vec[(l-Burnin)] = length(unique(z.vec))
      }
      
    }
    
    
    #mu.i amd sigma^2.i's posterior point estimate and it's posterior variance estimate from MCMC iterations. 
    mu.est = mu.vec.tot/Simsize*SS+MM
    mu.est.var = (mu.vec2.tot/Simsize-(mu.vec.tot/Simsize)^2)*SS^2
    sigma2.est = sigma2.vec.tot/Simsize*SS^2
    sigma2.est.var = (sigma2.vec2.tot/Simsize-(sigma2.vec.tot/Simsize)^2)*SS^4
    
    #Summarizing density estimate from MCMC iterations. 
    d.biv.grid.est = d.biv.grid.tot/Simsize
    d.mu.grid.est = d.mu.grid.tot/Simsize
    d.sigma2.grid.est = d.sigma2.grid.tot/Simsize
    
    #Because of level switching these estimates are not reliable.
    theta.hat = data.frame(theta.all/Simsize)
    theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
    
    names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    
    order.pi.vec = order(c(-theta.hat$pi.vec))
    theta.hat = theta.hat[order.pi.vec,]
    theta.hat.var = theta.hat.var[order.pi.vec,]
    
    #Clustering analysis
    find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
    find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans2=max(tab)/length(x); rm(tab); return(ans2)}
    
    #Estimate of cluster of (mu.i,sigma.i^2).
    z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
    z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
    
    
    theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    order.pi.vec.eff = order(-theta.eff$pi.vec.est)
    theta.eff = theta.eff[order.pi.vec.eff,]
    theta.eff.var = theta.eff.var[order.pi.vec.eff,]
    
    return(list(mu.est=mu.est, mu.est.var=mu.est.var, sigma2.est=sigma2.est, 
                sigma2.est.var=sigma2.est.var, 
                z.vec.est=z.vec.est, z.vec.est.prop=z.vec.est.prop, theta.hat=theta.hat,
                theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
                d.biv.grid.est=d.biv.grid.est,
                d.mu.grid.est=d.mu.grid.est, d.sigma2.grid.est=d.sigma2.grid.est,
                total.log.like.vec=total.log.like.vec,
                data.log.like.vec=data.log.like.vec, k.vec=k.vec))
    
  }
}

# Output: mu.est = estimated value of mu_i's based on MCMC iterations.
#         mu.est.var = estimated posterior variance of mu_i.
#         sigma2.est = estimated value of sigma_i^2's based on MCMC iterations.
#         sigma2.est.var = estimated posterior variance of sigma_i^2.
#         z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         d.sigma2.grid.est = estimated sigma2 density on sigma2.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 


#*******************************************************************************************************

# X_{i} = \mu_i+sqrt(\sigma_i^2) e_{i}. (\mu_i,\sigma_i^2) follows mixture of normal-inverse-gamma. \mu_i not observed.
# But \sigma_i^2 are observed. i=1,...,q.

#*******************************************************************************************************

# input: x.vec = all X_{i} in one vector of length q.
#        sigma2d = vector of of legth q where ith entry is \sigma_i^2 which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.nig.sknown = function(x.vec, sigma2d, gamma.pi = 0.1, k = 10,
                           Burnin = 5000, Simsize = 10000, mu.grid = c('Default'), sigma2.grid = c('Default'),
                           hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("invgamma")) install.packages("invgamma")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(invgamma)
  # library(truncnorm)
  
  # rinvgamma = invgamma::rinvgamma
  # dinvgamma = invgamma::dinvgamma
  
  #**********************************************************************************
  
  q = length(x.vec)
  
  #**********************************************************************************
  if(length(x.vec) != length(sigma2d)){
    print("length of x.vec and sigma2d does not match")
  } else{
    
    #centering and scaling the whole dataset
    MM = mean(x.vec)
    SS = sd(x.vec)
    # MM = 0
    # SS = 1
    x1.vec = (x.vec-MM)/SS
    sigma2d1 = sigma2d/SS
    
    #******************************End of Data Management************************************#
    
    #******************************Initialize MCMC parameters************************************#
    
    #intial value of mu.i vector which is muvec. Similarly, sigma2vec is the vector sigma2.i's of length q. 
    mu.vec = x1.vec
    sigma2.vec = sigma2d1
    
    #k-means clusters to determine initial values of cluster parameters.
    mscluster = kmeans(cbind(x.vec, sigma2d),k)
    
    #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
    #initial value of mvec
    m.vec = data.frame(mscluster$center)[,1]
    
    m.vec = (m.vec-MM)/SS
    
    #latent variable which indicates cluster number. zi=r means (mu.i,sigma.i^2|z.i=r) 
    #\sim NIG(m.r, lambda.r, alpha.r, beta.r)
    z.vec = mscluster$cluster
    
    #vector of M-H ratios of length k.
    mh.ratio1.vec = mh.ratio2.vec = rep(0,k)
    
    #matrix to store mu.i, sigma2.i and Z.i after each iteration.
    mu.vec.tot = rep(0,q)
    mu.vec2.tot = rep(0,q)
    zvec.mat = array(0,c(q,Simsize))
    
    #vector of length of unique z.i's in each iteration.
    k.vec = vector()
    
    #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters. 
    #Empty matrix to sum of each iteration.
    theta.all = array(0,c(k,5))
    theta2.all = array(0,c(k,5))
    
    #initial values alphavec, betavec, lambdavec.
    lambda.vec = rep(1,k)
    alpha_0 = (mean(sigma2d1))^2/var(sigma2d1)+2
    beta_0 = ((mean(sigma2d1))^2/var(sigma2d1)+1)*mean(sigma2d1)
    alpha.vec = rep(alpha_0,k)
    beta.vec = rep(beta_0,k)
    pi.vec = rep(1/k,k)
    
    
    if(sum(hyperparameters==c('Default'))>0){
      #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). So, setting 
      #tau2 = Var(\bar{X_i}) less informative than empirical bayes.
      m.0 = mean(x1.vec)
      tau2 = var(x1.vec)
      
      #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
      #lambda is generally less than 1.
      a.lambda = 1
      b.lambda = 1
      
      #alphavec follows Gamma(a.alpha,b.alpha).
      a.alpha = 1
      b.alpha = var(1/sigma2d1)/mean(1/sigma2d1)^2
      
      #betavec follows Gamma(a.beta,b.beta).
      a.beta = 1
      b.beta = var(1/sigma2d1)/mean(1/sigma2d1)
    }
    
    if (sum(hyperparameters!=c('Default'))>0){
      
      m.0 = (hyperparameters[1]-MM)/SS
      tau2 = hyperparameters[2]/SS^2
      a.lambda = hyperparameters[3]
      b.lambda = hyperparameters[4]
      a.alpha = hyperparameters[5]
      b.alpha = hyperparameters[6]
      a.beta = hyperparameters[7]
      b.beta = hyperparameters[8]*SS^2
    }
    
    
    #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
    #Dirichlet(gamma.pi/k,...,gamma.pi/k). 
    #Low value of gamma.pi non-informative. allows higher variation in pi.r's. gamma.pi = 0.1
    #[mu.grid,sigma2.grid] is a bivariate grid on which we are estimating the density.
    if (sum(mu.grid==c('Default'))>0){
      mu.grid = seq(min(x1.vec)-IQR(x1.vec), max(x1.vec)+IQR(x1.vec),
                    length.out=2*100+1)[seq(2,2*100,2)]
    }
    if (sum(sigma2.grid==c('Default'))>0){
      sigma2.grid = seq(max(min(sigma2d1)-IQR(sigma2d1),0.001), max(sigma2d1)+IQR(sigma2d1),
                        length.out=2*100+1)[seq(2,2*100,2)]
    }
    
    n1 = length(mu.grid)
    n2 = length(sigma2.grid)
    
    biv.grid = array(0,c(n1*n2,2))
    biv.grid[,1] = rep(mu.grid,n2)
    biv.grid[,2] = rep(sigma2.grid,each=n1)
    
    #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate, mu-marginal, 
    #and sigma2-marginal distribution.
    d.biv.grid.tot = rep(0,n1*n2)
    d.mu.grid.tot = rep(0,n1)
    d.sigma2.grid.tot = rep(0,n2)
    total.log.like.vec = data.log.like.vec = rep(0,Simsize)
    
    d.biv.mixture = function(mu.grid, sigma2.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dnorm(mu.grid, mean=m.vec[r], sd=sqrt(sigma2.grid/lambda.vec[r]))*
          dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
          dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
      }
      return(ans)
    }
    
    d.sigma2.mixture = function(sigma2.grid, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    loglike.alpha = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.alpha, b.alpha, r){
      ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
      return(ans)
    }
    
    loglike.beta = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.beta, b.beta, r){
      ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
      return(ans)
    }
    
    # parameters should be in rescaled scale. As it is on x1.vec
    total.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x1.vec, mean=mu.vec, sd=sqrt(sigma2.vec), log=TRUE))+
        sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
        sum(count.vec*log(pi.vec))+
        sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
        sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
        sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
        sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
        lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
      
      return(ans)
    }
    
    
    # parameters should be in original scale. As it is on x.vec
    data.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x.vec, mean=mu.vec, sd=sqrt(sigma2.vec), log=TRUE))+
        sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
        sum(count.vec*log(pi.vec))
      return(ans)
    }
    
    #*******************************MCMC algorithm*******************************************#
    for (l in 1:(Burnin+Simsize)){
      
      #updating mu.is
      mu.vec = rnorm(q, mean=(x1.vec+m.vec[z.vec]*lambda.vec[z.vec])/(1+lambda.vec[z.vec]), 
                     sd=sqrt(sigma2.vec/(1+lambda.vec[z.vec])))
      
      #not updating sigma2.is
      
      #updating z.is
      lognig.mat = array(0,dim=c(q,k))
      for (r in 1:k){
        lognig.mat[,r] = dnorm(mu.vec,mean=m.vec[r],sd=sqrt(sigma2.vec/lambda.vec[r]),log=TRUE) + 
          dinvgamma(sigma2.vec,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
      }
      
      lognig1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + lognig.mat
      
      #need it for computational purpose. Otherwise exp(small log values) going to be 0.
      lognig2.mat = lognig1.mat - apply(lognig1.mat,1,max)
      
      nig.mat = exp(lognig2.mat)
      
      cum.nig.mat = nig.mat %*% upper.tri(diag(k), diag = TRUE)
      
      random.q.vec = runif(q)*rowSums(nig.mat)
      
      z.vec = rowSums(random.q.vec > cum.nig.mat)+1
      
      #countvec is a vector of length k where each element is number of zis in each cluster
      count.vec = tabulate(z.vec,nbins=k)
      
      #updating pi.rs
      pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
      
      #updating m.ts
      for (r in 1:k){
        if(count.vec[r] > 0){
          m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mu.vec/sigma2.vec)[which(z.vec==r)])+m.0/tau2)/
                             (lambda.vec[r]*sum((1/sigma2.vec)[which(z.vec==r)])+1/tau2),
                           sd=sqrt(1/(lambda.vec[r]*sum((1/sigma2.vec)[which(z.vec==r)])+1/tau2)))
        }
        if(count.vec[r] == 0){
          m.vec[r] = rnorm(1, mean=m.0, sd=sqrt(tau2))
        }
      }
      
      #updating lambda.rs
      for (r in 1:k){
        if(count.vec[r] > 0){
          lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                                 rate=sum((mu.vec[which(z.vec==r)]-m.vec[r])^2/
                                            (2*sigma2.vec[which(z.vec==r)]))+b.lambda)
        }
        if(count.vec[r] == 0){
          lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
        }
      }
      
      #updating alpha.rs
      alphacand.vec = rnorm(k, mean=alpha.vec, sd=1)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = loglike.alpha(alphacand.vec,beta.vec,sigma2.vec,z.vec,a.alpha,b.alpha,r)-
            loglike.alpha(alpha.vec,beta.vec,sigma2.vec,z.vec,a.alpha,b.alpha,r)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
            dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
      }
      
      #updating beta.rs
      betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS^2)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = loglike.beta(alpha.vec,betacand.vec,sigma2.vec,z.vec,a.beta,b.beta,r)-
            loglike.beta(alpha.vec,beta.vec,sigma2.vec,z.vec,a.beta,b.beta,r)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
            dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
      }
      
      d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS+MM, 
                                  lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS+MM, lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.sigma2.iter1 = d.sigma2.mixture(sigma2.grid, alpha.vec, beta.vec*SS^2, pi.vec)
      
      theta = data.frame(cbind(m.vec*SS+MM,lambda.vec,alpha.vec,beta.vec*SS^2,pi.vec))
      names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
      
      # if(l%%1000==0){
      #   theta1 = theta
      #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
      #   theta1 = theta1[order(c(-theta1$pi.vec)),]
      #   print(l)
      #   print(theta1)
      # }
      
      if (l > Burnin){
        #1st stage parameters
        mu.vec.tot = mu.vec.tot+mu.vec
        mu.vec2.tot = mu.vec2.tot+mu.vec^2
        zvec.mat[,(l-Burnin)] = z.vec
        
        #nig parameters
        theta.all = theta.all+as.matrix(theta)
        theta2.all = theta2.all+as.matrix(theta^2)
        
        #estimated density
        d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
        d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
        d.sigma2.grid.tot = d.sigma2.grid.tot + d.sigma2.iter1
        
        #total.log.like is computed based on x1.vec
        total.log.like.vec[(l-Burnin)] = total.log.like(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, 
                                                        beta.vec, z.vec, pi.vec)
        #data.log.like is compared on x.vec
        data.log.like.vec[(l-Burnin)] = data.log.like(mu.vec*SS+MM, sigma2.vec*SS^2, m.vec*SS+MM, lambda.vec, 
                                                      alpha.vec, beta.vec*SS^2, z.vec, pi.vec)
        k.vec[(l-Burnin)] = length(unique(z.vec))
      }
      
    }
    
    
    #mu.i's posterior point estimate and it's posterior variance estimate from MCMC iterations.
    mu.est = mu.vec.tot/Simsize*SS+MM
    mu.est.var = (mu.vec2.tot/Simsize-(mu.vec.tot/Simsize)^2)*SS^2
    
    
    #Summarizing density estimate from MCMC iterations. 
    d.biv.grid.est = d.biv.grid.tot/Simsize
    d.mu.grid.est = d.mu.grid.tot/Simsize
    d.sigma2.grid.est = d.sigma2.grid.tot/Simsize
    
    #Because of level switching these estimates are not reliable.
    theta.hat = data.frame(theta.all/Simsize)
    theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
    
    names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    
    order.pi.vec = order(c(-theta.hat$pi.vec))
    theta.hat = theta.hat[order.pi.vec,]
    theta.hat.var = theta.hat.var[order.pi.vec,]
    
    #Clustering analysis
    find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
    find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans2=max(tab)/length(x); rm(tab); return(ans2)}
    
    #Estimate of cluster of (mu.i, sigma.i^2).
    z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
    z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
    
    
    theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    order.pi.vec.eff = order(-theta.eff$pi.vec.est)
    theta.eff = theta.eff[order.pi.vec.eff,]
    theta.eff.var = theta.eff.var[order.pi.vec.eff,]
    
    return(list(mu.est=mu.est, mu.est.var=mu.est.var,
                z.vec.est=z.vec.est, z.vec.est.prop=z.vec.est.prop, theta.hat=theta.hat,
                theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
                d.biv.grid.est=d.biv.grid.est,
                d.mu.grid.est=d.mu.grid.est, d.sigma2.grid.est=d.sigma2.grid.est,
                total.log.like.vec=total.log.like.vec,
                data.log.like.vec=data.log.like.vec, k.vec=k.vec))
    
  }
}

# Output: mu.est = estimated value of mu_i's based on MCMC iterations.
#         mu.est.var = estimated posterior variance of mu_i.
#         z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         d.sigma2.grid.est = estimated sigma2 density on sigma2.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 

#**********************************************************************************

# X_{ij} = \mu_i+1/sqrt(prec_i^2) e_{ij}.  
# (\mu_i,prec_i) follows mixture of normal-gamma. (\mu_i,prec_i) not observed.

#**********************************************************************************

# input: x.vec = all X_{ij} in one vector of length \sum_{i=1}^q n_i.
#        ids = vector of indicators which denote the the $i$ of the corresponding X_{ij}. 
#              Range is between 1,..,q. The number i appears n_i times.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        prec.grid = grid points over which the density of unknown precision parameter will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.ng = function(x.vec, ids, gamma.pi = 0.1, k = 10,
                   Burnin = 5000, Simsize = 10000, mu.grid = c('Default'), prec.grid = c('Default'),
                   hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(truncnorm)
  orderid = order(ids)
  x.vec = x.vec[orderid]
  ids = ids[orderid]
  
  #*********************************Data Management****************************************#
  #number of replictions of mu.i. nivec[i]=5 means 5 replicates for mu.i. For balanced data all elements nivec is same.
  #For unbalanced data we need each elements to greater than equal 2 for identifiability of the joint distribution.
  ni.vec = as.data.frame(table(ids))$Freq
  q = length(ni.vec)
  
  if(sum(as.vector(ni.vec)<2)>0){
    print("Not enough data to estimate mean variance simultaneously")
  } else{
    
    #centering and scaling the whole dataset
    MM = mean(x.vec)
    SS = sd(x.vec)
    # MM = 0
    # SS = 1
    x1.vec = (x.vec-MM)/SS
    
    #******************************End of Data Management************************************#
    
    #******************************Initialize MCMC parameters************************************#
    
    #intial value of mu.i vector which is muvec. Similarly, precvec is the vector prec.i's
    xbar.vec = mu.vec = tapply(x1.vec,ids,mean)
    S2.vec = tapply(x1.vec,ids,var)
    # precision parameter
    prec.vec = 1/S2.vec
    
    #k-means clusters to determine initial values of cluster parameters.
    mscluster = kmeans(cbind(xbar.vec*SS+MM, prec.vec/SS^2),k)
    
    #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
    #initial value of mvec
    m.vec = data.frame(mscluster$center)[,1]
    
    m.vec = (m.vec-MM)/SS
    
    #latent variable which indicates cluster number. zi=r means (mu.i,prec.i^2|z.i=r) 
    #\sim NorGamma(m.r, lambda.r, alpha.r, beta.r)
    z.vec = mscluster$cluster
    
    #vector of M-H ratios of length k.
    mh.ratio1.vec = mh.ratio2.vec = rep(0,k)
    
    #matrix to store mu.i, prec.i and Z.i after each iteration.
    mu.vec.tot = prec.vec.tot = rep(0,q)
    mu.vec2.tot = prec.vec2.tot = rep(0,q)
    zvec.mat = array(0,c(q,Simsize))
    
    #vector of length of unique z.i's in each iteration.
    k.vec = vector()
    
    #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
    theta.all = array(0,c(k,5))
    theta2.all = array(0,c(k,5))
    
    #initial values alphavec, betavec, lambdavec.
    lambda.vec = rep(1,k)
    alpha_0 = (mean(prec.vec))^2/var(prec.vec)
    beta_0 = mean(prec.vec)/var(prec.vec)
    alpha.vec = rep(alpha_0,k)
    beta.vec = rep(beta_0,k)
    pi.vec = rep(1/k,k)
    
    if(sum(hyperparameters==c('Default'))>0){
      #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). 
      #So, setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
      m.0 = (mean(x.vec)-MM)/SS
      tau2 = var(xbar.vec)
      
      #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
      #lambda is generally less than 1. 
      a.lambda = 1
      b.lambda = 1
      
      #alphavec follows Gamma(a.alpha,b.alpha).
      a.alpha = 1
      b.alpha = 1
      
      #betavec follows Gamma(a.beta,b.beta). 
      a.beta = 1
      b.beta = 1/(sd(x.vec)^2)*SS^2
    }
    
    if (sum(hyperparameters!=c('Default'))>0){
      
      m.0 = (hyperparameters[1]-MM)/SS
      tau2 = hyperparameters[2]/SS^2
      a.lambda = hyperparameters[3]
      b.lambda = hyperparameters[4]
      a.alpha = hyperparameters[5]
      b.alpha = hyperparameters[6]
      a.beta = hyperparameters[7]
      b.beta = hyperparameters[8]*SS^2
    }
    
    
    #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
    #Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
    #allows higher variation in pi.r's.  gamma.pi = 0.1
    #[mu.grid,prec.grid] is a bivariate grid on which we are estimating the density.
    if (sum(mu.grid==c('Default'))>0){
      mu.grid = seq(min(xbar.vec)-IQR(xbar.vec), max(xbar.vec)+IQR(xbar.vec),
                    length.out=2*100+1)[seq(2,2*100,2)]
    }
    if (sum(prec.grid==c('Default'))>0){
      prec.grid = seq(max(min(prec.vec)-IQR(prec.vec),0.001), max(prec.vec)+IQR(prec.vec),
                      length.out=2*100+1)[seq(2,2*100,2)]
    }
    
    n1 = length(mu.grid)
    n2 = length(prec.grid)
    
    biv.grid = array(0,c(n1*n2,2))
    biv.grid[,1] = rep(mu.grid,n2)
    biv.grid[,2] = rep(prec.grid,each=n1)
    
    #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate, mu-marginal, 
    #and sigma2-marginal distribution.
    d.biv.grid.tot = rep(0,n1*n2)
    d.mu.grid.tot = rep(0,n1)
    d.prec.grid.tot = rep(0,n2)
    total.log.like.vec = data.log.like.vec = rep(0,Simsize)
    
    d.biv.mixture = function(mu.grid, prec.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dnorm(mu.grid, mean=m.vec[r], sd=sqrt(1/(prec.grid*lambda.vec[r])))*
          dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
          dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
      }
      return(ans)
    }
    
    d.prec.mixture = function(prec.grid, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    loglike.alpha = function(alpha.vec, beta.vec, prec.vec, z.vec, a.alpha, b.alpha, r){
      ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
      return(ans)
    }
    
    loglike.beta = function(alpha.vec, beta.vec, prec.vec, z.vec, a.beta, b.beta, r){
      ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
      return(ans)
    }
    
    # parameters should be in rescaled scale. As it is on x1.vec
    total.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      mu.vec.large = rep(mu.vec,ni.vec)
      prec.vec.large = rep(prec.vec,ni.vec)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x1.vec, mean=mu.vec.large, sd=sqrt(1/prec.vec.large), log=TRUE))+
        sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
        sum(count.vec*log(pi.vec))+
        sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
        sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
        sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
        sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
        lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
      
      return(ans)
    }
    
    
    # parameters should be in original scale. As it is on x.vec
    data.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      mu.vec.large = rep(mu.vec,ni.vec)
      prec.vec.large = rep(prec.vec,ni.vec)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x.vec, mean=mu.vec.large, sd=sqrt(1/prec.vec.large), log=TRUE))+
        sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
        sum(count.vec*log(pi.vec))
      return(ans)
    }
    
    #*******************************MCMC algorithm*******************************************#
    for (l in 1:(Burnin+Simsize)){
      
      #updating mu.is
      mu.vec = rnorm(q, mean=(ni.vec*xbar.vec+m.vec[z.vec]*lambda.vec[z.vec])/(ni.vec+lambda.vec[z.vec]), 
                     sd=sqrt(1/(prec.vec*(ni.vec+lambda.vec[z.vec]))))
      
      #updating sigma.is
      mu.vec.large = rep(mu.vec,ni.vec)
      
      prec.vec.shape = 0.5*(ni.vec+1)+alpha.vec[z.vec]
      
      prec.vec.rate = 0.5*tapply((x1.vec-mu.vec.large)^2,ids,sum)+
        0.5*lambda.vec[z.vec]*(mu.vec-m.vec[z.vec])^2+beta.vec[z.vec]
      
      prec.vec = rgamma(q, shape=prec.vec.shape, rate=prec.vec.rate)
      
      #updating z.is
      logng.mat = array(0,dim=c(q,k))
      for (r in 1:k){
        logng.mat[,r] = dnorm(mu.vec,mean=m.vec[r],sd=sqrt(1/(prec.vec*lambda.vec[r])),log=TRUE) + 
          dgamma(prec.vec,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
      }
      
      logng1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + logng.mat
      
      #need it for computational purpose. Otherwise exp(small log values) going to be 0.
      logng2.mat = logng1.mat - apply(logng1.mat,1,max)
      
      ng.mat = exp(logng2.mat)
      
      cum.ng.mat = ng.mat %*% upper.tri(diag(k), diag = TRUE)
      
      random.q.vec = runif(q)*rowSums(ng.mat)
      
      z.vec = rowSums(random.q.vec > cum.ng.mat)+1
      
      #countvec is a vector of length k where each element is number of zis in each cluster
      count.vec = tabulate(z.vec,nbins=k)
      
      #updating pi.rs
      pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
      
      #updating m.ts
      for (r in 1:k){
        if(count.vec[r] > 0){
          m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mu.vec*prec.vec)[which(z.vec==r)])+m.0/tau2)/
                             (lambda.vec[r]*sum((prec.vec)[which(z.vec==r)])+1/tau2),
                           sd=sqrt(1/(lambda.vec[r]*sum((prec.vec)[which(z.vec==r)])+1/tau2)))
        }
        if(count.vec[r] == 0){
          m.vec[r] = rnorm(1,mean=m.0, sd=sqrt(tau2))
        }
      }
      
      #updating lambda.rs
      for (r in 1:k){
        if(count.vec[r] > 0){
          lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                                 rate=sum((mu.vec[which(z.vec==r)]-m.vec[r])^2*
                                            prec.vec[which(z.vec==r)]/2)+b.lambda)
        }
        if(count.vec[r] == 0){
          lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
        }
      }
      
      #updating alpha.rs
      alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
      
      for (r in 1:k){
        if(count.vec[r]>0 & alphacand.vec[r]>0){
          mh.ratio1.vec[r] = loglike.alpha(alphacand.vec,beta.vec,prec.vec,z.vec,a.alpha,b.alpha,r)-
            loglike.alpha(alpha.vec,beta.vec,prec.vec,z.vec,a.alpha,b.alpha,r)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
            dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
      }
      
      #updating beta.rs
      betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS^2)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = loglike.beta(alpha.vec,betacand.vec,prec.vec,z.vec,a.beta,b.beta,r)-
            loglike.beta(alpha.vec,beta.vec,prec.vec,z.vec,a.beta,b.beta,r)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
            dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
      }
      
      d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS+MM, 
                                  lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS+MM, lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.prec.iter1 = d.prec.mixture(prec.grid, alpha.vec, beta.vec*SS^2, pi.vec)
      
      theta = data.frame(cbind(m.vec*SS+MM,lambda.vec,alpha.vec,beta.vec*SS^2,pi.vec))
      names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
      
      # if(l%%1000==0){
      #   theta1 = theta
      #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
      #   theta1 = theta1[order(c(-theta1$pi.vec)),]
      #   print(l)
      #   print(theta1)
      # }
      
      if (l > Burnin){
        #1st stage parameters
        mu.vec.tot = mu.vec.tot+mu.vec
        mu.vec2.tot = mu.vec2.tot+mu.vec^2
        prec.vec.tot = prec.vec.tot+prec.vec
        prec.vec2.tot = prec.vec2.tot+prec.vec^2
        zvec.mat[,(l-Burnin)] = z.vec
        
        #nig parameters
        theta.all = theta.all+as.matrix(theta)
        theta2.all = theta2.all+as.matrix(theta^2)
        
        #estimated density
        d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
        d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
        d.prec.grid.tot = d.prec.grid.tot + d.prec.iter1
        
        #total.log.like is computed based on x1.vec
        total.log.like.vec[(l-Burnin)] = total.log.like(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, 
                                                        beta.vec, z.vec, pi.vec)
        #data.log.like is compared on x.vec
        data.log.like.vec[(l-Burnin)] = data.log.like(mu.vec*SS+MM, prec.vec/SS^2, m.vec*SS+MM, lambda.vec, 
                                                      alpha.vec, beta.vec*SS^2, z.vec, pi.vec)
        k.vec[(l-Burnin)] = length(unique(z.vec))
      }
      
    }
    
    #mu.i and sigma.i^2's posterior point estimate and it's posterior variance estimate from MCMC iterations.
    mu.est = mu.vec.tot/Simsize*SS+MM
    mu.est.var = (mu.vec2.tot/Simsize-(mu.vec.tot/Simsize)^2)*SS^2
    prec.est = (prec.vec.tot/Simsize)/SS^2
    prec.est.var = (prec.vec2.tot/Simsize-(prec.vec.tot/Simsize)^2)/SS^4
    
    #Summarizing density estimate from MCMC iterations. 
    d.biv.grid.est = d.biv.grid.tot/Simsize
    d.mu.grid.est = d.mu.grid.tot/Simsize
    d.prec.grid.est = d.prec.grid.tot/Simsize
    
    #Because of level switching these estimatea are not reliable.
    theta.hat = data.frame(theta.all/Simsize)
    theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
    
    names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    
    order.pi.vec = order(c(-theta.hat$pi.vec))
    theta.hat = theta.hat[order.pi.vec,]
    theta.hat.var = theta.hat.var[order.pi.vec,]
    
    #Clustering analysis
    find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
    find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans2=max(tab)/length(x); rm(tab); return(ans2)}
    
    #Estimate of cluster of (mu.i,sigma.i^2).
    z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
    z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
    
    theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    order.pi.vec.eff = order(-theta.eff$pi.vec.est)
    theta.eff = theta.eff[order.pi.vec.eff,]
    theta.eff.var = theta.eff.var[order.pi.vec.eff,]
    
    return(list(mu.est=mu.est, mu.est.var=mu.est.var, prec.est=prec.est, 
                prec.est.var=prec.est.var, 
                z.vec.est=z.vec.est, z.vec.est.prop=z.vec.est.prop, theta.hat=theta.hat,
                theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
                d.biv.grid.est=d.biv.grid.est,
                d.mu.grid.est=d.mu.grid.est, d.prec.grid.est=d.prec.grid.est,
                total.log.like.vec=total.log.like.vec,
                data.log.like.vec=data.log.like.vec, k.vec=k.vec))
    
  }
}


# Output: mu.est = estimated value of mu_i's based on MCMC iterations.
#         mu.est.var = estimated posterior variance of mu_i.
#         prec.est = estimated value of prec_i's based on MCMC iterations.
#         prec.est.var = estimated posterior variance of prec_i.
#         z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         d.prec.grid.est = estimated prec density on prec.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 

#**********************************************************************************

# X_{i} = \mu_i+1/sqrt(prec_i^2) e_{i}. (\mu_i,prec_i) follows mixture of normal-gamma. \mu_i not observed.
# prec_i is observed.

#**********************************************************************************

# input: x.vec = all X_{i} in one vector of length q.
#        precd = vector of known precision parameter of length q.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.ng.sknown = function(x.vec, precd, gamma.pi = 0.1, k = 10,
                          Burnin = 5000, Simsize = 10000, mu.grid = c('Default'), prec.grid = c('Default'),
                          hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(truncnorm)
  
  
  #*********************************Data Management****************************************#
  q = length(x.vec)
  
  if(length(x.vec) != length(precd)){
    print("Dimension of x.vec and sigma2d is not same")
  } else{
    
    #centering and scaling the whole dataset
    MM = mean(x.vec)
    SS = sd(x.vec)
    # MM = 0
    # SS = 1
    x1.vec = (x.vec-MM)/SS
    precd1 = precd*SS^2
    
    #******************************End of Data Management************************************#
    
    #******************************Initialize MCMC parameters************************************#
    
    #intial value of mu.i vector which is muvec. Similarly, precvec is the vector prec.i's
    mu.vec = x1.vec
    prec.vec = precd1
    # precision parameter
    
    #k-means clusters to determine initial values of cluster parameters.
    mscluster = kmeans(cbind(x.vec, precd),k)
    
    #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
    #initial value of mvec
    m.vec = data.frame(mscluster$center)[,1]
    
    m.vec = (m.vec-MM)/SS
    
    #latent variable which indicates cluster number. zi=r means (mu.i,prec.i^2|z.i=r) 
    #\sim NorGamma(m.r, lambda.r, alpha.r, beta.r)
    z.vec = mscluster$cluster
    
    #vector of M-H ratios of length k.
    mh.ratio1.vec = mh.ratio2.vec = rep(0,k)
    
    #matrix to store mu.i, prec.i and Z.i after each iteration.
    mu.vec.tot = rep(0,q)
    mu.vec2.tot = rep(0,q)
    zvec.mat = array(0,c(q,Simsize))
    
    #vector of length of unique z.i's in each iteration.
    k.vec = vector()
    
    #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
    theta.all = array(0,c(k,5))
    theta2.all = array(0,c(k,5))
    
    #initial values alphavec, betavec, lambdavec.
    lambda.vec = rep(1,k)
    alpha_0 = (mean(prec.vec))^2/var(prec.vec)
    beta_0 = mean(prec.vec)/var(prec.vec)
    alpha.vec = rep(alpha_0,k)
    beta.vec = rep(beta_0,k)
    pi.vec = rep(1/k,k)
    
    if(sum(hyperparameters==c('Default'))>0){
      #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). 
      #So, setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
      m.0 = (mean(x.vec)-MM)/SS
      tau2 = var(x.vec)
      
      #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
      #lambda is generally less than 1. 
      a.lambda = 1
      b.lambda = 1
      
      #alphavec follows Gamma(a.alpha,b.alpha).
      a.alpha = 1
      b.alpha = var(precd1)/mean(precd1)^2
      
      #betavec follows Gamma(a.beta,b.beta).
      a.beta = 1
      b.beta = var(precd1)/mean(precd1)
    }
    
    if (sum(hyperparameters!=c('Default'))>0){
      
      m.0 = (hyperparameters[1]-MM)/SS
      tau2 = hyperparameters[2]/SS^2
      a.lambda = hyperparameters[3]
      b.lambda = hyperparameters[4]
      a.alpha = hyperparameters[5]
      b.alpha = hyperparameters[6]
      a.beta = hyperparameters[7]
      b.beta = hyperparameters[8]*SS^2
    }
    
    
    #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k 
    #follows Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
    #allows higher variation in pi.r's.  gamma.pi = 0.1
    #[mu.grid,prec.grid] is a bivariate grid on which we are estimating the density.
    if (sum(mu.grid==c('Default'))>0){
      mu.grid = seq(min(x.vec)-IQR(x.vec), max(x.vec)+IQR(x.vec),
                    length.out=2*100+1)[seq(2,2*100,2)]
    }
    if (sum(prec.grid==c('Default'))>0){
      prec.grid = seq(max(min(prec.vec)-IQR(prec.vec),0.001), max(prec.vec)+IQR(prec.vec),
                      length.out=2*100+1)[seq(2,2*100,2)]
    }
    
    n1 = length(mu.grid)
    n2 = length(prec.grid)
    
    biv.grid = array(0,c(n1*n2,2))
    biv.grid[,1] = rep(mu.grid,n2)
    biv.grid[,2] = rep(prec.grid,each=n1)
    
    #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate,
    #mu-marginal, and sigma2-marginal distribution.
    d.biv.grid.tot = rep(0,n1*n2)
    d.mu.grid.tot = rep(0,n1)
    d.prec.grid.tot = rep(0,n2)
    total.log.like.vec = data.log.like.vec = rep(0,Simsize)
    
    d.biv.mixture = function(mu.grid, prec.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dnorm(mu.grid, mean=m.vec[r], sd=sqrt(1/(prec.grid*lambda.vec[r])))*
          dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
          dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
      }
      return(ans)
    }
    
    d.prec.mixture = function(prec.grid, alpha.vec, beta.vec, pi.vec){
      ans = 0
      k.comp = length(pi.vec)
      for (r in 1:k.comp){
        ans = ans + pi.vec[r]*dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
      }
      return(ans)
    }
    
    loglike.alpha = function(alpha.vec, beta.vec, prec.vec, z.vec, a.alpha, b.alpha, r){
      ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
      return(ans)
    }
    
    loglike.beta = function(alpha.vec, beta.vec, prec.vec, z.vec, a.beta, b.beta, r){
      ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
        dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
      return(ans)
    }
    
    # parameters should be in rescaled scale. As it is on x1.vec
    total.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x1.vec, mean=mu.vec, sd=sqrt(1/prec.vec), log=TRUE))+
        sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
        sum(count.vec*log(pi.vec))+
        sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
        sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
        sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
        sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
        lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
      
      return(ans)
    }
    
    
    # parameters should be in original scale. As it is on x.vec
    data.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
      count.vec = tabulate(z.vec,nbins=k)
      m.vec.large = m.vec[z.vec]
      lambda.vec.large = lambda.vec[z.vec]
      alpha.vec.large = alpha.vec[z.vec]
      beta.vec.large = beta.vec[z.vec]
      
      ans = sum(dnorm(x.vec, mean=mu.vec, sd=sqrt(1/prec.vec), log=TRUE))+
        sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
        sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
        sum(count.vec*log(pi.vec))
      return(ans)
    }
    
    #*******************************MCMC algorithm*******************************************#
    for (l in 1:(Burnin+Simsize)){
      
      #updating mu.is
      mu.vec = rnorm(q, mean=(x1.vec+m.vec[z.vec]*lambda.vec[z.vec])/(1+lambda.vec[z.vec]), 
                     sd=sqrt(1/(prec.vec*(1+lambda.vec[z.vec]))))
      
      
      
      #updating z.is
      logng.mat = array(0,dim=c(q,k))
      for (r in 1:k){
        logng.mat[,r] = dnorm(mu.vec,mean=m.vec[r],sd=sqrt(1/(prec.vec*lambda.vec[r])),log=TRUE) + 
          dgamma(prec.vec,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
      }
      
      logng1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + logng.mat
      
      #need it for computational purpose. Otherwise exp(small log values) going to be 0.
      logng2.mat = logng1.mat - apply(logng1.mat,1,max)
      
      ng.mat = exp(logng2.mat)
      
      cum.ng.mat = ng.mat %*% upper.tri(diag(k), diag = TRUE)
      
      random.q.vec = runif(q)*rowSums(ng.mat)
      
      z.vec = rowSums(random.q.vec > cum.ng.mat)+1
      
      #countvec is a vector of length k where each element is number of zis in each cluster
      count.vec = tabulate(z.vec,nbins=k)
      
      #updating pi.rs
      pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
      
      #updating m.ts
      for (r in 1:k){
        if(count.vec[r] > 0){
          m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mu.vec*prec.vec)[which(z.vec==r)])+m.0/tau2)/
                             (lambda.vec[r]*sum((prec.vec)[which(z.vec==r)])+1/tau2),
                           sd=sqrt(1/(lambda.vec[r]*sum((prec.vec)[which(z.vec==r)])+1/tau2)))
        }
        if(count.vec[r] == 0){
          m.vec[r] = rnorm(1,mean=m.0, sd=sqrt(tau2))
        }
      }
      
      #updating lambda.rs
      for (r in 1:k){
        if(count.vec[r] > 0){
          lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                                 rate=sum((mu.vec[which(z.vec==r)]-m.vec[r])^2*
                                            prec.vec[which(z.vec==r)]/2)+b.lambda)
        }
        if(count.vec[r] == 0){
          lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
        }
      }
      
      #updating alpha.rs
      alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
      
      for (r in 1:k){
        if(count.vec[r]>0 & alphacand.vec[r]>0){
          mh.ratio1.vec[r] = loglike.alpha(alphacand.vec,beta.vec,prec.vec,z.vec,a.alpha,b.alpha,r)-
            loglike.alpha(alpha.vec,beta.vec,prec.vec,z.vec,a.alpha,b.alpha,r)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & alphacand.vec[r] > 0){
          mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
            dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
          if(mh.ratio1.vec[r] > log(runif(1))){
            alpha.vec[r] = alphacand.vec[r]
          }
        }
      }
      
      #updating beta.rs
      betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS^2)
      
      for (r in 1:k){
        if(count.vec[r] > 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = loglike.beta(alpha.vec,betacand.vec,prec.vec,z.vec,a.beta,b.beta,r)-
            loglike.beta(alpha.vec,beta.vec,prec.vec,z.vec,a.beta,b.beta,r)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
        if(count.vec[r] == 0 & betacand.vec[r] > 0){
          mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
            dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
          if(mh.ratio2.vec[r] > log(runif(1))){
            beta.vec[r] = betacand.vec[r]
          }
        }
      }
      
      d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS+MM, 
                                  lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS+MM, lambda.vec, alpha.vec, beta.vec*SS^2, pi.vec)
      
      d.prec.iter1 = d.prec.mixture(prec.grid, alpha.vec, beta.vec*SS^2, pi.vec)
      
      theta = data.frame(cbind(m.vec*SS+MM,lambda.vec,alpha.vec,beta.vec*SS^2,pi.vec))
      names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
      
      # if(l%%1000==0){
      #   theta1 = theta
      #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
      #   theta1 = theta1[order(c(-theta1$pi.vec)),]
      #   print(l)
      #   print(theta1)
      # }
      
      if (l > Burnin){
        #1st stage parameters
        mu.vec.tot = mu.vec.tot+mu.vec
        mu.vec2.tot = mu.vec2.tot+mu.vec^2
        zvec.mat[,(l-Burnin)] = z.vec
        
        #nig parameters
        theta.all = theta.all+as.matrix(theta)
        theta2.all = theta2.all+as.matrix(theta^2)
        
        #estimated density
        d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
        d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
        d.prec.grid.tot = d.prec.grid.tot + d.prec.iter1
        
        #total.log.like is computed based on x1.vec
        total.log.like.vec[(l-Burnin)] = total.log.like(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, 
                                                        beta.vec, z.vec, pi.vec)
        #data.log.like is compared on x.vec
        data.log.like.vec[(l-Burnin)] = data.log.like(mu.vec*SS+MM, prec.vec*SS^2, m.vec*SS+MM, lambda.vec, 
                                                      alpha.vec, beta.vec*SS^2, z.vec, pi.vec)
        k.vec[(l-Burnin)] = length(unique(z.vec))
      }
      
    }
    
    #mu.i, sigma.i^2's posterior point estimate and it's posterior variance estimate from MCMC iterations.
    mu.est = mu.vec.tot/Simsize*SS+MM
    mu.est.var = (mu.vec2.tot/Simsize-(mu.vec.tot/Simsize)^2)*SS^2
    
    
    #Summarizing density estimate from MCMC iterations. 
    d.biv.grid.est = d.biv.grid.tot/Simsize
    d.mu.grid.est = d.mu.grid.tot/Simsize
    d.prec.grid.est = d.prec.grid.tot/Simsize
    
    #Because of level switching these estimatea are not reliable.
    theta.hat = data.frame(theta.all/Simsize)
    theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
    
    names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
    
    order.pi.vec = order(c(-theta.hat$pi.vec))
    theta.hat = theta.hat[order.pi.vec,]
    theta.hat.var = theta.hat.var[order.pi.vec,]
    
    #Clustering analysis
    find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
    find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
    ans2=max(tab)/length(x); rm(tab); return(ans2)}
    
    #Estimate of cluster of (mu.i,sigma.i^2).
    z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
    z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
    
    theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
    order.pi.vec.eff = order(-theta.eff$pi.vec.est)
    theta.eff = theta.eff[order.pi.vec.eff,]
    theta.eff.var = theta.eff.var[order.pi.vec.eff,]
    
    return(list(mu.est=mu.est, mu.est.var=mu.est.var,
                z.vec.est=z.vec.est, z.vec.est.prop=z.vec.est.prop, theta.hat=theta.hat,
                theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
                d.biv.grid.est=d.biv.grid.est,
                d.mu.grid.est=d.mu.grid.est, d.prec.grid.est=d.prec.grid.est,
                total.log.like.vec=total.log.like.vec,
                data.log.like.vec=data.log.like.vec, k.vec=k.vec))
    
  }
}

# Output: mu.est = estimated value of mu_i's based on MCMC iterations.
#         mu.est.var = estimated posterior variance of mu_i.
#         z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         d.prec.grid.est = estimated prec density on sigma2.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 

#**********************************************************************************

# (\mu_i,\sigma_i^2) follows mixture of normal-inverse-gamma. (\mu_i,\sigma_i^2) observed.

#**********************************************************************************

# input: mud = vector of length q where ith element is \mu_i which is known.
#        sigma2d = vector of length q where ith element is \sigma_i^2 which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)


DPMM.nig.msknown = function(mud, sigma2d, gamma.pi = 0.1, k = 10, Burnin = 5000, Simsize = 10000, 
                            mu.grid=c('Default'), sigma2.grid=c('Default'), hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("invgamma")) install.packages("invgamma")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(invgamma)
  # library(truncnorm)
  
  # rinvgamma = invgamma::rinvgamma
  # dinvgamma = invgamma::dinvgamma
  
  #*********************************Data Management****************************************#
  
  q = length(mud)
  
  #centering and scaling the whole dataset
  MM1 = mean(mud)
  SS1 = sd(mud)
  SS2 = sd(sigma2d)
  # MM1 = 0
  # SS1 = 1
  # SS2 = 1
  
  mud1 = (mud-MM1)/SS1
  sigma2d1 = sigma2d/SS2
  
  #******************************End of Data Management************************************#
  
  #******************************Initialize MCMC parameters************************************#
  
  #k-means clusters to determine initial values of cluster parameters.
  mscluster = kmeans(cbind(mud,sigma2d),k)
  
  #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
  #initial value of mvec
  m.vec = (data.frame(mscluster$center)$mud-MM1)/SS1
  
  #latent variable which indicates cluster number. zi=r means (mu.i,sigma.i^2|z.i=r) 
  #\sim NIG(m.r, lambda.r, alpha.r, beta.r)
  z.vec = mscluster$cluster
  
  #vector of M-H ratios of length k.
  mh.ratio1.vec = mh.ratio2.vec = vector()
  
  #matrix to store Z.i after each iteration.
  zvec.mat = array(0,c(q,Simsize))
  
  #vector of length of unique z.i's in each iteration.
  k.vec = vector()
  
  #matrix of m.r's, alpha.r's, beta.r's, lambda.r's where r=1,...,k. Empty matrix to store each iteration.
  #mvec.mat = lambdavec.mat = pivec.mat = alphavec.mat = betavec.mat = array(0,c(k,Simsize))
  
  #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
  theta.all = array(0,c(k,5))
  theta2.all = array(0,c(k,5))
  
  #initial values alphavec, betavec, lambdavec all taken as 1, without trying to estimate from data.
  lambda.vec = rep(1,k)
  alpha_0 = (mean(sigma2d1))^2/var(sigma2d1)+2
  beta_0 = ((mean(sigma2d1))^2/var(sigma2d1)+1)*mean(sigma2d1)
  alpha.vec = rep(alpha_0,k)
  beta.vec = rep(beta_0,k)
  pi.vec = rep(1/k,k)
  
  if(sum(hyperparameters==c('Default'))>0){
    #Hyperprior parameters. mvec follows N(m.0, \tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). So, 
    #setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
    # m.0 = mean(x1.vec)
    m.0 = (mean(mud)-MM1)/SS1
    tau2 = var(mud)/SS1^2
    
    #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
    #lambda is generally less than 1. 
    a.lambda = 1
    b.lambda = 1
    
    #alphavec follows Gamma(a.alpha,b.alpha). 
    a.alpha = 1
    b.alpha = 1
    
    #betavec follows Gamma(a.beta,b.beta). 
    a.beta = 1
    b.beta = 1/(sd(sigma2d)^2)*SS2^2
  }
  
  if (sum(hyperparameters!=c('Default'))>0){
    
    m.0 = (hyperparameters[1]-MM1)/SS1
    tau2 = hyperparameters[2]/SS1^2
    a.lambda = hyperparameters[3]
    b.lambda = hyperparameters[4]
    a.alpha = hyperparameters[5]
    b.alpha = hyperparameters[6]
    a.beta = hyperparameters[7]
    b.beta = hyperparameters[8]*SS2^2
  }
  
  #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
  #Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
  #allows higher variation in pi.r's. gamma.pi = 0.1
  if (sum(mu.grid==c('Default'))>0){
    mu.grid = seq(min(mud)-IQR(mud), max(mud)+IQR(mud),
                  length.out=2*100+1)[seq(2,2*100,2)]
  }
  if (sum(sigma2.grid==c('Default'))>0){
    sigma2.grid = seq(max(min(sigma2d)-IQR(sigma2d),0.001), max(sigma2d)+IQR(sigma2d),
                      length.out=2*100+1)[seq(2,2*100,2)]
  }
  
  #[mu.grid,sigma2.grid] is a bivariate grid on which we are estimating the density.
  n1 = length(mu.grid)
  n2 = length(sigma2.grid)
  
  biv.grid = array(0,c(n1*n2,2))
  biv.grid[,1] = rep(mu.grid,n2)
  biv.grid[,2] = rep(sigma2.grid,each=n1)
  
  #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate, 
  #mu-marginal, and sigma2-marginal distribution.
  d.biv.grid.tot = rep(0,n1*n2)
  d.mu.grid.tot = rep(0,n1)
  d.sigma2.grid.tot = rep(0,n2)
  total.log.like.vec = data.log.like.vec = vector()
  
  d.biv.mixture = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dnorm(mu.vec, mean=m.vec[r], sd=sqrt(sigma2.vec/lambda.vec[r]))*
        dinvgamma(sigma2.vec, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
        dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
    }
    return(ans)
  }
  
  d.sigma2.mixture = function(sigma2.grid, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  loglike.alpha = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.alpha, b.alpha, r){
    ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
    return(ans)
  }
  
  loglike.beta = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.beta, b.beta,r){
    ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
    return(ans)
  }
  
  # parameters should be in rescaled scale. As it is on mud1 and sigma2d1
  total.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
      sum(count.vec*log(pi.vec))+
      sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
      sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
      sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
      sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
      lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
    
    return(ans)
  }
  
  
  # parameters should be in original scale. As it is on mud and sigma2d
  data.log.like = function(mu.vec, sigma2.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(sigma2.vec/lambda.vec.large), log=TRUE))+
      sum(count.vec*log(pi.vec))
    return(ans)
  }
  
  #*******************************MCMC algorithm*******************************************#
  
  for (l in 1:(Burnin+Simsize)){
    
    #updating z.is
    lognig.mat = array(0,dim=c(q,k))
    for (r in 1:k){
      lognig.mat[,r] = dnorm(mud1,mean=m.vec[r],sd=sqrt(sigma2d1/lambda.vec[r]),log=TRUE) + 
        dinvgamma(sigma2d1,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
    }
    
    lognig1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + lognig.mat
    
    #need it for computational purpose. Otherwise exp(small log values) going to be 0.
    lognig2.mat = lognig1.mat - apply(lognig1.mat,1,max)
    
    nig.mat = exp(lognig2.mat)
    
    cum.nig.mat = nig.mat %*% upper.tri(diag(k), diag = TRUE)
    
    random.q.vec = runif(q)*rowSums(nig.mat)
    
    z.vec = rowSums(random.q.vec>cum.nig.mat)+1
    
    #countvec is a vector of length k where each element is number of zis in each cluster
    count.vec = tabulate(z.vec,nbins=k)
    
    #updating pi.rs
    pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
    
    #updating m.ts
    for (r in 1:k){
      if(count.vec[r] > 0){
        m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mud1/sigma2d1)[which(z.vec==r)])+m.0/tau2)/
                           (lambda.vec[r]*sum((1/sigma2d1)[which(z.vec==r)])+1/tau2),
                         sd=sqrt(1/(lambda.vec[r]*sum((1/sigma2d1)[which(z.vec==r)])+1/tau2)))
      }
      if(count.vec[r] == 0){
        m.vec[r] = rnorm(1,mean=m.0,sd=sqrt(tau2))
      }
    }
    
    #updating lambda.rs
    for (r in 1:k){
      if(count.vec[r] > 0){
        lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                               rate=sum((mud1[which(z.vec==r)]-m.vec[r])^2/
                                          (2*sigma2d1[which(z.vec==r)]))+b.lambda)
      }
      if(count.vec[r] == 0){
        lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
      }
    }
    
    #updating alpha.rs
    alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
    
    for (r in 1:k){
      if(count.vec[r] > 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = loglike.alpha(alphacand.vec,beta.vec,sigma2d1,z.vec,a.alpha,b.alpha,r)-
          loglike.alpha(alpha.vec,beta.vec,sigma2d1,z.vec,a.alpha,b.alpha,r)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
          dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
    }
    
    #updating beta.rs
    betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS2^2)
    
    for (r in 1:k){
      if(length(which(z.vec==r)) > 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = loglike.beta(alpha.vec,betacand.vec,sigma2d1,z.vec,a.beta,b.beta,r)-
          loglike.beta(alpha.vec,beta.vec,sigma2d1,z.vec,a.beta,b.beta,r)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
          dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
    }
    
    d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS1+MM1, 
                                lambda.vec, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS1+MM1, lambda.vec, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    d.sigma2.iter1 = d.sigma2.mixture(sigma2.grid, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    theta = data.frame(cbind(m.vec*SS1+MM1,lambda.vec*(SS2/SS1)^2,alpha.vec,beta.vec*SS2^2,pi.vec))
    names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
    
    # if(l%%1000==0){
    #   theta1 = theta
    #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
    #   theta1 = theta1[order(c(-theta1$pi.vec)),]
    #   print(l)
    #   print(theta1)
    # }
    
    if (l > Burnin){
      #1st stage parameters
      zvec.mat[,(l-Burnin)] = z.vec
      #main parameters
      theta.all = theta.all+as.matrix(theta)
      theta2.all = theta.all+as.matrix(theta^2)
      
      #nig parameters
      d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
      d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
      d.sigma2.grid.tot = d.sigma2.grid.tot + d.sigma2.iter1
      
      #total.log.like is computed based on mud1 and sigma2d1
      total.log.like.vec[(l-Burnin)] = total.log.like(mud1, sigma2d1, m.vec, lambda.vec, alpha.vec, 
                                                      beta.vec, z.vec, pi.vec)
      #data.log.like is compared on mud and sigma2d
      data.log.like.vec[(l-Burnin)] = data.log.like(mud, sigma2d, m.vec*SS1+MM1, lambda.vec*(SS2/SS1)^2, 
                                                    alpha.vec, beta.vec*SS2^2, z.vec, pi.vec)
      k.vec[(l-Burnin)] = length(unique(z.vec))
    }
    
  }
  
  #Summarizing density estimate from MCMC iterations.
  d.biv.grid.est = d.biv.grid.tot/Simsize
  d.mu.grid.est = d.mu.grid.tot/Simsize
  d.sigma2.grid.est = d.sigma2.grid.tot/Simsize
  
  #Because of level switching these estimatea are not reliable.
  theta.hat = data.frame(theta.all/Simsize)
  theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
  
  names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
  names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
  
  order.pi.vec = order(c(-theta.hat$pi.vec))
  theta.hat = theta.hat[order.pi.vec,]
  theta.hat.var = theta.hat.var[order.pi.vec,]
  
  #Clustering analysis
  find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
  find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans2=max(tab)/length(x); rm(tab); return(ans2)}
  
  z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
  z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
  
  theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  order.pi.vec.eff = order(-theta.eff$pi.vec.est)
  theta.eff = theta.eff[order.pi.vec.eff,]
  theta.eff.var = theta.eff.var[order.pi.vec.eff,]
  
  return(list(z.vec.est=z.vec.est, 
              z.vec.est.prop=z.vec.est.prop, 
              theta.hat=theta.hat, theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
              d.biv.grid.est=d.biv.grid.est,
              d.mu.grid.est=d.mu.grid.est, d.sigma2.grid.est=d.sigma2.grid.est, 
              total.log.like.vec=total.log.like.vec, 
              data.log.like.vec=data.log.like.vec, k.vec=k.vec))
  
}

# Output: z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 


#**********************************************************************************

#(\mu_i,prec_i) follows mixture of normal-gamma. (\mu_i,prec_i) not observed.

#**********************************************************************************

# input: mud = vector of length q where ith element is \mu_i which is known.
#        precd = vector of length q where ith element is \prec_i which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        prec.grid = grid points over which the density of precision will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.ng.msknown = function(mud, precd, gamma.pi = 0.1, k = 10, Burnin = 5000, Simsize = 10000, 
                           mu.grid = c('Default'), prec.grid = c('Default'), hyperparameters = c('Default')){
  
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(truncnorm)
  
  #*********************************Data Management****************************************#
  
  q = length(mud)
  
  #centering and scaling the whole dataset
  MM1 = mean(mud)
  SS1 = sd(mud)
  SS2 = sd(precd)
  # MM1 = 0
  # SS1 = 1
  # SS2 = 1
  
  mud1 = (mud-MM1)/SS1
  precd1 = precd/SS2
  
  #******************************End of Data Management************************************#
  
  #******************************Initialize MCMC parameters************************************#
  
  #k-means clusters to determine initial values of cluster parameters.
  mscluster = kmeans(cbind(mud,precd),k)
  
  #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
  #initial value of mvec
  m.vec = (data.frame(mscluster$center)$mud-MM1)/SS1
  
  #latent variable which indicates cluster number. zi=r means (mu.i,prec.i|z.i=r) 
  #\sim NorGamma(m.r, lambda.r, alpha.r, beta.r)
  z.vec = mscluster$cluster
  
  #vector of M-H ratios of length k.
  mh.ratio1.vec = mh.ratio2.vec = vector()
  
  #matrix to store Z.i after each iteration.
  zvec.mat = array(0,c(q,Simsize))
  
  #vector of length of unique z.i's in each iteration.
  k.vec = vector()
  
  #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
  theta.all = array(0,c(k,5))
  theta2.all = array(0,c(k,5))
  
  #initial values alphavec, betavec, lambdavec all taken as 1, without trying to estimate from data.
  lambda.vec = rep(1,k)
  alpha_0 = (mean(precd1))^2/var(precd1)
  beta_0 = mean(precd1)/var(precd1)
  alpha.vec = rep(alpha_0,k)
  beta.vec = rep(beta_0,k)
  pi.vec = rep(1/k,k)
  
  if(sum(hyperparameters==c('Default'))>0){
    #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). So, 
    #setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
    # m.0 = mean(x1.vec)
    m.0 = (mean(mud)-MM1)/SS1
    tau2 = var(mud)/SS1^2
    
    #lambdavec follows Gamma(a.lambda,b.lambda). lambda is ratio of within group variance/ between group variance. 
    #lambda is generally less than 1. 
    a.lambda = 1
    b.lambda = 1
    
    #alphavec follows Gamma(a.alpha,b.alpha). 
    a.alpha = 1
    b.alpha = 1
    
    #betavec follows Gamma(a.beta,b.beta). 
    a.beta = 1
    b.beta = 1/(sd(precd)^2)*SS2^2
  }
  
  if (sum(hyperparameters!=c('Default'))>0){
    
    m.0 = (hyperparameters[1]-MM1)/SS1
    tau2 = hyperparameters[2]/SS1^2
    a.lambda = hyperparameters[3]
    b.lambda = hyperparameters[4]
    a.alpha = hyperparameters[5]
    b.alpha = hyperparameters[6]
    a.beta = hyperparameters[7]
    b.beta = hyperparameters[8]*SS2^2
  }
  
  #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
  #Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
  #allows higher variation in pi.r's. gamma.pi = 0.1
  if (sum(mu.grid==c('Default'))>0){
    mu.grid = seq(min(mud)-IQR(mud), max(mud)+IQR(mud),
                  length.out=2*100+1)[seq(2,2*100,2)]
  }
  if (sum(prec.grid==c('Default'))>0){
    prec.grid = seq(max(min(precd)-IQR(precd),0.001), max(precd)+IQR(precd),
                    length.out=2*100+1)[seq(2,2*100,2)]
  }
  
  #[mu.grid,prec.grid] is a bivariate grid on which we are estimating the density.
  n1 = length(mu.grid)
  n2 = length(prec.grid)
  
  biv.grid = array(0,c(n1*n2,2))
  biv.grid[,1] = rep(mu.grid,n2)
  biv.grid[,2] = rep(prec.grid,each=n1)
  
  #arrays to store density value at each grid points. biv, mu, sigma2 represents bivariate, 
  #mu-marginal, and sigma2-marginal distribution.
  d.biv.grid.tot = rep(0,n1*n2)
  d.mu.grid.tot = rep(0,n1)
  d.prec.grid.tot = rep(0,n2)
  
  total.log.like.vec = data.log.like.vec = vector()
  
  d.biv.mixture = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dnorm(mu.vec, mean=m.vec[r], sd=sqrt(1/(prec.vec*lambda.vec[r])))*
        dgamma(prec.vec, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  d.mu.mixture = function(mu.grid, m.vec, lambda.vec, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*
        dt(sqrt(lambda.vec[r]*alpha.vec[r]/beta.vec[r])*(mu.grid-m.vec[r]), df=2*alpha.vec[r])
    }
    return(ans)
  }
  
  d.prec.mixture = function(prec.grid, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  loglike.alpha = function(alpha.vec, beta.vec, prec.vec, z.vec, a.alpha, b.alpha, r){
    ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
    return(ans)
  }
  
  loglike.beta = function(alpha.vec, beta.vec, prec.vec, z.vec, a.beta, b.beta,r){
    ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
    return(ans)
  }
  
  #parameters should be in rescaled scale. As it is on mud1 and precd1
  total.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
      sum(count.vec*log(pi.vec))+
      sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
      sum(dgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
      sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
      sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
      lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
    
    return(ans)
  }
  
  
  #parameters should be in original scale. As it is on mud and precd
  data.log.like = function(mu.vec, prec.vec, m.vec, lambda.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(1/(prec.vec*lambda.vec.large)), log=TRUE))+
      sum(count.vec*log(pi.vec))
    return(ans)
  }
  
  #*******************************MCMC algorithm*******************************************#
  
  for (l in 1:(Burnin+Simsize)){
    
    #updating z.is
    logng.mat = array(0,dim=c(q,k))
    for (r in 1:k){
      logng.mat[,r] = dnorm(mud1,mean=m.vec[r],sd=sqrt(1/(precd1*lambda.vec[r])),log=TRUE) + 
        dgamma(precd1,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
    }
    
    logng1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + logng.mat
    
    #need it for computational purpose. Otherwise exp(small log values) going to be 0.
    logng2.mat = logng1.mat - apply(logng1.mat,1,max)
    
    ng.mat = exp(logng2.mat)
    
    cum.ng.mat = ng.mat %*% upper.tri(diag(k), diag = TRUE)
    
    random.q.vec = runif(q)*rowSums(ng.mat)
    
    z.vec = rowSums(random.q.vec>cum.ng.mat)+1
    
    #countvec is a vector of length k where each element is number of zis in each cluster
    count.vec = tabulate(z.vec,nbins=k)
    
    #updating pi.rs
    pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
    
    #updating m.ts
    for (r in 1:k){
      if(count.vec[r] > 0){
        m.vec[r] = rnorm(1,mean=(lambda.vec[r]*sum((mud1*precd1)[which(z.vec==r)])+m.0/tau2)/
                           (lambda.vec[r]*sum(precd1[which(z.vec==r)])+1/tau2),
                         sd=sqrt(1/(lambda.vec[r]*sum(precd1[which(z.vec==r)])+1/tau2)))
      }
      if(count.vec[r] == 0){
        m.vec[r] = rnorm(1,mean=m.0,sd=sqrt(tau2))
      }
    }
    
    #updating lambda.rs
    for (r in 1:k){
      if(count.vec[r] > 0){
        lambda.vec[r] = rgamma(1,shape=count.vec[r]/2+a.lambda, 
                               rate=sum((mud1[which(z.vec==r)]-m.vec[r])^2*
                                          prec.vec[which(z.vec==r)]/2)+b.lambda)
      }
      if(count.vec[r] == 0){
        lambda.vec[r] = rgamma(1,shape=a.lambda, rate=b.lambda)
      }
    }
    
    #updating alpha.rs
    alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
    
    for (r in 1:k){
      if(count.vec[r] > 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = loglike.alpha(alphacand.vec, beta.vec, precd1, z.vec, a.alpha, b.alpha, r)-
          loglike.alpha(alpha.vec, beta.vec, precd1, z.vec, a.alpha, b.alpha, r)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
          dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
    }
    
    #updating beta.rs
    betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS2^2)
    
    for (r in 1:k){
      if(length(which(z.vec==r)) > 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = loglike.beta(alpha.vec, betacand.vec, precd1,z.vec, a.beta, b.beta, r)-
          loglike.beta(alpha.vec, beta.vec, precd1, z.vec, a.beta, b.beta, r)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
          dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
    }
    
    d.biv.iter1 = d.biv.mixture(biv.grid[,1], biv.grid[,2], m.vec*SS1+MM1, 
                                lambda.vec, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS1+MM1, lambda.vec, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    d.prec.iter1 = d.prec.mixture(prec.grid, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    theta = data.frame(cbind(m.vec*SS1+MM1, lambda.vec*(SS2/SS1)^2, alpha.vec, beta.vec*SS2^2, pi.vec))
    names(theta) = c("m.vec","lambda.vec","alpha.vec","beta.vec","pi.vec")
    
    # if(l%%1000==0){
    #   theta1 = theta
    #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
    #   theta1 = theta1[order(c(-theta1$pi.vec)),]
    #   print(l)
    #   print(theta1)
    # }
    
    if (l > Burnin){
      #1st stage parameters
      zvec.mat[,(l-Burnin)] = z.vec
      #main parameters
      theta.all = theta.all+as.matrix(theta)
      theta2.all = theta.all+as.matrix(theta^2)
      #estimated density
      d.biv.grid.tot = d.biv.grid.tot + d.biv.iter1
      
      d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
      
      d.prec.grid.tot = d.prec.grid.tot + d.prec.iter1
      
      #total.log.like is computed based on mud1 and precd1
      total.log.like.vec[(l-Burnin)] = total.log.like(mud1, precd1, m.vec, lambda.vec, alpha.vec, 
                                                      beta.vec, z.vec, pi.vec)
      #data.log.like is compared on mud and precd
      data.log.like.vec[(l-Burnin)] = data.log.like(mud, precd, m.vec*SS1+MM1, lambda.vec*(SS2/SS1)^2, 
                                                    alpha.vec, beta.vec*SS2^2, z.vec, pi.vec)
      k.vec[(l-Burnin)] = length(unique(z.vec))
    }
    
  }
  
  
  d.biv.grid.est = d.biv.grid.tot/Simsize
  d.mu.grid.est = d.mu.grid.tot/Simsize
  d.prec.grid.est = d.prec.grid.tot/Simsize
  
  #Because of level switching these estimatea are not reliable.
  theta.hat = data.frame(theta.all/Simsize)
  theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
  
  names(theta.hat) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
  names(theta.hat.var) = c("m.vec.est","lambda.vec.est","alpha.vec.est","beta.vec.est","pi.vec.est")
  
  order.pi.vec = order(c(-theta.hat$pi.vec))
  theta.hat = theta.hat[order.pi.vec,]
  theta.hat.var = theta.hat.var[order.pi.vec,]
  
  #Clustering analysis
  find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
  find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans2=max(tab)/length(x); rm(tab); return(ans2)}
  
  z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
  z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
  
  theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  order.pi.vec.eff = order(-theta.eff$pi.vec.est)
  theta.eff = theta.eff[order.pi.vec.eff,]
  theta.eff.var = theta.eff.var[order.pi.vec.eff,]
  
  return(list(z.vec.est=z.vec.est, 
              z.vec.est.prop=z.vec.est.prop, 
              theta.hat=theta.hat, theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
              d.biv.grid.est=d.biv.grid.est,
              d.mu.grid.est=d.mu.grid.est, d.prec.grid.est=d.prec.grid.est, 
              total.log.like.vec=total.log.like.vec, 
              data.log.like.vec=data.log.like.vec, k.vec=k.vec))
  
}

# Output: z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 


#**********************************************************************************

# \sigma_i^2 follows mixture of inverse-gamma. \sigma_i^2 is observed.

#**********************************************************************************

# input: sigma2d = vector of length q where ith element is \sigma_i^2 which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.ig.sknown = function(sigma2d, gamma.pi = 0.1, k = 10, Burnin = 5000, Simsize = 10000, 
                          sigma2.grid = c('Default'), hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("invgamma")) install.packages("invgamma")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(invgamma)
  # library(truncnorm)
  
  # rinvgamma = invgamma::rinvgamma
  # dinvgamma = invgamma::dinvgamma
  
  #*********************************Data Management****************************************#
  
  q = length(sigma2d)
  
  #centering and scaling the whole dataset
  SS2 = sd(sigma2d)
  # SS2 = 1
  sigma2d1 = sigma2d/SS2
  
  #******************************End of Data Management************************************#
  
  #******************************Initialize MCMC parameters************************************#
  
  #k-means clusters to determine initial values of cluster parameters.
  mscluster = kmeans(sigma2d,k)
  
  #latent variable which indicates cluster number. zi=r means sigma.i^2|z.i=r \sim IG(alpha.r, beta.r).
  z.vec = mscluster$cluster
  
  #vector of alpha.r's, beta.r's where r=1,...,k. alphavec, betavec, pivec is vector of 
  #alpha.r's, beta.r's and pi.r's with length k.
  alpha.vec = beta.vec = mh.ratio1.vec = mh.ratio2.vec = vector()
  
  #matrix to store Z.i after each iteration.
  zvec.mat = array(0,c(q,Simsize))
  
  #vector of length of unique z.i's in each iteration.
  k.vec = vector()
  
  
  #matrix of m.r's, lambda.r's, alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
  theta.all = array(0,c(k,3))
  theta2.all = array(0,c(k,3))
  
  #initial values alphavec, betavec, lambdavec.
  alpha_0 = (mean(sigma2d1))^2/var(sigma2d1)+2
  beta_0 = ((mean(sigma2d1))^2/var(sigma2d1)+1)*mean(sigma2d1)
  alpha.vec = rep(alpha_0,k)
  beta.vec = rep(beta_0,k)
  pi.vec = rep(1/k,k)
  
  if(sum(hyperparameters==c('Default'))>0){
    #alphavec follows Gamma(a.alpha,b.alpha). 
    a.alpha = 1
    b.alpha = 1
    
    #betavec follows Gamma(a.beta,b.beta). 
    a.beta = 1
    b.beta = 1/(sd(sigma2d)^2)*SS2^2
  }
  
  if (sum(hyperparameters!=c('Default'))>0){
    
    a.alpha = hyperparameters[1]
    b.alpha = hyperparameters[2]
    a.beta = hyperparameters[3]
    b.beta = hyperparameters[4]*SS2^2
  }
  
  #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
  #Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
  #allows higher variation in pi.r's. gamma.pi = 0.1
  
  if (sum(sigma2.grid==c('Default'))>0){
    sigma2.grid = seq(max(min(sigma2d)-IQR(sigma2d),0.001), max(sigma2d)+IQR(sigma2d),
                      length.out=2*100+1)[seq(2,2*100,2)]
  }
  
  #sigma2.grid is a univariate grid on which we are estimating the density.
  n2 = length(sigma2.grid)
  
  #array to store density on a grid
  d.sigma2.grid.tot = rep(0,n2)
  total.log.like.vec = data.log.like.vec = vector()
  
  d.sigma2.mixture = function(sigma2.grid, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dinvgamma(sigma2.grid, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  loglike.alpha = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.alpha, b.alpha, r){
    ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
    return(ans)
  }
  
  loglike.beta = function(alpha.vec, beta.vec, sigma2.vec, z.vec, a.beta, b.beta,r){
    ans = sum(dinvgamma(sigma2.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
    return(ans)
  }
  
  # parameters should be in rescaled scale. As it is on sigma2d1
  total.log.like = function(sigma2.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(count.vec*log(pi.vec))+
      sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
      sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
      lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
    
    return(ans)
  }
  
  
  # parameters should be in original scale. As it is on sigma2d
  data.log.like = function(sigma2.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dinvgamma(sigma2.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(count.vec*log(pi.vec))
    return(ans)
  }
  
  
  #*******************************MCMC algorithm*******************************************#
  
  for (l in 1:(Burnin+Simsize)){
    #updating z.is
    logig.mat = array(0,dim=c(q,k))
    for (r in 1:k){
      logig.mat[,r] = dinvgamma(sigma2d1,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
    }
    
    logig1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + logig.mat
    
    #need it for computational purpose. Otherwise exp(small log values) going to be 0.
    logig2.mat = logig1.mat - apply(logig1.mat,1,max)
    
    ig.mat = exp(logig2.mat)
    
    cum.ig.mat = ig.mat %*% upper.tri(diag(k), diag = TRUE)
    
    random.q.vec = runif(q)*rowSums(ig.mat)
    
    z.vec = rowSums(random.q.vec > cum.ig.mat)+1
    
    #countvec is a vector of length k where each element is number of zis in each cluster
    count.vec = tabulate(z.vec,nbins=k)
    
    #updating pi.rs
    pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
    
    #updating alpha.rs
    alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
    
    for (r in 1:k){
      if(count.vec[r] > 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = loglike.alpha(alphacand.vec, beta.vec, sigma2d1, z.vec, a.alpha, b.alpha, r)-
          loglike.alpha(alpha.vec, beta.vec, sigma2d1, z.vec, a.alpha, b.alpha, r)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
          dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
    }
    
    #updating beta.rs
    betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS2^2)
    
    for (r in 1:k){
      if(count.vec[r] > 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = loglike.beta(alpha.vec, betacand.vec, sigma2d1, z.vec, a.beta, b.beta, r)-
          loglike.beta(alpha.vec, beta.vec, sigma2d1, z.vec, a.beta, b.beta, r)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
          dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
    }
    
    d.sigma2.iter1 = d.sigma2.mixture(sigma2.grid, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    theta = data.frame(cbind(alpha.vec,beta.vec*SS2^2,pi.vec))
    names(theta) = c("alpha.vec","beta.vec","pi.vec")
    
    
    # if(l%%1000==0){
    #   theta1 = theta
    #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
    #   theta1 = theta1[order(c(-theta1$pi.vec)),]
    #   print(l)
    #   print(theta1)
    # }
    
    if (l > Burnin){
      #1st stage parameters
      zvec.mat[,(l-Burnin)] = z.vec
      #main parameters
      theta.all = theta.all+as.matrix(theta)
      theta2.all = theta.all+as.matrix(theta^2)
      #estimated density
      d.sigma2.grid.tot = d.sigma2.grid.tot + d.sigma2.iter1
      
      #total.log.like is computed based on mud1 and sigma2d1
      total.log.like.vec[(l-Burnin)] = total.log.like(sigma2d1, alpha.vec, 
                                                      beta.vec, z.vec, pi.vec)
      #data.log.like is compared on mud and sigma2d
      data.log.like.vec[(l-Burnin)] = data.log.like(sigma2d, 
                                                    alpha.vec, beta.vec*SS2^2, z.vec, pi.vec)
      k.vec[(l-Burnin)] = length(unique(z.vec))
    }
    
  }
  
  d.sigma2.grid.est = d.sigma2.grid.tot/Simsize
  
  #Because of level switching these estimatea are not reliable.
  theta.hat = data.frame(theta.all/Simsize)
  theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
  
  names(theta.hat) = c("alpha.vec.est","beta.vec.est","pi.vec.est")
  names(theta.hat.var) = c("alpha.vec.est","beta.vec.est","pi.vec.est")
  
  order.pi.vec = order(c(-theta.hat$pi.vec))
  theta.hat = theta.hat[order.pi.vec,]
  theta.hat.var = theta.hat.var[order.pi.vec,]
  
  #Clustering analysis
  find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
  find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans2=max(tab)/length(x); rm(tab); return(ans2)}
  
  z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
  z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
  
  theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  order.pi.vec.eff = order(-theta.eff$pi.vec.est)
  theta.eff = theta.eff[order.pi.vec.eff,]
  theta.eff.var = theta.eff.var[order.pi.vec.eff,]
  
  return(list(z.vec.est=z.vec.est, 
              z.vec.est.prop=z.vec.est.prop, 
              theta.hat=theta.hat, theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
              d.sigma2.grid.est=d.sigma2.grid.est, 
              total.log.like.vec=total.log.like.vec, 
              data.log.like.vec=data.log.like.vec, k.vec=k.vec))
  
}

# Output: z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         d.sigma2.grid.est = estimated sigma2 density on sigma2.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 

#**********************************************************************************

# \mu_i follows mixture of normal. \mu_i is observed.

#**********************************************************************************

# input: mud = vector of length q where ith element is \mu_i which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.nor.mknown = function(mud, gamma.pi = 0.1, k = 10, Burnin = 5000, Simsize = 10000, 
                           mu.grid = c('Default'), hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(truncnorm)
  
  
  #*********************************Data Management****************************************#
  
  q = length(mud)
  
  #centering and scaling the whole dataset
  MM1 = mean(mud)
  SS1 = sd(mud)
  # MM1 = 0
  # SS1 = 1
  # SS2 = 1
  mud1 = (mud-MM1)/SS1
  
  #******************************End of Data Management************************************#
  
  #******************************Initialize MCMC parameters************************************#
  
  #k-means clusters to determine initial values of cluster parameters.
  mscluster = kmeans(mud,k)
  
  #vector of m.r's where r=1,...,k. mvec is vector of m.r's with length k.
  #initial value of mvec
  m.vec = (mscluster$center-MM1)/SS1
  
  #latent variable which indicates cluster number. zi=r means sigma.i^2|z.i=r \sim N(m.r, lambda.r)
  z.vec = mscluster$cluster
  
  lambda.vec = mscluster$withinss/(SS1^2*tabulate(z.vec,nbins = k))
  
  
  #lambda.r's where r=1,...,k. lambdavec is vector of lambda.r's with length k.
  mh.ratio1.vec = mh.ratio2.vec = vector()
  
  #matrix to store Z.i after each iteration.
  zvec.mat = array(0,c(q,Simsize))
  
  #vector of length of unique z.i's in each iteration.
  k.vec = vector()
  
  
  #matrix of m.r's, lambda.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
  theta.all = array(0,c(k,3))
  theta2.all = array(0,c(k,3))
  
  pi.vec = rep(1/k,k)
  
  if(sum(hyperparameters==c('Default'))>0){
    #Hyperprior parameters. mvec follows N(m.0,\tau^2). m.0=E(X.ij), and tau2 < Var(\bar{X_i}). 
    #So, setting tau2 = Var(\bar{X_i}) less informative than empirical bayes.
    m.0 = (mean(mud)-MM1)/SS1
    tau2 = var(mud)/SS1^2
    
    #lambdavec follows InvGamma(a.lambda=1,b.lambda=1).
    a.lambda = 1
    b.lambda = 1
  }
  
  if (sum(hyperparameters!=c('Default'))>0){
    
    m.0 = (hyperparameters[1]-MM1)/SS1
    tau2 = hyperparameters[2]/SS1^2
    a.lambda = hyperparameters[3]
    b.lambda = hyperparameters[4]/SS1^2
  }
  
  #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
  #Dirichlet(gamma.pi/k,...,gamma.pi/k). 
  #Low value of gamma.pi non-informative which allows higher variation in pi.r's. gamma.pi = 0.1
  if (sum(mu.grid==c('Default'))>0){
    mu.grid = seq(min(mud)-IQR(mud), max(mud)+IQR(mud),
                  length.out=2*100+1)[seq(2,2*100,2)]
  }
  
  #mu.grid is a univariate grid on which we are estimating the density.
  n1 = length(mu.grid)
  
  # array to store
  d.mu.grid.tot = rep(0,n1)
  total.log.like.vec = data.log.like.vec = vector()
  
  
  
  d.mu.mixture = function(mu.grid, m.vec, lambda.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]* dnorm(mu.grid, mean=m.vec, sd=sqrt(lambda.vec))
    }
    return(ans)
  }
  
  
  #parameters should be in rescaled scale. As it is on mud1
  total.log.like = function(mu.vec, m.vec, lambda.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    
    
    ans = sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(lambda.vec.large), log=TRUE))+
      sum(count.vec*log(pi.vec))+
      sum(dnorm(m.vec, mean=m.0, sd=sqrt(tau2)))+
      sum(dinvgamma(lambda.vec, shape=a.lambda, rate=b.lambda, log=TRUE))+
      lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
    
    return(ans)
  }
  
  
  #parameters should be in original scale. As it is on mud
  data.log.like = function(mu.vec, m.vec, lambda.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    m.vec.large = m.vec[z.vec]
    lambda.vec.large = lambda.vec[z.vec]
    
    ans = sum(dnorm(mu.vec, mean=m.vec.large, sd=sqrt(lambda.vec.large), log=TRUE))+
      sum(count.vec*log(pi.vec))
    return(ans)
  }
  
  #*******************************MCMC algorithm*******************************************#
  
  for (l in 1:(Burnin+Simsize)){
    
    #updating z.is
    lognor.mat = array(0,dim=c(q,k))
    for (r in 1:k){
      lognor.mat[,r] = dnorm(mud1,mean=m.vec[r],sd=sqrt(lambda.vec[r]),log=TRUE)
    }
    
    lognor1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + lognor.mat
    
    #need it for computational purpose. Otherwise exp(small log values) going to be 0.
    lognor2.mat = lognor1.mat - apply(lognor1.mat,1,max)
    
    nor.mat = exp(lognor2.mat)
    
    cum.nor.mat = nor.mat %*% upper.tri(diag(k), diag = TRUE)
    
    random.q.vec = runif(q)*rowSums(nor.mat)
    
    z.vec = rowSums(random.q.vec > cum.nor.mat)+1
    
    #countvec is a vector of length k where each element is number of zis in each cluster
    count.vec = tabulate(z.vec,nbins=k)
    
    #updating pi.rs
    pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
    
    #updating m.ts
    for (r in 1:k){
      if(count.vec[r] > 0){
        m.vec[r] = rnorm(1, mean=(sum(mud1[which(z.vec==r)])/lambda.vec[r]+m.0/tau2)/
                           (count.vec[r]/lambda.vec[r]+1/tau2),
                         sd=sqrt(1/(count.vec[r]/lambda.vec[r]+1/tau2)))
      }
      if(count.vec[r] == 0){
        m.vec[r] = rnorm(1, mean=m.0, sd=sqrt(tau2))
      }
    }
    
    #updating lambda.rs
    for (r in 1:k){
      if(count.vec[r] > 0){
        lambda.vec[r] = rinvgamma(1, shape=count.vec[r]/2+a.lambda, 
                                  rate=sum((mud1[which(z.vec==r)]-m.vec[r])^2)/2+b.lambda)
      }
      if(count.vec[r] == 0){
        lambda.vec[r] = rinvgamma(1, shape=a.lambda, rate=b.lambda)
      }
    }
    
    
    d.mu.iter1 = d.mu.mixture(mu.grid, m.vec*SS1+MM1, lambda.vec*SS1^2, pi.vec)
    
    theta = data.frame(cbind(m.vec*SS1+MM1, lambda.vec*SS1^2, pi.vec))
    names(theta) = c("m.vec","lambda.vec","pi.vec")
    
    # if(l%%1000==0){
    #   theta1 = theta
    #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
    #   theta1 = theta1[order(c(-theta1$pi.vec)),]
    #   print(l)
    #   print(theta1)
    # }
    
    if (l > Burnin){
      #1st stage parameters
      zvec.mat[,(l-Burnin)] = z.vec
      #main parameters
      theta.all = theta.all+as.matrix(theta)
      theta2.all = theta.all+as.matrix(theta^2)
      
      #estimated density
      d.mu.grid.tot = d.mu.grid.tot + d.mu.iter1
      
      #total.log.like is computed based on mud1
      total.log.like.vec[(l-Burnin)] = total.log.like(mud1, m.vec, lambda.vec, z.vec, pi.vec)
      #data.log.like is compared on mud
      data.log.like.vec[(l-Burnin)] = data.log.like(mud, m.vec*SS1+MM1, lambda.vec*SS1^2, z.vec, pi.vec)
      k.vec[(l-Burnin)] = length(unique(z.vec))
    }
    
  }
  
  d.mu.grid.est = d.mu.grid.tot/Simsize
  
  #Because of level switching these estimatea are not reliable.
  theta.hat = data.frame(theta.all/Simsize)
  theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
  
  names(theta.hat) = c("m.vec.est","lambda.vec.est","pi.vec.est")
  names(theta.hat.var) = c("m.vec.est","lambda.vec.est","pi.vec.est")
  
  order.pi.vec = order(c(-theta.hat$pi.vec))
  theta.hat = theta.hat[order.pi.vec,]
  theta.hat.var = theta.hat.var[order.pi.vec,]
  
  #Clustering analysis
  find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
  find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans2=max(tab)/length(x); rm(tab); return(ans2)}
  
  z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
  z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
  
  theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  order.pi.vec.eff = order(-theta.eff$pi.vec.est)
  theta.eff = theta.eff[order.pi.vec.eff,]
  theta.eff.var = theta.eff.var[order.pi.vec.eff,]
  
  return(list(z.vec.est=z.vec.est, 
              z.vec.est.prop=z.vec.est.prop, 
              theta.hat=theta.hat, theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
              d.mu.grid.est=d.mu.grid.est, 
              total.log.like.vec=total.log.like.vec, 
              data.log.like.vec=data.log.like.vec, k.vec=k.vec))
  
}

# Output: z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,sigma_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.biv.grid.est = estimated bivariate density of the 2-D grid defined by mu.grid and sigma2.grid.
#         d.mu.grid.est = estimated mu density on mu.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 

#**********************************************************************************

# precision = 1/\sigma_i^2 follows mixture of gamma. \sigma_i^2 is observed.

#**********************************************************************************

# input: precd = vector of length q where ith element is prec.i = 1/\sigma_i^2 which is known.
#        gamma.pi = Dirichlet process parameter.
#        k = number of maxumum cluster. A paramter for truncated dirichlet process.
#        Burnin = number of initial MCMC samples which will be discarded.
#        Simsize = number of MCMC samples that will be used to estimate parameters.
#        mu.grid = grid points over which the density of mu will be estimated.
#        sigma2.grid = grid points over which the density of sigma^2 will be estimated.
#        hyperparameters = vector of hyperparameters (m.0, tau2, a.lambda, b.lambda, a.alpha, b.alpha)

DPMM.g.sknown = function(precd, gamma.pi = 0.1, k = 10, Burnin = 5000, Simsize = 10000, 
                         prec.grid = c('Default'), hyperparameters = c('Default')){
  
  # if (!require("MCMCpack")) install.packages("MCMCpack")
  # if (!require("truncnorm")) install.packages("truncnorm")
  
  # library(MCMCpack)
  # library(truncnorm)
  
  #*********************************Data Management****************************************#
  
  q = length(precd)
  
  #centering and scaling the whole dataset
  SS2 = sd(precd)
  # SS1 = 1
  # SS2 = 1
  precd1 = precd/SS2
  
  #******************************End of Data Management************************************#
  
  #******************************Initialize MCMC parameters************************************#
  
  #k-means clusters to determine initial values of cluster parameters.
  mscluster = kmeans(precd,k)
  
  #latent variable which indicates cluster number. zi=r means prec.i|z.i=r \sim Gamma(alpha.r, beta.r)
  z.vec = mscluster$cluster
  
  #vector of M-H ratios of length k.
  alpha.vec = beta.vec = mh.ratio1.vec = mh.ratio2.vec = vector()
  
  #matrix to store Z.i after each iteration.
  zvec.mat = array(0,c(q,Simsize))
  
  #vector of length of unique z.i's in each iteration.
  k.vec = vector()
  
  #matrix of alpha.r's, beta.r's, pi.r's  where r=1,...,k. theta is the matrix of paramaters.
  theta.all = array(0,c(k,3))
  theta2.all = array(0,c(k,3))
  
  #initial values alphavec, betavec.
  alpha_0 = (mean(precd1))^2/var(precd1)
  beta_0 = mean(precd1)/var(precd1)
  alpha.vec = rep(alpha_0,k)
  beta.vec = rep(beta_0,k)
  pi.vec = rep(1/k,k)
  
  if(sum(hyperparameters==c('Default'))>0){
    #alphavec follows Gamma(a.alpha,b.alpha). 
    a.alpha = 1
    b.alpha = 1
    
    #betavec follows Gamma(a.beta,b.beta). 
    a.beta = 1
    b.beta = 1/(sd(precd)^2)*SS2^2
  }
  
  if (sum(hyperparameters!=c('Default'))>0){
    a.alpha = hyperparameters[1]
    b.alpha = hyperparameters[2]
    a.beta = hyperparameters[3]
    b.beta = hyperparameters[4]*SS2^2
  }
  
  #gamma.pi is dirichlet process prior on pi.r. pi.r=Prob(z.i=r). pivec vector of length k follows 
  #Dirichlet(gamma.pi/k,...,gamma.pi/k). Low value of gamma.pi non-informative.
  #allows higher variation in pi.r's. gamma.pi = 0.1
  
  if (sum(prec.grid==c('Default'))>0){
    prec.grid = seq(max(min(precd)-IQR(precd),0.001), max(precd)+IQR(precd),
                    length.out=2*100+1)[seq(2,2*100,2)]
  }
  
  #prec.grid is a univariate grid on which we are estimating the density.
  n2 = length(prec.grid)
  
  # array to store density on a grid
  d.prec.grid.tot = rep(0,n2)
  
  total.log.like.vec = data.log.like.vec = vector()
  
  d.prec.mixture = function(prec.grid, alpha.vec, beta.vec, pi.vec){
    ans = 0
    k.comp = length(pi.vec)
    for (r in 1:k.comp){
      ans = ans + pi.vec[r]*dgamma(prec.grid, shape=alpha.vec[r], rate=beta.vec[r])
    }
    return(ans)
  }
  
  loglike.alpha = function(alpha.vec, beta.vec, prec.vec, z.vec, a.alpha, b.alpha, r){
    ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
    return(ans)
  }
  
  loglike.beta = function(alpha.vec, beta.vec, prec.vec, z.vec, a.beta, b.beta,r){
    ans = sum(dgamma(prec.vec[which(z.vec==r)], shape=alpha.vec[r], rate=beta.vec[r], log=TRUE))+
      dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
    return(ans)
  }
  
  #parameters should be in rescaled scale. As it is on precd1
  total.log.like = function(prec.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(count.vec*log(pi.vec))+
      sum(dgamma(alpha.vec, shape=a.alpha, rate=b.alpha, log=TRUE))+
      sum(dgamma(beta.vec, shape=a.beta, rate=b.beta, log=TRUE))+
      lgamma(gamma.pi)-k*lgamma(gamma.pi/k)+(gamma.pi/k-1)*sum(log(pi.vec))
    
    return(ans)
  }
  
  
  #parameters should be in original scale. As it is on precd
  data.log.like = function(prec.vec, alpha.vec, beta.vec, z.vec, pi.vec){
    count.vec = tabulate(z.vec,nbins=k)
    alpha.vec.large = alpha.vec[z.vec]
    beta.vec.large = beta.vec[z.vec]
    
    ans = sum(dgamma(prec.vec, shape=alpha.vec.large, rate=beta.vec.large, log=TRUE))+
      sum(count.vec*log(pi.vec))
    return(ans)
  }
  
  #*******************************MCMC algorithm*******************************************#
  
  for (l in 1:(Burnin+Simsize)){
    
    #updating z.is
    logng.mat = array(0,dim=c(q,k))
    for (r in 1:k){
      logng.mat[,r] = dgamma(precd1,shape=alpha.vec[r],rate=beta.vec[r],log=TRUE)
    }
    
    logng1.mat = matrix(rep(log(pi.vec),each=q),nrow=q) + logng.mat
    
    #need it for computational purpose. Otherwise exp(small log values) going to be 0.
    logng2.mat = logng1.mat - apply(logng1.mat,1,max)
    
    ng.mat = exp(logng2.mat)
    
    cum.ng.mat = ng.mat %*% upper.tri(diag(k), diag = TRUE)
    
    random.q.vec = runif(q)*rowSums(ng.mat)
    
    z.vec = rowSums(random.q.vec>cum.ng.mat)+1
    
    #countvec is a vector of length k where each element is number of zis in each cluster
    count.vec = tabulate(z.vec,nbins=k)
    
    #updating pi.rs
    pi.vec = c(rdirichlet(1,count.vec+gamma.pi/k))
    
    #updating alpha.rs
    alphacand.vec = rnorm(k,mean=alpha.vec,sd=1)
    
    for (r in 1:k){
      if(count.vec[r] > 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = loglike.alpha(alphacand.vec, beta.vec, precd1, z.vec, a.alpha, b.alpha, r)-
          loglike.alpha(alpha.vec, beta.vec, precd1, z.vec, a.alpha, b.alpha, r)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & alphacand.vec[r] > 0){
        mh.ratio1.vec[r] = dgamma(alphacand.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)-
          dgamma(alpha.vec[r], shape=a.alpha, rate=b.alpha, log=TRUE)
        if(mh.ratio1.vec[r] > log(runif(1))){
          alpha.vec[r] = alphacand.vec[r]
        }
      }
    }
    
    #updating beta.rs
    betacand.vec = rnorm(k,mean=beta.vec,sd=4/SS2^2)
    
    for (r in 1:k){
      if(length(which(z.vec==r)) > 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = loglike.beta(alpha.vec, betacand.vec, precd1,z.vec, a.beta, b.beta, r)-
          loglike.beta(alpha.vec, beta.vec, precd1, z.vec, a.beta, b.beta, r)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
      if(count.vec[r] == 0 & betacand.vec[r] > 0){
        mh.ratio2.vec[r] = dgamma(betacand.vec[r], shape=a.beta, rate=b.beta, log=TRUE)-
          dgamma(beta.vec[r], shape=a.beta, rate=b.beta, log=TRUE)
        if(mh.ratio2.vec[r] > log(runif(1))){
          beta.vec[r] = betacand.vec[r]
        }
      }
    }
    
    d.prec.iter1 = d.prec.mixture(prec.grid, alpha.vec, beta.vec*SS2^2, pi.vec)
    
    theta = data.frame(cbind(alpha.vec, beta.vec*SS2^2, pi.vec))
    names(theta) = c("alpha.vec","beta.vec","pi.vec")
    
    # if(l%%1000==0){
    #   theta1 = theta
    #   theta1 = theta1[as.vector(data.frame(table(z.vec))$z.vec),]
    #   theta1 = theta1[order(c(-theta1$pi.vec)),]
    #   print(l)
    #   print(theta1)
    # }
    
    if (l > Burnin){
      #1st stage parameters
      zvec.mat[,(l-Burnin)] = z.vec
      #main parameters
      theta.all = theta.all+as.matrix(theta)
      theta2.all = theta.all+as.matrix(theta^2)
      #estimated density
      d.prec.grid.tot = d.prec.grid.tot + d.prec.iter1
      
      #total.log.like is computed based on mud1 and precd1
      total.log.like.vec[(l-Burnin)] = total.log.like(precd1, alpha.vec, 
                                                      beta.vec, z.vec, pi.vec)
      #data.log.like is compared on mud and precd
      data.log.like.vec[(l-Burnin)] = data.log.like(precd,
                                                    alpha.vec, beta.vec*SS2^2, z.vec, pi.vec)
      k.vec[(l-Burnin)] = length(unique(z.vec))
    }
    
  }
  
  d.prec.grid.est = d.prec.grid.tot/Simsize
  
  #Because of level switching these estimatea are not reliable.
  theta.hat = data.frame(theta.all/Simsize)
  theta.hat.var = data.frame(theta2.all/Simsize-(theta.all/Simsize)^2)
  
  names(theta.hat) = c("alpha.vec.est","beta.vec.est","pi.vec.est")
  names(theta.hat.var) = c("alpha.vec.est","beta.vec.est","pi.vec.est")
  
  order.pi.vec = order(c(-theta.hat$pi.vec))
  theta.hat = theta.hat[order.pi.vec,]
  theta.hat.var = theta.hat.var[order.pi.vec,]
  
  #Clustering analysis
  find.discrete.mode = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans1=ux[tab==max(tab)][1]; rm(tab); return(ans1)}
  find.discrete.mode.freq = function(x){ux=unique(x); tab=tabulate(match(x,ux)); 
  ans2=max(tab)/length(x); rm(tab); return(ans2)}
  
  z.vec.est = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode))
  z.vec.est.prop = unlist(apply(zvec.mat[,(Simsize-min(1000,Simsize/2)):Simsize],1,find.discrete.mode.freq))
  
  theta.eff = theta.hat[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  theta.eff.var = theta.hat.var[as.numeric(data.frame(table(z.vec.est))$z.vec.est),]
  order.pi.vec.eff = order(-theta.eff$pi.vec.est)
  theta.eff = theta.eff[order.pi.vec.eff,]
  theta.eff.var = theta.eff.var[order.pi.vec.eff,]
  
  return(list(z.vec.est=z.vec.est, 
              z.vec.est.prop=z.vec.est.prop, 
              theta.hat=theta.hat, theta.hat.var=theta.hat.var, theta.eff=theta.eff, theta.eff.var=theta.eff.var,
              d.prec.grid.est=d.prec.grid.est, 
              total.log.like.vec=total.log.like.vec, 
              data.log.like.vec=data.log.like.vec, k.vec=k.vec))
  
}

# Output: z.vec.est = an indicator vector to indicates the corresponding cluster of (mu_i,prec_i^2).
#         z.vec.est.prop = indicates reliabilty of z.vec.est. Proportions of MCMC iterations agree with z.vec.est.
#         theta.hat = estimated of 1st level hyperparameters. Not very reliable because of level switching.
#         theta.hat.var = variance of theta.hat.
#         theta.eff = part of theta.hat. Only active components of theta.hat.
#         theta.eff.var = part of theta.hat.var. Only active components of theta.hat.var.
#         d.prec.grid.est = estimated prec density on prec.grid.
#         total.log.like.vec = total loglikelihood as defined in a function, calculated for each iteration.
#         data.log.like.vec = data loglikelihood as defined in a function, calculated for each iteration.
#         k.vec = number active clusters of each iteration. 


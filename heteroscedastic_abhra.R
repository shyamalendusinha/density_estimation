  ####################################################################################
   ### 	R program for fitting Bayesian semiparametric density deconvolution models  ### 
   ###  This program assumes that the measurement errors are heteroscedastic        ###
   ####################################################################################   
   
   
   	#################
   	### Author(s) ###
   	#################
   
# Author of the codes: Abhra Sarkar
# Authors of the paper: Abhra Sarkar, Bani K. Mallick, John Staudenmayer, Debdeep Pati, Raymond J. Carroll
# Author for correspondence: Raymond J. Carroll (carroll@stat.tamu.edu)

UNIVARIATE_DECON_HETEROSCEDASTIC <- function(ws,mis,error.fit=c("normal"),simsize=5000,burnin=3000,plot_results=TRUE,x.grid = x.grid)
{

  variable.name <- c("varname")

  ######################################
  ### Add Libraries and Source Files ###
  ######################################
  
  # install.packages("mvtnorm")
  # install.packages("msm")
  # install.packages("MCMCpack")
  library(mvtnorm) 
  library(msm)
  library(MCMCpack)
  
  n = length(mis)
  
  inds <- rep(1:n,times=mis)
  
##############################Functions###################################################

#############################################################

# Function to create Double difference matrix which we needed to penalize 2nd deference of B-spline parameters for smoothing. (Eilers and Marx 1996)

P.mat <- function(K)
{
  # penalty matrix for density
  D <- diag(rep(1,K))
  D <- diff(diff(D))
  P <- t(D)%*%D 
  return(P)
}

#P.mat(4)

#############################################################

# Defining B-spline basis function using recursion defined in De Boor page 90. 

B.basis <- function(x,knots)
{
  # delta <- knots[2]-knots[1]
  n <- length(x)
  K <- length(knots)
  B <- matrix(0,n,K+1)
  for (jj in 1:(K-1))
  {
    act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
    act.x <- x[act.inds]
    resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
    
    B[act.inds,jj] <- (1/2)*(1-resc.x)^2
    B[act.inds,jj+1] <- -(resc.x^2)+resc.x+1/2
    B[act.inds,jj+2] <- (resc.x^2)/2
  }
  return(B)
}

knots = seq(0,1,length.out=10)
x = sample(seq(0,1,by=0.001), 100, replace = FALSE, prob = NULL)
round(B.basis(x,knots),3)

#############################################################

# Not used in this paper. This function does normalize B-spline Basis function so that it can be used in B-spline density estimation. Staudenmayer(2008)

B.basis.normalized <- function(x,knots)		# equidistant knots
{
  B <- B.basis(x,knots)
  delta <- knots[2]-knots[1]
  K <- length(knots)
  area.tot <- (K-1)*delta
  return(B/area.tot)
}

knots = seq(0,1,length.out=10)
x = sample(seq(0,1,by=0.001), 100, replace = FALSE, prob = NULL)
round(B.basis.normalized(x,knots),3)

#############################################################


rnormbspline <- function(n,coeffs,minx,maxx)
{
  grid <- seq(minx,maxx,len=500)
  delta <- (grid[2]-grid[1])
  knots <- seq(minx,maxx,len=(length(coeffs)-1))
  y <- B.basis.normalized(grid,knots)%*%(coeffs*length(knots))
  y <- y/(sum(y)*delta)
  Fy <- cumsum(y)*delta
  u <- runif(n)	
  x <- numeric(n)
  for(ii in 1:n)
    x[ii] <- grid[min(which(u[ii]<Fy))+1]
  x <- x + runif(n,0,delta)
  return(x)
}


knots = seq(0,1,length.out=10)
x = sample(seq(0,1,by=0.001), 100, replace = FALSE, prob = NULL)
round(rnormbspline(100,c(0.5,0.75,0.8),0,1),3)

#############################################################


dnormbspline <- function(x,coeffs,minx,maxx)	# x must be equidistant grid
{
  knots <- seq(minx,maxx,len=(length(coeffs)-1))
  y <- B.basis.normalized(x,knots)%*%(coeffs*length(knots))
  delta <- (x[2]-x[1])
  y <- y/(sum(y)*delta)
  return(y)
}


#############################################################

# Random number generation for skew normal

rskewnorm <- function(n,mean,sd,skewness)
{
  c = 2/3.1415926
  delta = skewness/sqrt(1+skewness^2)
  xi = - (sqrt(c) * skewness)/sqrt(1 + skewness^2 * (1-c))
  omega = sqrt(1+xi^2)
  
  y = delta * abs(rnorm(n,0,1)) + (1-delta^2)^(1/2) * rnorm(n,0,1)
  y = xi + omega * y
  return(mean+sd*y)
}

# Density for skew normal

dskewnorm <- function(x,mean,sd,skewness)
{
  c = 2/3.1415926
  delta = skewness/sqrt(1+skewness^2)
  zeta1 = delta * sqrt(c)
  zeta2 = sqrt(1 - c*delta^2)
  y = numeric(length(x))
  xmod = zeta1 + zeta2*(x-mean)/sd
  for(i in 1:length(x))
    y[i] = (2*zeta2/sd[i]) * dnorm(xmod[i]) * pnorm(skewness*xmod[i])
  return(y)
}


#############################################################


d.scaled.restricted.mix.norm <- function(x,mean,sd,pi,params)
{
  if(is.matrix(params))
  {
    p = params[,1]
    mu_curl = params[,2]
    sigmasq1 = params[,3]
    sigmasq2 = params[,4]
  }
  else
  {
    p = params[1]
    mu_curl = params[2]
    sigmasq1 = params[3]
    sigmasq2 = params[4]
  }
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  mu1 = c1*mu_curl
  mu2 = c2*mu_curl
  
  sd.e = sqrt(var.e.fn(pi,params))
  y = matrix(0,nrow=length(p),ncol=length(x))
  for(kk in 1:length(p))
    y[kk,] = pi[kk]*(p[kk]*dnorm(x,(mean+sd*mu1[kk])/sd.e,sd*sqrt(sigmasq1[kk])/sd.e) + (1-p[kk])*dnorm(x,(mean+sd*mu2[kk])/sd.e,sd*sqrt(sigmasq2[kk])/sd.e))
  y = colSums(y)
  return(y)
}



d.restricted.mix.norm <- function(x,mean,sd,params)
{
  if(is.matrix(params))
  {
    p = params[,1]
    mu_curl = params[,2]
    sigmasq1 = params[,3]
    sigmasq2 = params[,4]
  }
  else
  {
    p = params[1]
    mu_curl = params[2]
    sigmasq1 = params[3]
    sigmasq2 = params[4]
  }
  
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  mu1 = c1*mu_curl
  mu2 = c2*mu_curl
  
  
  y = p*dnorm(x,mean+sd*mu1,sd*sqrt(sigmasq1)) + (1-p)*dnorm(x,mean+sd*mu2,sd*sqrt(sigmasq2))
  return(y)
}



mixnorm <- function(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,plot=TRUE)
{
  m = length(pi)
  y = numeric(n)
  density <- numeric(length(e.grid))
  inds = sample(1:m,n,TRUE,prob=pi)
  for(ii in 1:m)
  {
    temp = which(inds==ii)
    y[temp] = r.restricted.mix.norm(length(temp),c(p[ii],mu_curl[ii],sigmasq1[ii],sigmasq2[ii]))
  }
  params = cbind(p,mu_curl,sigmasq1,sigmasq2)
  y = y/sqrt(var.e.fn(pi,params))
  density <- d.scaled.restricted.mix.norm(e.grid,mean=0,sd=1,pi,params)
  if(plot)
  {
    hist(y,xlim=c(min(e.grid),max(e.grid)),breaks=30,freq=FALSE)
    points(e.grid,density,type="l")
  }
  
  return(list(es=y,density=density))
}


r.restricted.mix.norm <- function(n,params)
{
  p = params[1]
  mu_curl = params[2]
  sigmasq1 = params[3]
  sigmasq2 = params[4]
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  mu1 = c1*mu_curl
  mu2 = c2*mu_curl
  
  
  inds = sample(0:1,n,TRUE,prob=c(p,(1-p)))
  y = numeric(n)
  temp = which(inds==0)
  y[temp] = rnorm(length(temp),mu1,sqrt(sigmasq1)) 
  temp = which(inds==1)
  y[temp] = rnorm(length(temp),mu2,sqrt(sigmasq2))
  return(y)
}



var.e.fn <- function(pi,params)
{
  if(is.matrix(params))
  {
    p = params[,1]
    mu_curl = params[,2]
    sigmasq1 = params[,3]
    sigmasq2 = params[,4]
  }
  else
  {
    p = params[1]
    mu_curl = params[2]
    sigmasq1 = params[3]
    sigmasq2 = params[4]
  }
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  mu1 = c1*mu_curl
  mu2 = c2*mu_curl
  
  
  y = p*(mu1^{2}+sigmasq1) + (1-p)*(mu2^{2}+sigmasq2)
  y = sum(pi*y)
  return(y)
}



gamma.e.fn <- function(pi,params)
{
  if(is.matrix(params))
  {
    p = params[,1]
    mu_curl = params[,2]
    sigmasq1 = params[,3]
    sigmasq2 = params[,4]
  }
  else
  {
    p = params[1]
    mu_curl = params[2]
    sigmasq1 = params[3]
    sigmasq2 = params[4]
  }
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  mu1 = c1*mu_curl
  mu2 = c2*mu_curl
  
  
  m2 = p*(mu1^{2}+sigmasq1) + (1-p)*(mu2^{2}+sigmasq2)
  m2 = sum(pi*m2)
  
  m3 = p*(mu1^{3}+3*mu1*sigmasq1) + (1-p)*(mu2^{3}+3*mu2*sigmasq2)
  m3 = sum(pi*m3)
  
  m4 = p*(mu1^{4}+6*mu1^{2}*sigmasq1+3*sigmasq1^{2}) + (1-p)*(mu2^{3}+6*mu2^{2}*sigmasq2+3*sigmasq2^{2})
  m4 = sum(pi*m4)
  
  gamma1 = m3/m2^{3/2}
  gamma2 = m4/m2^{2}-3
  
  return(c(gamma1,gamma2))
}


r.proposal.params.restricted.mix.norm <- function(p.a,p.b,sigmasq.mu_curl,s11,s12,s21,s22)
{
  p = rbeta(1,p.a,p.b)
  
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  mu_curl = rnorm(1,0,sqrt(sigmasq.mu_curl))
  sigmasq1 = 1/rgamma(1,s11,s12)
  sigmasq2 = 1/rgamma(1,s21,s22)
  
  y = c(p,mu_curl,sigmasq1,sigmasq2)
  return(y)
}    


r.tnorm.proposal.params.restricted.mix.norm <- function(params)
{
  current.p <- params[1]
  current.mu_curl <- params[2]
  current.sigmasq1 <- params[3]
  current.sigmasq2 <- params[4]
  
  p = rtnorm(1,current.p,0.01,lower=0,upper=1)
  c1 = (1-p)/(p^2+(1-p)^2)
  c2 = -p/(p^2+(1-p)^2)
  
  
  sigmasq1 = rtnorm(1,current.sigmasq1,0.1,lower=0,upper=Inf)
  sigmasq2 = rtnorm(1,current.sigmasq2,0.1,lower=0,upper=Inf)
  
  mu_curl = rnorm(1,current.mu_curl,0.1)
  
  y = c(p,mu_curl,sigmasq1,sigmasq2)
  return(y)
}    




#############################################################

dlaplace <- function(x,mu,b)
{
  y = exp(-abs(x-mu)/b)/(2*b)
  return(y)
}



rlaplace <- function(n,mu,b)
{
  u = runif(n)
  y = mu - b*sign(u-0.5)*log(1-2*abs(u-0.5))
  return(y)
}



mixlaplace <- function(n,pi,mu,b,e.grid,plot=TRUE)		# Produces mixtures of Laplace
{
  m = length(pi)
  y = numeric(n)
  density = numeric(length(e.grid))
  inds = sample(1:m,n,TRUE,pi)
  mu.e = sum(pi*mu)
  sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
  for(ii in 1:m)
  {
    temp = which(inds==ii)
    y[temp] = rlaplace(length(temp),mu[ii],b[ii])
    density = density + pi[ii]*dlaplace(e.grid,(mu[ii]-mu.e)/sd.e,b[ii]/sd.e)
  }
  y = (y-mu.e)/sd.e
  if(plot)
  {
    hist(y,xlim=c(min(e.grid),max(e.grid)),freq=FALSE,breaks=30)
    points(e.grid,density,type="l")
  }
  return(list(es=y,density=density))
}


gamma.e.Laplace.fn <- function(pi,mu,b)
{
  m1.dash = sum(pi*mu) 
  m2.dash = sum(pi*(2*b^{2}+mu^{2}))
  m3.dash = sum(pi*(mu^{3}+6*b^{2}*mu))
  m4.dash = sum(pi*(mu^{4}+12*b^{2}*mu^{2}+24*b^{4}))
  
  m2 = m2.dash-m1.dash^{2}
  m3 = m3.dash - 3*m2.dash*m1.dash + 2*m1.dash^{3}
  m4 = m4.dash - 4*m3.dash*m1.dash + 6*m2.dash*m1.dash^{2} - 3*m1.dash^{4}
  
  gamma1 = m3/m2^{3/2}
  gamma2 = m4/m2^{2}-3
  
  #print(rbind(m1.dash,m2.dash,m3.dash,m4.dash))
  #print(rbind(m2,m3,m4))
  
  return(c(gamma1,gamma2))
}



#############################################################


fr <- function(thetas)			# function to be maximized
{
  vars = B.basis(xs,knots.t)%*%exp(thetas)
  y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum(rep(log(vars),times=mis))/2 - sum(us^2/rep(vars,times=mis))/2
  return(-y)
}


############################################################


gr <- function(thetas)			# gradient function of fr
{
  vars = B.basis(xs,knots.t)%*%exp(thetas)
  y = - P.t %*% thetas / s2t
  B.basis.components = matrix(0,nrow=K.t+1,ncol=n)
  for(kk in 1:(K.t+1))
  {
    thetas.new = rep(0,K.t+1)
    thetas.new[kk] = 1
    B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
  }
  for(kk in 1:(K.t+1))
    for(ii in 1:n)
      for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
        y[kk] = y[kk] - (1 - (us[jj]^2)/vars[ii]) * B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
  return(-y)
}


############################################################


prop.sig.thetas.fn <- function(thetas,s2t)
{
  vars = B.basis(xs,knots.t)%*%exp(thetas)
  prop.sig.thetas = matrix(0,nrow=K.t+1,ncol=K.t+1)
  B.basis.components = matrix(0,nrow=K.t+1,ncol=n);
  for(kk in 1:(K.t+1))
  {
    thetas.new = rep(0,K.t+1)
    thetas.new[kk] = 1
    B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
  }
  
  for(kk in 1:(K.t+1))
  {
    for(ll in kk:(K.t+1))
    {
      if(kk==ll)
        for(ii in 1:n)
          for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
            prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2) + (1-(us[jj]^2)/vars[ii])*B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
      else
        for(ii in 1:n)
          for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
            prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2)
    } 
  }
  for(kk in 2:(K.t+1))
    for(ll in 1:(kk-1))
      prop.sig.thetas[kk,ll] = prop.sig.thetas[ll,kk]
  prop.sig.thetas = prop.sig.thetas + P.t/s2t
  prop.sig.thetas = round(solve(prop.sig.thetas),4)
  
  return(prop.sig.thetas)
}


#############################################################


pi.fn <- function(x,alpha,beta,xstar,K)
{
  w <- pnorm(alpha-beta*abs(x-xstar)^2)
  pi <- numeric(K)
  pi[1] <- w[1]
  for(k in 2:(K-1))
    pi[k] <- w[k]*prod(1-w[1:(k-1)])
  pi[K] <- prod(w[1:(K-1)])
  return(pi)
}


##########################################################################################







   






	#################################
	### Priors and Initial Values ###
	#################################

### Initialization and prior of xs and us

wbars <- tapply(ws,inds,"mean") 
s2is <- as.vector(tapply(ws,inds,var))


xs <- current.xs <- start.xs <- as.vector(wbars)
us <- ws - rep(xs,times=mis)
# x.grid <- seq(min(wbars),max(wbars),length=500)



alpha.x=0.1	# robust but should not be set too small 

# Normal-Inverse-Scaled-Chisquare
nu0.x = 1/5
gama0.x = 3	                  			# Note the spelling, also note that this must be greater than 2
mu0.x = mean(xs)
sigmasq0.x = var(xs)/((1+1/nu0.x)*(gama0.x/(gama0.x-2)))
student.x.grid = (x.grid-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))

z.x <- 1:n								# one separate clusters for each obs
mu.x <- xs								# unique values
sigmasq.x <- rep(var(xs)/5,n)			# unique values



### Prior and initialization of s2t and thetas

alpha.t <- .1
beta.t <- .1

s2t <- current.s2t <- start.s2t <- 0.1

   K.t <- 10   			# num of knots, results are robust to variations in values
   P.t <- P.mat(K.t+1) 	# penalty matrix
   range.start.xs <- diff(range(start.xs))
   knots.t <- seq(min(start.xs)-.1*range.start.xs,
                  max(start.xs)+.1*range.start.xs,
                  length=(K.t))
   optim_results <- optim(rep(1,K.t+1), fr, gr, method = "BFGS")
   thetas <- current.thetas <- start.thetas <- optim_results$par
   sig.tune.thetas.1 <- 0		# results are very robust to these two tuning parameters
   sig.tune.thetas.2 <- 0.1	
   prop.sig.thetas <- prop.sig.thetas.fn(start.thetas,start.s2t)
   
   var.grid <- seq(min(wbars),max(wbars),length=100)
   vars <- current.vars <- B.basis(xs,knots.t)%*%exp(current.thetas)
   vfix <- exp(thetas[2])+exp(thetas[3])
   
   
   ### Prior and initialization of skewness
   e.grid <- seq(-4,4,length=100)
   if(error.fit=="skewnormal")
   {
     mu0.skewness <- 0
     sigmasq0.skewness <- 16
     skewness <- current.skewness <- start.skewness <- 0 
   }
   
   ### Prior and initialization for mixture
   if(error.fit=="mixture")
   {	
     a.u <- 1
     b.u <- 1
     sigma0.u <- 3
     simsize.mh.u <- 20
     z.u <- rep(1,sum(mis))
     k.u <- 1
     alpha.u <- 0.1	# robust but should not be set too large or small
     params.u <- matrix(c(0.5,0,1,1),nrow=20,ncol=4,byrow=T)		# unique values
     density.e.est <- numeric(length(e.grid))
     accepted <- numeric(simsize)
   }
   
   
   
   
   
   
   
   
   ###############################
   ### Storage for MCMC Output ###
   ###############################
   
   
   density.x.est <- numeric(length(x.grid))
   density.x0.est <- matrix(0,length(x.grid),simsize-burnin)
   var.est <- numeric(length(var.grid))
   e.grid <- seq(-4,4,length=100)
   density.e.est <- numeric(length(e.grid))
   proposed.log.prior <- current.log.prior <- 0
   proposed.log.proposal <- current.log.proposal <- 0
   no.clusters <- numeric(simsize)
   
   
   
   
   
   
   
   
   
   
   
   
   
   ##################
   ### Start MCMC ###
   ##################
   
   
   
   for (iii in 1:simsize)
   {
#     if(iii%%1000==0)
#       print(iii)
     
     
     
     ### Updating z.x
     
     student.x = (xs-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))
     marginal.tpdfs = dt(student.x,df=gama0.x)
     
     for(ii in 1:n)
     {
       prob.x = tabulate(z.x) 
       k.x = length(unique(z.x[-ii]))
       if((prob.x[z.x[ii]]==1)&&(z.x[ii]<max(z.x)))    # z.x[ii] appears only once in z.x
       {
         temp.mu.x = mu.x[z.x[ii]]
         temp.sigmasq.x = sigmasq.x[z.x[ii]]
         for(jj in z.x[ii]:k.x)
         {
           mu.x[jj] = mu.x[jj+1]
           sigmasq.x[jj] = sigmasq.x[jj+1]
         }
         z.x[z.x > z.x[ii]] = z.x[z.x > z.x[ii]]-1
         z.x[ii] = k.x+1         # Necessary steps
         mu.x[k.x+1] = temp.mu.x
         sigmasq.x[k.x+1] = temp.sigmasq.x  
       }
       prob.minus.x = tabulate(z.x[-ii])
       for(kk in 1:k.x)
       {
         likelihood = dnorm(xs[ii],mu.x[kk],sqrt(sigmasq.x[kk]))
         prob.x[kk] = prob.minus.x[kk] * likelihood
       }
       prob.x[k.x+1] = alpha.x * marginal.tpdfs[ii]
       prob.x = prob.x/sum(prob.x)    		# NOT really necessary        
       z.x[ii] = sample(k.x+1,1,TRUE,prob.x)   	# New z.x[ii] drawn
       if(z.x[ii]==(k.x+1))
       {
         #sigmasq.x[k.x+1] = 1/rgamma(1,shape=gama0.x,rate=sigmasq0.x)			# Normal-Inverse-Gamma
         sigmasq.x[k.x+1] = 1/rgamma(1,shape=gama0.x/2,rate=gama0.x*sigmasq0.x/2)	# Normal-Inverse-Scaled-Chisquare
         mu.x[k.x+1] = rnorm(1,mu0.x,sqrt(sigmasq.x[k.x+1]/nu0.x))
       }
     }
     
     
     ### Updating mu.x, sigmsq.x
     
     k.x = no.clusters[iii] = max(z.x)                # Number of clusters
     for(kk in 1:k.x)
     {
       temp = which(z.x==kk)
       xspool = xs[temp]
       
       nutemp.x = nu0.x + length(xspool)
       mutemp.x = (nu0.x*mu0.x + sum(xspool)) / (nutemp.x)
       gamatemp.x = gama0.x + length(xspool)				# Normal-Inverse-Scaled-Chisquare
       sigmasqtemp.x = (gama0.x*sigmasq0.x + sum(xspool^2) + nu0.x*mu0.x^2 - nutemp.x*mutemp.x^2) / gamatemp.x	
       sigmasq.x[kk] = 1/rgamma(1,shape=gamatemp.x/2,rate=(gamatemp.x*sigmasqtemp.x)/2)
       mu.x[kk] = rnorm(1,mutemp.x,sqrt(sigmasq.x[kk]/nutemp.x))
     }
     
     
     
     
     if(error.fit=="normal")
     {
       ### Updating xs (and us)
       
       proposed.xs = rtnorm(n,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
       proposed.vars = B.basis(proposed.xs,knots.t)%*%exp(thetas)
       
       prob.x = (tabulate(z.x)+alpha.x)/(alpha.x+n)
       proposed.prior = current.prior = numeric(n)
       for(kk in 1:k.x)
       {
         proposed.prior = proposed.prior + prob.x[kk]*dnorm(proposed.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
         current.prior = current.prior + prob.x[kk]*dnorm(current.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
       }
       
       mh.ratio = proposed.prior/current.prior
       mh.ratio = (   mh.ratio * dtnorm(current.xs,mean=proposed.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
                      / dtnorm(proposed.xs,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))   )
       
       temp.current.likelihood = dnorm(ws,mean=rep(current.xs,times=mis),sd=rep(sqrt(current.vars),times=mis))
       temp.proposed.likelihood = dnorm(ws,mean=rep(proposed.xs,times=mis),sd=rep(sqrt(proposed.vars),times=mis))
       
       current.likelihood = tapply(temp.current.likelihood,inds,"prod")
       proposed.likelihood = tapply(temp.proposed.likelihood,inds,"prod")
       
       mh.ratio = mh.ratio * proposed.likelihood/current.likelihood	
       mh.ratio[is.nan(mh.ratio)] = 0
       
       u = runif(n)
       inds.to.replace =(1:n)[u<mh.ratio]
       xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
       vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
       
       us = ws - rep(xs,times=mis)
       
       
       ### Updating theta
       
       proposed.thetas = rmvnorm(1,current.thetas,diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas)
       proposed.thetas = t(proposed.thetas)
       proposed.vars = B.basis(xs,knots.t)%*%exp(proposed.thetas)
       
       current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
       proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
       
       temp.current.likelihood = dnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(current.vars),times=mis))
       temp.proposed.likelihood = dnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(proposed.vars),times=mis))
       
       current.log.likelihood = sum(log(temp.current.likelihood))
       proposed.log.likelihood = sum(log(temp.proposed.likelihood))
       
       log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
       
       log.u = log(runif(1))
       if(log.u<log.mh.ratio)
       {
         thetas = current.thetas = proposed.thetas
         vars = current.vars = proposed.vars
       }
     }
     
     
     
     
     
     if(error.fit=="skewnormal")
     {
       ### Updating xs (and us)
       
       proposed.xs = rtnorm(n,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
       proposed.vars = B.basis(proposed.xs,knots.t)%*%exp(thetas)
       
       prob.x = tabulate(z.x)/(n)
       proposed.prior = current.prior = numeric(n)
       for(kk in 1:k.x)
       {
         proposed.prior = proposed.prior + prob.x[kk]*dnorm(proposed.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
         current.prior = current.prior + prob.x[kk]*dnorm(current.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
       }
       
       mh.ratio = proposed.prior/current.prior
       mh.ratio = (   mh.ratio * dtnorm(current.xs,mean=proposed.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
                      / dtnorm(proposed.xs,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))   )
       
       temp.current.likelihood = dskewnorm(ws,mean=rep(current.xs,times=mis),sd=rep(sqrt(current.vars),times=mis),skewness)
       temp.proposed.likelihood = dskewnorm(ws,mean=rep(proposed.xs,times=mis),sd=rep(sqrt(proposed.vars),times=mis),skewness)
       
       current.likelihood = tapply(temp.current.likelihood,inds,"prod")
       proposed.likelihood = tapply(temp.proposed.likelihood,inds,"prod")
       
       mh.ratio = mh.ratio * proposed.likelihood/current.likelihood	
       mh.ratio[is.nan(mh.ratio)] = 0
       
       u = runif(n)
       inds.to.replace =(1:n)[u<mh.ratio]
       xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
       vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
       
       us = ws - rep(xs,times=mis)
       
       
       ### Updating theta
       
       proposed.thetas = rmvnorm(1,current.thetas,(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas))
       proposed.thetas = t(proposed.thetas)
       proposed.vars = B.basis(xs,knots.t)%*%exp(proposed.thetas)
       
       current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
       proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
       
       temp.current.likelihood = dskewnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(current.vars),times=mis),skewness)
       temp.proposed.likelihood = dskewnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(proposed.vars),times=mis),skewness)
       
       current.log.likelihood = sum(log(temp.current.likelihood))
       proposed.log.likelihood = sum(log(temp.proposed.likelihood))
       
       log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
       
       log.u = log(runif(1))
       if(log.u<log.mh.ratio)
       {
         thetas = current.thetas = proposed.thetas
         vars = current.vars = proposed.vars
       }
       
       
       ### Updating skewness
       
       proposed.skewness <- rnorm(1,current.skewness,sd=1)
       
       current.log.prior = - (current.skewness-mu0.skewness)^2/(2*sigmasq0.skewness)
       proposed.log.prior = - (proposed.skewness-mu0.skewness)^2/(2*sigmasq0.skewness)
       
       temp.current.likelihood = dskewnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(vars),times=mis),current.skewness)
       temp.proposed.likelihood = dskewnorm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(vars),times=mis),proposed.skewness)
       
       current.log.likelihood = sum(log(temp.current.likelihood))
       proposed.log.likelihood = sum(log(temp.proposed.likelihood))
       
       log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
       
       log.u = log(runif(1))
       if(log.u<log.mh.ratio)
         skewness = current.skewness = proposed.skewness
     }
     
     
     
     
     
     if(error.fit=="mixture")
     {
       ### Updating xs (and us)
       
       proposed.xs = rtnorm(n,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
       proposed.vars = B.basis(proposed.xs,knots.t)%*%exp(thetas)
       
       prob.x = tabulate(z.x)/(n)
       proposed.prior = current.prior = numeric(n)
       for(kk in 1:k.x)
       {
         proposed.prior = proposed.prior + prob.x[kk]*dnorm(proposed.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
         current.prior = current.prior + prob.x[kk]*dnorm(current.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
       }
       
       mh.ratio = proposed.prior/current.prior
       mh.ratio = (   mh.ratio * dtnorm(current.xs,mean=proposed.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))
                      / dtnorm(proposed.xs,mean=current.xs,sd=diff(range(start.xs))/6,lower=min(knots.t),upper=max(knots.t))   )
       
       temp.current.likelihood = d.restricted.mix.norm(ws,mean=rep(current.xs,times=mis),sd=rep(sqrt(current.vars),times=mis),params.u[z.u,])
       temp.proposed.likelihood = d.restricted.mix.norm(ws,mean=rep(proposed.xs,times=mis),sd=rep(sqrt(proposed.vars),times=mis),params.u[z.u,])
       
       current.likelihood = tapply(temp.current.likelihood,inds,"prod")
       proposed.likelihood = tapply(temp.proposed.likelihood,inds,"prod")
       
       mh.ratio = mh.ratio * proposed.likelihood/current.likelihood	
       mh.ratio[is.nan(mh.ratio)] = 0
       
       u = runif(n)
       inds.to.replace =(1:n)[u<mh.ratio]
       xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
       vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
       
       us = ws - rep(xs,times=mis)
       
       
       ### Updating theta
       
       proposed.thetas = rmvnorm(1,current.thetas,prop.sig.thetas)
       while(proposed.thetas[3]>log(vfix))
         proposed.thetas = rmvnorm(1,current.thetas,sig.tune.thetas.1 * prop.sig.thetas)
       proposed.thetas[2] = log(vfix-exp(proposed.thetas[3]))
       proposed.thetas = t(proposed.thetas)
       proposed.vars = B.basis(xs,knots.t)%*%exp(proposed.thetas)
       
       current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
       proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
       
       temp.current.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(current.vars),times=mis),params.u[z.u,])
       temp.proposed.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=sum(mis)),sd=rep(sqrt(proposed.vars),times=mis),params.u[z.u,])
       
       current.log.likelihood = sum(log(temp.current.likelihood))
       proposed.log.likelihood = sum(log(temp.proposed.likelihood))
       
       log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
       
       log.u = log(runif(1))
       if(log.u<log.mh.ratio)
       {
         thetas = current.thetas = proposed.thetas
         vars = current.vars = proposed.vars
       }
       
       
       
       ### Updating z.u
       
       for(ii in 1:sum(mis))
       {
         prob.u = tabulate(z.u) 
         k.u = length(unique(z.u[-ii]))
         if((prob.u[z.u[ii]]==1)&&(z.u[ii]<max(z.u)))    # z.u[ii] appears only once in z.u
         {
           temp.params.u = params.u[z.u[ii],]
           for(jj in z.u[ii]:k.u)         	# Necessary step			
             params.u[jj,] = params.u[jj+1,]
           z.u[z.u > z.u[ii]] = z.u[z.u > z.u[ii]]-1
           z.u[ii] = k.u+1         		# Necessary step 
           params.u[k.u+1,] = temp.params.u
         }
         
         prob.minus.u <- tabulate(z.u[-ii])
         prob.minus.u[k.u+1] = alpha.u
         prob.u = prob.minus.u/sum(prob.minus.u)    # NOT really necessary
         
         proposed.z.u = sample(k.u+1,1,TRUE,prob.u)   # z.u[ii] proposed
         
         if(z.u[ii]!=proposed.z.u)
         {
           if(proposed.z.u==(k.u+1))
             proposed.params.u <- r.proposal.params.restricted.mix.norm(1,1,3,3,2.5,3,2.5)
           else
             proposed.params.u <- params.u[proposed.z.u,]
           
           proposed.likelihood = d.restricted.mix.norm(us[ii],mean=0,sd=sqrt(current.vars[inds[ii]]),proposed.params.u)
           current.likelihood = d.restricted.mix.norm(us[ii],mean=0,sd=sqrt(current.vars[inds[ii]]),params.u[z.u[ii],])
           acc.prob <- proposed.likelihood/current.likelihood
           if(runif(1)<acc.prob)
           {
             z.u[ii] <- proposed.z.u
             params.u[z.u[ii],] <- proposed.params.u
           }
         }
       }
       
       
       
       ### Updating params.u
       
       k.u = max(z.u)                # Number of clusters
       if(iii>2000)
         simsize.mh.u = 1
       for(rr in 1:simsize.mh.u)
       {
         for(kk in 1:k.u)
         {
           temp = which(z.u==kk)
           uspool = us[temp]
           varspool = vars[inds[temp]]
           
           proposed.params.u <- r.tnorm.proposal.params.restricted.mix.norm(params.u[kk,])
           
           proposed.log.likelihood = sum(log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.u)))
           current.log.likelihood = sum(log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.u[kk,])))
           
           proposed.log.prior = log(dnorm(proposed.params.u[2],mean=0,sd=sigma0.u))+sum(log(dgamma(1/proposed.params.u[3:4],shape=a.u,scale=b.u)))
           current.log.prior = log(dnorm(params.u[kk,2],mean=0,sd=sigma0.u))+sum(log(dgamma(1/params.u[kk,3:4],shape=a.u,scale=b.u)))
           
           proposed.log.proposal = log(dtnorm(params.u[kk,1],mean=proposed.params.u[1],sd=0.01))+sum(log(dtnorm(params.u[kk,3:4],mean=proposed.params.u[3:4],sd=0.1,lower=0,upper=Inf)))
           current.log.proposal = log(dnorm(proposed.params.u[1],mean=params.u[kk,1],sd=0.01))+sum(log(dtnorm(proposed.params.u[3:4],mean=params.u[kk,3:4],sd=0.1,lower=0,upper=Inf)))
           
           log.acc.prob <- (proposed.log.likelihood-current.log.likelihood) + (proposed.log.prior-current.log.prior) + (proposed.log.proposal-current.log.proposal)
           # Only the likelihood term dominates. To improve speed, the second and the third terms can be neglected. 
           
           if(log(runif(1))<log.acc.prob)
           {
             params.u[kk,] <- proposed.params.u
             accepted[iii] <- 1
           }
         } 
       }
     }
     
     
     
     
     
     
     ### Updating s2t
     
     s2t = 1/rgamma(1,shape=alpha.t+(K.t+1)/2,rate=beta.t+t(thetas)%*%P.t%*%thetas/2)
     
     
     
     
     
     
     if(iii>burnin)
     {
       prob.x = tabulate(z.x)/(alpha.x+n)
       prob.x[k.x+1] = alpha.x/(alpha.x+n)
       for(kk in 1:k.x)
         density.x0.est[,iii-burnin] = density.x0.est[,iii-burnin] + prob.x[kk]*dnorm(x.grid,mu.x[kk],sqrt(sigmasq.x[kk]))
       density.x0.est[,iii-burnin] = density.x0.est[,iii-burnin] + prob.x[k.x+1]*dt(student.x.grid,df=gama0.x)
       
       if(error.fit!="mixture")
         var.est = var.est + B.basis(var.grid,knots.t)%*%exp(thetas)
       
       if(error.fit=="skewnormal")
       {
         density.e.est = density.e.est + dskewnorm(e.grid,mean=rep(0,times=length(e.grid)),rep(1,times=length(e.grid)),skewness)
       }
       
       if(error.fit=="mixture")
       {
         k.u = max(z.u)                # Number of clusters
         prob.e <- tabulate(z.u)/sum(mis)
         var.e = var.e.fn(prob.e[1:k.u],params.u[1:k.u,])
         var.est = var.est + B.basis(var.grid,knots.t)%*%exp(thetas) * var.e
         density.e.est <- density.e.est + d.scaled.restricted.mix.norm(e.grid,0,1,prob.e[1:k.u],params.u[1:k.u,])
       }
     }
     
   }
   
   
   
   
   
   
   density.x.est <- rowMeans(density.x0.est)
   var.est <- var.est/(simsize-burnin)
   if(error.fit=="normal")
     density.e.est <- dnorm(e.grid)
   if(error.fit!="normal")
     density.e.est <- density.e.est/(simsize-burnin)
   


if(plot_results=="TRUE")
	{
	par(mfrow=c(3,1))
	plot(x.grid[1:250],density.x.est[1:250],ylim=c(0,max(density.x.est)),type="l",col="blue",lwd=3,xlab="",ylab="")
	#plot(x.grid,density.x.est,ylim=c(0,max(density.x.est)),type="l")
	title(main=paste("Simulated ",variable.name," True X Estimated Density",sep=""))

	plot(e.grid,density.e.est,type="l",col="blue",lwd=3,xlab="",ylab="")
	points(e.grid,dnorm(e.grid),ylim=c(0,max(density.e.est)),type="l",lty=2)
	if (error.fit == "normal"){
	  title(main=paste("Simulated ",variable.name," Estimated Assumed Normal Error Density",sep=""))
	}
	if (error.fit == "skewnormal"){
	  title(main=paste("Simulated ",variable.name," Estimated Assumed Skew-Normal Error Density",sep=""))
	}
	if (error.fit == "mixture"){
	  title(main=paste("Simulated ",variable.name," Estimated Nonparametric Error Density",sep=""))
	}
	
	plot(var.grid[1:50],var.est[1:50],ylim=c(0,max(var.est[1:50])),type="l",col="blue",lwd=3,xlab="",ylab="")
	#plot(var.grid,var.est,ylim=c(0,max(var.est)),type="l")
	points(wbars,s2is,xlim=c(min(var.grid),max(var.grid)),pch="*")
	title(main=paste("Simulated ",variable.name," Estimated Variance Function",sep=""))
	
	par(mfrow=c(1,1))
	}


	filename <- paste("Data",variable.name,error.fit,simsize,burnin,sep="_")
	filename <- paste(filename,".RData",sep="")

return(list(x.grid=x.grid,density.x.est=density.x.est,e.grid,density.e.est,var.grid,var.est,wbars,s2is,file=filename))
}







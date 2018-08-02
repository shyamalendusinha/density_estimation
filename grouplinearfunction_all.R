
# The whole code was combined version of 3 codes downloaded from https://github.com/MaZhuang/grouplinear/(functions.R,functions_XKB.R,dynamic_sure.R)

# Group-linear Functions

# "2014-11-13 11:33:34 EST"

## spherically symmetric estimator with c_n = c^*_n
spher <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    x. - bhat*(x. - mean(x.))
  }
}


## spherically symmetric estimator with c_n = c^*_n, shrinkage toward zero
spher.zero <- function(x.,v.){
  n. <- length(x.)
  cstar <- max( 1-2*( max(v.)/mean(v.) )/n., 0)
  bhat <- min( cstar*mean(v.)/mean(x.^2), 1 )
  (1- bhat)*x.
}

## function that returns the common bhat (replicated)
spher.bhat <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    return(rep(bhat,n.))
  }
}

## group-linear estimator

grouplinear <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

## group-linear estimator with shrinkage toward zero

grouplinear.zero <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher.zero,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

# # ## example
# n <- 300
# v <- runif(n,.1,1)
# theta <- v-mean(v)
# x <- rnorm(n,theta,sd=sqrt(v))
# grouplinear(x,v)
# mean( (grouplinear(x,v)-theta)^2 )   
# mean( (grouplinear.zero(x,v)-theta)^2 )   


## sure for grouplinear estimator (these are estimates of risk -- not of theta)
sure.spher <- function(x.,v.){
  n. <- length(x.)
  # cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0) ##modified
  if (n.==0) {0 
  } else if ( (n.<3) ) {sum(v.)  #| (var(x.)==0) 
  }
  else if (max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)==0){
    sum(v.) 
  }
  else if (var(x.)==0){
    (2-n.)/n.*sum(v.)+sum((x.-mean(x.))^2)
  }
  else {	# can set sure to an arbitrary value if var(x.)=0, since this event is of measure zero
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0) ##modified
    b <- cstar * mean(v.)/var(x.)
    b <- min(1,b)
    db <- -cstar * mean(v.)/(var(x.))^2 * as.numeric( cstar * mean(v.)/var(x.) < 1 )##
    sum(   v. + ( b * (x.-mean(x.)) )^2 - 2 * v. * (  (1-1/n.) * b + 2 * (x.-mean(x.))^2 * db/(n.-1)  )   )
  }
}



sure.spher.zero <- function(x.,v.){
  n. <- length(x.)
  if (n.==0) {0 
  }else if ( (n.==1) | (var(x.)==0) ) {sum(v.) 
  }else {	# can set sure to an arbitrary value if var(x.)=0, since this event is of measure zero
    cstar <- 1-2*max(v.)/sum(v.) * ( 1-2*max(v.)/sum(v.) > 0 )
    b <- cstar * mean(v.)/mean(x.^2)
    b <- min(1,b)
    db <- -cstar * mean(v.)/mean(x.^2)^2 * ( cstar * mean(v.)/mean(x.^2) < 1 )
    sum(   v. + ( b * (x.) )^2 - 2 * v. * (  b + 2 * x.^2 * db/n.  )   )
  }
}


sure <- function(nbreak,x,v){ #nbreak=num of bins
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak)
  xsub <- split(x,splitby,drop=T)
  vsub <- split(v,splitby,drop=T)
  suresub <- mapply(sure.spher,xsub,vsub)   
  sum(suresub)/n
}

sure.zero <- function(nbreak,x,v){ #nbreak=num of bins
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak)
  xsub <- split(x,splitby,drop=T)
  vsub <- split(v,splitby,drop=T)
  suresub <- mapply(sure.spher.zero,xsub,vsub)
  sum(suresub)/n
}


## sure grouplinear (equal-bins)
grouplinear.sure <- function(x,v,kmax=10){ #kmax is maximum bin-count to search over
  sure.vec <- c(sure.spher(x,v), sapply(X=2:kmax,FUN=sure,x,v))
  khat.sure <- which.min(sure.vec)
  est <- if(khat.sure>1) grouplinear( x,v,nbreak=khat.sure) else spher(x,v)
}

grouplinear.sure.zero <- function(x,v,kmax=10){ #kmax is maximum bin-count to search over
  sure.vec <- c(sure.spher.zero(x,v), sapply(X=2:kmax,FUN=sure.zero,x,v))
  khat.sure <- which.min(sure.vec)
  est <- if(khat.sure>1) grouplinear( x,v,nbreak=khat.sure) else spher(x,v)
}

# XKB (2012) code

library(isotone)

# bayes rule for fixed lambda,mu
thetahat <- function(X,A,lambda,mu){
  lambda/(lambda+A) * X + A/(lambda+A) * mu
}

# d/dlambda{SURE(lambda,mu=mu.hat.SURE(lambda))} (proportional to)
g <- function(lambda,X,A){  
  sum( A^2/(lambda+A)^3 * (X-(  sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )  ))^2 - A^2/(lambda+A)^2 )
  #equivalent to the following(which is just easier to read):
  #mu <- sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )
  #sum( A^2/(lambda+A)^3 * (X-mu)^2 - A^2/(lambda+A)^2 )
}

# SURE(lambda,mu)
f <- function(par,X,A){  
  lambda <- par[1]
  mu <- par[2]
  sum(  A/(lambda+A)^2 * ( A * (X-mu)^2 + lambda^2 - A^2 )  )
}

# SURE.G(lambda)
f.G <- function(lambda,X,A){  
  sum(  ( A/(A+lambda) )^2 * (X - mean(X))^2 + A/(A+lambda) * (lambda - A + 2/p * A)  )
}


thetahat.M <- function(X,A){
  lambda.sure <- ifelse( g(0,X=X,A=A)*g(max(A)*1000,X=X,A=A) < 0, uniroot(g,c(0,max(A)*1000),X=X,A=A, tol=1e-9)$root, optim(c(mean(   pmax(  ( X-mean(X) )^2 - 	A,0  )   ), mean(X)),f,X=X,A=A, method = "L-BFGS-B",lower=c(0,-Inf))$par[1] )
  mu.sure <- sum( A^2/(lambda.sure+A)^2 * X ) / sum( A^2/(lambda.sure+A)^2 )
  thetahat(X,A,lambda.sure,mu.sure)
}

thetahat.G <- function(X,A){
  lambda <- optimize(f.G,lower=0,upper=1000,X=X,A=A)$minimum
  thetahat(X,A,lambda,mean(X))
}

thetahat.SG <- function(X,A){
  p <- length(X)
  fit <- gpava( z = A, y = A * (1-1/p) / (X-mean(X))^2, weights = (X-mean(X))^2, solver = weighted.mean, ties="primary" )
  bhat <- pmin(  pmax( fit$x,0 ),1  )
  (1-bhat) * X + bhat * mean(X)
}

#example
#p <- 100
#A <- runif(p,.1,1)
#theta <- rnorm(p)
#X <- rnorm(p,theta,sqrt(A))
#thetahat.M(X,A)
#thetahat.G(X,A)
#thetahat.SG(X,A)


# plot(X,thetahat.M(X,A), cex=.5, pch=16, ylim = range(thetahat.M(X,A), thetahat.SG(X,A)))
# points(X,thetahat.SG(X,A), cex=.5, pch=16, col='blue')


# SURE^G(lambda)

sure.G <- function(lambda,x,v){
  x.bar <- mean(x)
  mean(  v^2/(v+lambda)^2 * (x-x.bar)^2 + v/(v+lambda) * (lambda - v + 2/n * v)  )  
}

# lambda.vec <- ppoints(100)*100
# sure.G(lambda.vec,X,A)
# y <- sapply(lambda.vec,sure.G,x=x,v=v)
# plot(lambda.vec,y,pch=16,cex=.5)

# extended James-Stein (equation (7.3) in Xie et al, 2012)
JS <- function(X,A){
  p <- length(X)
  muhat <- sum( X*(1/A) )/sum(1/A)
  b <- 1 - (p-3) / sum( (X-muhat)^2/A )
  return( muhat + max(0,b) * (X-muhat) )
}



DynamicSure=function(x,v){
  n <- length(x)
  a=matrix(rep(0,n*n),ncol=n) ##separation
  b=a ##value
  for(i in 1:n){
    a[i,i]=i
    b[i,i]=v[i]
  }
  for (l in 1:(n-1)){
    # if (l %% 100==0){
    #    print(l)
    # }
    for (i in 1:(n-l)){
      j=l+i
      sure=sure.spher(x[i:j], v[i:j])
      #	print(sure)
      a[i,j]=j
      b[i,j]=sure
      for (k in i:(j-1)){
        #temp=sure.spher(x[i:k], v[i:k])+sure.spher(x[(k+1):j], v[(k+1):j])
        temp=b[i,k]+b[k+1,j]
        if (b[i,j]>temp){
          #	print(sure)
          a[i,j]=k
          b[i,j]=temp
        }
      }
      
    }
  }
  list(a,b)
}

DynamicSureMin=function(x,v,d=40){
  d=floor(d)
  n <- length(x)
  a=matrix(rep(0,n*n),ncol=n) ##separation
  b=a ##value
  for (i in 1:(n-d+1)){
    j=i+d-1
    a[i,j]=j
    b[i,j]=sure.spher(x[i:j], v[i:j])
  }
  for (l in d:(n-1)){
    # if (l %% 100==0){
    #    print(l)
    # }
    for (i in 1:(n-l)){
      j=l+i
      sure=sure.spher(x[i:j], v[i:j])
      #	print(sure)
      a[i,j]=j
      b[i,j]=sure
      if ((i+d-1)<=(j-d)){
        for (k in (i+d-1):(j-d)){
          #temp=sure.spher(x[i:k], v[i:k])+sure.spher(x[(k+1):j], v[(k+1):j])
          temp=b[i,k]+b[k+1,j]
          if (b[i,j]>temp){
            #	print(sure)
            a[i,j]=k
            b[i,j]=temp
          }
        }}
      
    }
  }
  list(a,b)
}



partition=function(position,i,j){
  if (position[i,j]==j){
    return(j)
  }else if (position[i,j]==i){
    return(i)
  }
  else{
    a=partition(position,i,position[i,j])
    b=partition(position,position[i,j],j)
    return(c(a,position[i,j],b))
  }
}




#################shrink towards the mean in the bin

dynamic.grouplinear <- function(x,v,group){ #nbreak=num of bins
  ngroup <- length(group)
  n=length(x)
  est=rep(0,n)
  for (i in 1:(ngroup-1)){
    est[(group[i]+1):group[i+1]]=spher(x[(group[i]+1):group[i+1]],v[(group[i]+1):group[i+1]])
  }
  est
}

GroupSure<- function(x,v){ 
  c=DynamicSure(x,v)
  position=c[[1]]
  n=dim(position)[1]
  group=partition(position,1,n)
  group=c(0, group,n)
  group=unique(group)
  est=dynamic.grouplinear(x,v,group)
  return(est)
}

GroupSureMin<- function(x,v,d){ 
  d=floor(d)
  c=DynamicSureMin(x,v,d)
  position=c[[1]]
  n=dim(position)[1]
  group=partition(position,1,n)
  group=c(0, group,n)
  group=unique(group)
  est=dynamic.grouplinear(x,v,group)
  return(est)
}

#################shrink towards zero
dynamic.grouplinear.zero <- function(x,v,group){ #nbreak=num of bins
  ngroup <- length(group)
  n=length(x)
  est=rep(0,n)
  for (i in 1:(ngroup-1)){
    est[(group[i]+1):group[i+1]]=spher.zero(x[(group[i]+1):group[i+1]],v[(group[i]+1):group[i+1]])
  }
  est
}


GroupSure.zero<- function(x,v){ 
  c=DynamicSure(x,v)
  position=c[[1]]
  n=dim(position)[1]
  group=partition(position,1,n)
  group=c(0, group,n)
  group=unique(group)
  est=dynamic.grouplinear.zero(x,v,group)
  return(est)
}


GroupSureMin.zero<- function(x,v,d){ 
  d=floor(d)
  c=DynamicSureMin(x,v,d)
  position=c[[1]]
  n=dim(position)[1]
  group=partition(position,1,n)
  group=c(0, group,n)
  group=unique(group)
  est=dynamic.grouplinear.zero(x,v,group)
  return(est)
}

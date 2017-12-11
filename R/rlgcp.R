

.set.cov <- function(separable,model,param,sigma2)
  {
    mods <- 0
    modt <- 0
    mod <- 0    

    models <- c("exponential","cauchy","stable","wave","gneiting","cesare","matern","none")

    for(i in 1:length(model))
      {
        M <- which(models==model[i])
        if (length(M)==0) stop("the model is not implemented")
      }

    for (i in 1:length(unique(model)))
      {
        if (((isTRUE(separable)) & ((model[i]==models[5]) | (model[i]==models[6]))) | ((!(isTRUE(separable))) & ((model[i]==models[1]) | (model[i]==models[2]) | (model[i]==models[3]) | (model[i]==models[4]) | (model[i]==models[7])))) stop("'stcov' does not match with 'model'")
      }

    if (isTRUE(separable))
      {
        if ((length(model)!=1) & (length(model)!=2))
          stop("for separable covariance functions, 'model' must be of length 1 or 2")
        if (length(model)==1)
          {
            if (model=="none")
              {
                mods <- 0
                modt <- 0
              }
            if (model=="exponential")
              {
                mods <- 1
                modt <- 1
              }
            if (model=="stable")
              {
                mods <- 2
                if ((param[1] >1) | (param[1]<0)) stop("Stable model parameter must lie in (0,1]")
                modt <- 2
                if ((param[2] >1) | (param[2]<0)) stop("Stable model parameter must lie in (0,1]")
              }
            if (model=="cauchy")
              {
                mods <- 3
                if (param[1]<=0) stop("Cauchy model parameter must be strictly positive")
                modt <- 3
                if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
              }
            if (model=="wave")
              {
                mods <- 4
                modt <- 4
                }
            if (model=="matern") 
			{
                  mods <- 7
			if (param[2]<=0 | param[1]<=0) stop("Matern model parameters must be strictly positive")
                  modt <- 7
			if (param[3]<=0 | param[4]<=0) stop("Matern model parameters must be strictly positive")
                  }
          }
            if (length(model)==2)
              {
                if (model[1]=="none")
                    mods <- 0
                if (model[2]=="none")
                  modt <- 0
                if (model[1]=="exponential")
                    mods <- 1
                if (model[2]=="exponential")
                  modt <- 1
                if (model[1]=="stable")
                    {
                      mods <- 2
                      if ((param[1] >1) | (param[1]<0)) stop("Stable model parameter must lie in (0,1]")
                    }
                if (model[2]=="stable")
                  {
                    modt <- 2
                    if ((param[2] >1) | (param[2]<0)) stop("Stable model parameter must lie in (0,1]")
                  }
                if (model[1]=="cauchy")
                  {
                    mods <- 3
                    if (param[1]<=0) stop("Cauchy model parmaeter must be strictly positive")
                  }
                if (model[2]=="cauchy")
                  {
                    modt <- 3
                    if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
                  }
                if (model[1]=="wave")
                    mods <- 4
                if (model[2]=="wave")
                  modt <- 4
                if (model[1]=="matern")
			{
                  mods <- 7
			if (param[3]<=0 | param[1]<=0) stop("Matern model parameters must be strictly positive")
                  }
		    if (model[2]=="matern")
			{
                  modt <- 7
			if (param[2]<=0 | param[4]<=0) stop("Matern model parameters must be strictly positive")
                  }
              }
      }
    if (!(isTRUE(separable)))
      {
        if (length(model)!=1)
          stop("for non-separable covariance functions, 'model' must be of length 1")
        if (model=="gneiting")
          {
            mod <- 5
            if (param[6]<1) stop("for Gneiting's covariance function, the sixth parameter must be greater than 1")
            if ((param[3]<=0) | (param[3]>1)) stop("for Gneiting's covariance function, the third parameter must lie in (0,1]")
            if ((param[4]<=0) | (param[4]>1)) stop("for Gneiting's covariance function, the fourth parameter must lie in (0,1]")
            if ((param[5]!=1) & (param[5]!=2) & (param[5]!=3)) stop("for Gneiting's covariance function, the fifth parameter must be 1, 2 or 3")
            if ((param[2]!=1) & (param[2]!=2) & (param[2]!=3)) stop("for Gneiting's covariance function, the second parameter must be 1, 2 or 3")
            if ((param[2]==1) & ((param[1]<0) | (param[1]>2))) stop("for Gneiting's covariance function, if the second parameter equals 1, the first parameter must lie in [0,2]") 
            if ((param[2]==2) & (param[1]<=0)) stop("for Gneiting's covariance function, if the second parameter equals 2, the first parameter must be strictly positive")            
          }
        if (model=="cesare")
          {
            mod <- 6
            if (((param[1]>2) | (param[1]<1)) | ((param[2]>2) | (param[2]<1))) stop("for De Cesare's model, the first and second parameters must lie in [1,2]")
            if (param[3]<3/2) stop("for De Cesare's model, the third parameter must be greater than 3/2")
          }
      }

    return(model=c(mods,modt,mod))
  }

.matern = function (d, scale = 1, alpha = 1, nu = 0.5) 
{
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    d <- d * alpha
    d[d == 0] <- 1e-10
    k <- 1/((2^(nu - 1)) * gamma(nu))
    res <- scale * k * (d^nu) * besselK(d, nu)
    return(res)
}


.covst <- function(x,y,times,separable=TRUE,model,param=c(1,1,1,1,1,1),sigma2=1,scale=c(1,1),plot=TRUE,nlevels=10,aniso=0,ani.pars=c(0,1))
{

  nt <- length(times)
  nx = length(x)
  ny = length(y)

  model <- .set.cov(separable,model,param,sigma2)
  
  gs <- array(0, dim = c(nx,ny,nt))
  storage.mode(gs) <- "double"

  gs <- .Fortran("covst",
                 (gs),
                 as.double(x),
                 as.double(y),
                 as.double(times),
                 as.integer(nx),
                 as.integer(ny),
                 as.integer(nt),
                 as.integer(model),
                 as.double(param),
                 as.double(sigma2),
                 as.double(scale),
                 as.double(aniso),
                 as.double(ani.pars))[[1]]


  if (plot==TRUE)
    {
      image(x,y,gs[,,2],col=grey((1000:1)/1000),xlab="x",ylab="y",cex.axis=1.5,cex.lab=2,font=2)
      contour(x,y,gs[,,2],add=T,col=4,labcex=1.5,nlevels=nlevels)

    }
  
  invisible(return(list(x=x,y=y,times=times,gs=gs)))

}


.gauss3D <- function(nx=100,ny=100,nt=100,xlim,ylim,tlim,separable=TRUE,model="exponential",param=c(1,1,1,1,1,1),scale=c(1,1),var.grf=1,mean.grf=0,exact=TRUE,aniso=0,ani.pars=c(0,1))
{

  N <- c(nx,ny,nt)

  mod <- .set.cov(separable,model,param,var.grf)
  
  g <- floor(log(2*(N-1))/log(2))+1
  M <- 2^g

  count <- 1
  changG <- TRUE
  while(changG==TRUE)
    {
      L <- rep(-9999,M[1]*M[2]*M[3])  
      
      storage.mode(L) <- "double"
      
      res <- .Fortran("circ",
                      (L),
                      as.integer(M),
                      as.integer(N),
                      as.double(xlim),
                      as.double(ylim),
                      as.double(tlim),
                      as.integer(mod), 
                      as.double(param),
                      as.double(var.grf),
                      as.double(scale),
                      as.double(aniso),
                      as.double(ani.pars))[[1]]
 
      L <- array(res,dim=M)
      FTL <- fft(L)

      if (isTRUE(exact))
        {      
          if (min(Re(FTL))<0)
            {
              g <- g+1
              M <- 2^g
              changG <- TRUE
              count <- count+1
            }
          else
            changG <- FALSE
        }
      else
        {
          FTL[Re(FTL)<0] <- 0
          changG <- FALSE
        }

    }
  print(count)
  
#  X <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#  X <- fft(X,inverse=TRUE)
#  A <- sqrt(FTL)*X
#  G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
#  G <- G+mean.grf

## Remarks
# --------
#
#  The right way is:
    X1 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
    X2 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
    X <- fft(X1+1i*X2,inverse=TRUE)
    A <- sqrt(FTL)*X
    G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
#  but it doesn't change the results and it is slightly faster. 
#  
# Taking the first elements of FTL and then computing G (changing M by N)
# is faster than taking the first elements of G, but it does not provide
# the same results
  
  invisible(return(list(G=G,L=L)))
}

.anisotropy.rotation <- function(theta=0,dzeta=1,inverse=TRUE) {
  rotate <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
  scale <- diag(c(1,1/dzeta))
  Q.inv <- scale%*%rotate
  if(inverse) return(Q.inv)
  solve(Q.inv)
}


rlgcp <- function(s.region, t.region, replace=TRUE, npoints=NULL, nsim=1, nx=100, ny=100, nt=100,separable=TRUE,model="exponential",param=c(1,1,1,1,1,1),scale=c(1,1),var.grf=1,mean.grf=0,lmax=NULL,discrete.time=FALSE,exact=FALSE,anisotropy=FALSE,ani.pars=NULL)
{
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)
  if (anisotropy==TRUE & length(ani.pars) != 2) 
    stop("Argument ani.pars must be a vector with 2 elements: the anisotropy angle and anisotropy ratio")
  aniso = ifelse(anisotropy,1,0)

  ###############
  
  stcov=.covst(x=seq(-1,1,length=nx),y=seq(-1,1,length=ny),times=seq(0,1,length=nt),separable=separable,model,param=param,sigma2=var.grf,scale=scale,plot=FALSE,nlevels=10,aniso=aniso,ani.pars=ani.pars)

  
  ###############
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]
  tau <- c(start=t.region[1],end=t.region[2],step=(t.region[2]-t.region[1])/(nt-1))
#  bpoly <- bbox(s.region)

  lambdamax <- lmax
  pattern <- list()
  Lambdafin <- list()
  ni <- 1

  s.grid <- .make.grid(nx,ny,s.region)
  s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)

  if (discrete.time==TRUE)
    {
      vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
      if (nt>length(vect))
        {
          nt <- length(vect)
          warning("nt used is less than the one given in argument")
          t.grid <- list(times=vect,tinc=1)
        }
      else
        {
          vect <- round(seq(floor(t.region[1]),ceiling(t.region[2]),length=nt))
          t.grid <- list(times=vect,tinc=round(t.area/(nt-1)))
        }
    }
  else
    t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=(t.area/(nt-1)))

  while(ni<=nsim)
    {
      S <- .gauss3D(nx=nx,ny=ny,nt=nt,xlim=range(s.region[,1]),ylim=range(s.region[,2]),tlim=range(t.region),separable=separable,model=model,param=param,scale=scale,var.grf=var.grf,mean.grf=mean.grf,exact=exact,aniso=aniso,ani.pars=ani.pars)
      covar = S$L
      S = S$G
      
      Lambda <- exp(S)

      mut <- rep(0,nt)
      for (it in 1:nt)
        {
          Lambda[,,it][s.grid$mask==FALSE] <- NaN
          mut[it] <- sum(Lambda[,,it],na.rm=TRUE)
        }
      
      if (is.null(npoints))
        {
          en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc*t.grid$tinc
          npoints <- round(rpois(n=1,lambda=en),0)
        }

      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=TRUE)
  
      npts <- round(lambdamax/(s.area*t.area),0)
      if (npts==0) stop("there is no data to thin")

      if ((replace==FALSE) & (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
          t.grid$times = times.init
        }
###      else
###        times.init <- runif(nt,min=t.region[1],max=t.region[2])

      XX=rep(s.grid$X,nt)
      YY=rep(s.grid$Y,nt)
      TT=rep(t.grid$times,each=length(s.grid$X))
      
#      df=NULL
#      for(nl in 1:nt)
#      {
#        lambdal=as.im(list(x=s.grid$x,y=s.grid$y,z=Lambda[,,nl]))
#        df <- rbind(df,as.data.frame(lambdal))
#      }
      
      df=NULL
      for(nl in 1:nt)
      {
        LL = Lambda
        LL[is.na(LL)] = -999
        lambdal=as.im(list(x=s.grid$x,y=s.grid$y,z=LL[,,nl]))
        df <- rbind(df,as.data.frame(lambdal))
      }
      df$value[df$value==-999]=0
      
      samp=sample.int(length(XX),npoints,replace=TRUE,prob=df$value)
      xx <- XX[samp] + runif(npoints, -s.grid$xinc/2, s.grid$xinc/2)
      yy <- YY[samp] + runif(npoints, -s.grid$yinc/2, s.grid$yinc/2)
      
      if (discrete.time==TRUE)
      tt <- floor(TT[samp] + runif(npoints, -t.grid$tinc/2, t.grid$tinc/2))
      else
      tt <- TT[samp] + runif(npoints, -t.grid$tinc/2, t.grid$tinc/2)
        
        
      xyt.init=cbind(x=xx,y=yy,t=tt)

      retain.eq.F <- FALSE
      while(retain.eq.F==FALSE)
      {
         pts <- inpip(pts=cbind(xyt.init[,1],xyt.init[,2]),poly=s.region)
         xyt.init <- xyt.init[pts,]
         
         ptt <- xyt.init[,3]>=t.region[1] & xyt.init[,3]<=t.region[2]
         xyt.init <- xyt.init[ptt,]
         
         npts = dim(xyt.init)[1]
         if (npts == npoints) retain.eq.F <- TRUE
         else
         {
           samp=sample.int(length(XX),npoints-npts,replace=TRUE,prob=df$value)
           
           xx <- XX[samp] + runif(npoints-npts, -s.grid$xinc/2, s.grid$xinc/2)
           yy <- YY[samp] + runif(npoints-npts, -s.grid$yinc/2, s.grid$yinc/2)
           
           if (discrete.time==TRUE)
             tt <- floor(TT[samp] + runif(npoints-npts, -t.grid$tinc/2, t.grid$tinc/2))
           else
             tt <- TT[samp] + runif(npoints-npts, -t.grid$tinc/2, t.grid$tinc/2)
           
           xyt.init1=cbind(x=xx,y=yy,t=tt)
           xyt.init=rbind(xyt.init,xyt.init1)
         }}
      

      index.times <- sort(xyt.init[,3],index.return=TRUE)$ix
    
#      if (anisotropy)
#      {
#      Q=.anisotropy.rotation(theta=ani.pars[1],dzeta=ani.pars[2],inverse=FALSE)
#      resani=Q%*%rbind(x,y)
#      x=resani[1,]
#      y=resani[2,]
#      }
      pattern.interm <- xyt.init[index.times,]
      pattern.interm <- xyt.init

      if (nsim==1)
        {
          pattern <- as.3dpoints(pattern.interm)
          Lambdafin <- Lambda
        }
      else
        {
          pattern[[ni]] <- as.3dpoints(pattern.interm)
          Lambdafin[[ni]] <- Lambda
        }
      ni <- ni+1
    }

  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region,Lambda=Lambdafin,aniso=list(anisotropy=anisotropy,ani=ani.pars),stcov=stcov,covar=covar,grid=list(s.grid,t.grid))))
}










gte <- function(xyt,t.region,t.lambda,dt,kt="epanech",ht,correction="none",approach="simplified"){
  
  verifyclass(xyt,"stpp")
  
  correc <- c("none","isotropic","border","modified.border","translate","setcovf")
  id <- match(correction,correc,nomatch=NA)
  if (any(nbg <- is.na(id))){
    messnbg <- paste("unrecognised correction method:",paste(dQuote(correction[nbg]),collapse=","))
    stop(messnbg,call.=FALSE)
  }
  id <- unique(id)	
  correc2 <- rep(0,6)
  correc2[id] <- 1	
  
  appro <- c("simplified","standardised")
  im <- match(approach,appro,nomatch=NA)
  if (any(nbm <- is.na(im))){
    messnbm <- paste("unrecognised type of estimator:",paste(dQuote(approach[nbm]),collapse=","))
    stop(messnbm,call.=FALSE)
  }
  im <- unique(im)
  appro2 <- rep(0,2)
  appro2[im] <- 1
  
  ker <- c("box","epanech","biweight")
  ik <- match(kt,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(kt[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  ik <- unique(ik)
  ker2 <- rep(0,3)
  ker2[ik] <- 1
  
  dup <- duplicated(data.frame(xyt[,1],xyt[,2],xyt[,3]),fromLast = TRUE)[1]
  if (dup == TRUE){
    messnbd <- paste("spatio-temporal data contain duplicated points")
    warning(messnbd,call.=FALSE)
  }
  
  if (missing(t.region)){
    tr <- range(xyt[,3],na.rm=TRUE)
    tw <- diff(tr)
    t.region <- c(tr[1]-0.01*tw,tr[2]+0.01*tw)
  }
  
  xyt.inside <- .intim(xyt,t.region)

  if (missing(ht)){
    d <- dist(xyt.inside[,3])
    ht <- dpik(d,kernel=kt,range.x=c(min(d),max(d)))
  }
  
  bsupt <- max(t.region)
  binft <- min(t.region)
  W <- sbox(xyt.inside[,1:2], xfrac=0.01, yfrac=0.01)
  a <- diff(range(W[,1]))
  b <- diff(range(W[,2]))
  
  if (missing(dt)) {
    maxt <- (bsupt-binft)/4
    dt <- seq(ht,maxt,len=100)[-1]
  }
  if(dt[1]==0){
    dt <- dt[-1]
  }
  
  kernel <- c(kt=kt,ht=ht)
  gtetheo <- ((a^2)+(b^2))/12
  
  pts <- xyt.inside[,1:2]
  xytimes <- xyt.inside[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  npt <- length(ptsx)
  ndt <- length(dt)
  gtet <- rep(0,ndt)
  
  storage.mode(gtet) <- "double"
  
  if (appro2[1]==1){
    gteout <- .Fortran("gtecore",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),as.double(dt),
                        as.integer(ndt),as.integer(ker2),as.double(ht),(gtet))
    gtet <- gteout[[9]]
    
    dtf <- rep(0,ndt+2)
    dtf[3:(ndt+2)] <- dt
    dt <- dtf 
    
    egte <- rep(0,ndt+2)
    egte[2] <- gtet[1]
    egte[3:(ndt+2)] <- gtet
    
    invisible(return(list(egte=egte,dt=dt,kernel=kernel,gtetheo=gtetheo)))
  } else {
    
    if(missing(t.lambda)){
      misl <- 1
      t.lambda <- rep(npt/(bsupt-binft),npt)
    } else {
      misl <- 0
      if (length(t.lambda)==1){
        t.lambda <- rep(t.lambda,npt)
      }
    }
    
    wrt <- array(0,dim=c(npt,npt))
    wtt <- array(0,dim=c(npt,npt))
    wbit <- array(0,dim=c(npt,ndt))
    wbimodt <- array(0,dim=c(npt,ndt))
    wst <- rep(0,ndt)
    
    # correction="isotropic"
    if(correction=="isotropic"){
      wist <- tedgeRipley(ptst,binft,bsupt)
      wrt <- 1/wist
    }
    
    # correction="translate"
    if(correction=="translate"){
      wtrat <- tedgeTrans(ptst,t.region)
      wtt <- 1/wtrat
    }
    
    #  correction=="border" or "modified border"
    if(any(correction=="border")|any(correction=="modified.border")){
      bj = .bdist.times(xytimes, t.region)
      for(j in 1:ndt) { 
        wbit[,j] <- (bj>dt[j])/sum((bj>dt[j])/t.lambda)
        wbimodt[,j] <- (bj>dt[j])/.eroded.areat(t.region,dt[j])
      }
      wbit[is.na(wbit)] <- 0
      wbimodt[is.na(wbimodt)] <- 0
    }
    
    # correction="setcovf"
    if(correction=="setcovf"){
      wsett <- tsetcovf(dt,ndt,bsupt-binft)
      wst <- 1/wsett
    }
    
    gteout <- .Fortran("gtecoreinh",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),
                        as.double(dt),as.integer(ndt),as.double(t.lambda),as.integer(ker2),
                        as.double(ht),as.double(wrt),as.double(wtt),as.double(wbit),
                        as.double(wbimodt),as.double(wst),as.integer(correc2),(gtet))
    
    gtet <- gteout[[16]]
    
    dtf <- rep(0,ndt+2)
    dtf[3:(ndt+2)] <- dt
    dt <- dtf 
    
    egte <- rep(0,ndt+2)
    egte[2] <- gtet[1]
    egte[3:(ndt+2)] <- gtet
    
    invisible(return(list(egte=egte,dt=dt,kernel=kernel,gtetheo=gtetheo)))
  }
}

tedgeRipley <- function(times,binft,bsupt){
  ntimes <- length(times)
  wrt <- matrix(1,ncol=ntimes,nrow=ntimes)
  
  for(i in 1:ntimes){
    ti <- times[i]
    for(j in 1:ntimes){
      tij <- abs(ti-times[j])
      if (i!=j){
        bsup <- ti+tij
        binf <- ti-tij
        if ((bsup<=bsupt)&(binf>=binft)){
          wrt[i,j] <- 1}
        else {
          wrt[i,j] <- 2
        }
      }
    }
  }
  invisible(return(wrt))
}

tedgeTrans <- function(times,t.region){
  if (missing(t.region)){
    t.region <- range(times)
  }
  
  ntimes <- length(times)
  a <- diff(range(t.region))
  wtt <- matrix(a,ncol=ntimes,nrow=ntimes)
  
  for(i in 1:ntimes){
    for(j in 1:ntimes){
      if (i!=j){
        b <- a-abs(times[i]-times[j])
        wtt[i,j] <- a/b
      }
    }
  }
  invisible(return(wtt))
}

tsetcovf <- function(times,ntimes,longit){
  wst <- rep(0,ntimes)
  for (i in 1:ntimes){
    wst[i] <- longit-times[i]}
  invisible(return(wst=wst))
}

.bdist.times=function(times,t.region){
  if (missing(t.region)){
    t.region <- range(times)
  }
  ntimes <- length(times)
  a <- min(t.region)
  b <- max(t.region)
  
  bj <- NULL
  for(j in 1:ntimes)
    bj <- c(bj,min(c(abs(times[j]-a),abs(times[j]-b))))
  
  invisible(return(bj))
}
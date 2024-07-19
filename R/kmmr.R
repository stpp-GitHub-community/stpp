kmmr <- function(xyt,s.region,s.lambda,ds,ks="epanech",hs,correction="none",approach="simplified"){
  
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
  ik <- match(ks,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(ks[nbk]),collapse=","))
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
  
  options(warn = -1) 
  
  if (missing(s.region)) s.region <- sbox(xyt[,1:2], xfrac=0.01, yfrac=0.01)
  
  xp <- s.region[,1]
  yp <- s.region[,2]
  nedges <- length(xp)
  yp <- yp - min(yp) 
  nxt <- c(2:nedges, 1)
  dx <- xp[nxt] - xp
  ym <- (yp + yp[nxt])/2
  Areaxy <- -sum(dx * ym)
  
  if (Areaxy > 0){
    bsw <- owin(poly = list(x = s.region[,1], y = s.region[,2]))
  }else{
    bsw <- owin(poly = list(x = s.region[,1][length(s.region[,1]):1], y = s.region[,2][length(s.region[,1]):1]))
  }
  
  area <- area(bsw)
  pert <- perimeter(bsw)
  
  ok <- inside.owin(xyt[,1],xyt[,2],w=bsw)
  xyt.inside <- data.frame(x=xyt[,1][ok],y=xyt[,2][ok],t=xyt[,3][ok])
  
  ptsx <- xyt.inside[,1]
  ptsy <- xyt.inside[,2]
  ptst <- xyt.inside[,3]
  
  pxy <- ppp(x=ptsx,y=ptsy,window=bsw)  
  
  if (missing(hs)){
   hs <- bw.stoyan(pxy)
  }
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(0, maxd,len=100)
    ds <- sort(ds)
  }
  if(ds[1]==0){ds <- ds[-1]
  }
  
  kernel <- c(ks=ks,hs=hs)
  kmmrtheo <- 1
  npt <- pxy$n[1]
  nds <- length(ds)
  mummr <- mean(ptst)
  ekmmr <- rep(0,nds)

  storage.mode(ekmmr) <- "double"
  
  if (appro2[1]==1){
    kmmrout <- .Fortran("kmmrcore",as.double(ptsx),as.double(ptsy),as.double(ptst),
                        as.integer(npt),as.double(ds),as.integer(nds),as.integer(ker2)
                        ,as.double(hs),(ekmmr),PACKAGE="msfstpp")
    
    ekmmr <- kmmrout[[9]]/(mummr^2)
    
    invisible(return(list(ekmmr=ekmmr,ds=ds,kernel=kernel,kmmrtheo=kmmrtheo)))
  } else {
    
  if(missing(s.lambda)){
    misl <- 1
    s.lambda <- rep(npt/area,npt)
  } else {
    misl <- 0
  if (length(s.lambda)==1){
    s.lambda <- rep(s.lambda,npt)
  }
    }
    
  wrs <- array(0,dim=c(npt,npt))
  wts <- array(0,dim=c(npt,npt))
  wbi <- array(0,dim=c(npt,nds))
  wbimod <- array(0,dim=c(npt,nds))
  wss <- rep(0,nds)

  # correction="isotropic"
  if(correction=="isotropic"){
    wisot <- edge.Ripley(pxy,pairdist(pxy))
    wrs <- 1/wisot
    }
  
  # correction="translate"
  if(correction=="translate"){
    wtras <- edge.Trans(pxy)
    wts <- 1/wtras
    }
  
  #  correction=="border" or "modified border
  if(any(correction=="border")|any(correction=="modified.border")){
    bi <- bdist.points(pxy)
    for(i in 1:nds) { 
      wbi[,i] <- (bi>ds[i])/sum((bi>ds[i])/s.lambda)
      wbimod[,i] <- (bi>ds[i])/eroded.areas(bsw,ds[i])
     }
    wbi[is.na(wbi)] <- 0
    wbimod[is.na(wbimod)] <- 0
  }
  
  # correction="setcovf"
  if(correction=="setcovf"){
    for (i in 1:nds){
      wss[i] <- area-((pert*ds[i])/pi)
      }
    wss <- 1/wss
  }
  
  options(warn = 0)
  
  kmmrout <- .Fortran("kmmrcoreinh",as.double(ptsx),as.double(ptsy),as.double(ptst),
                     as.integer(npt),as.double(ds),as.integer(nds),as.double(s.lambda),
                     as.integer(ker2),as.double(hs),as.double(wrs),as.double(wts),
                     as.double(wbi),as.double(wbimod),as.double(wss),as.integer(correc2),
                     (ekmmr),PACKAGE="msfstpp")
  
   ekmmr <- kmmrout[[16]]/(mummr^2)
   
   invisible(return(list(ekmmr=ekmmr,ds=ds,kernel=kernel,kmmrtheo=kmmrtheo,s.lambda=s.lambda)))
   }
}

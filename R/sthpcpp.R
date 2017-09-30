sthpcpp <- function(lambp, r, mu, s.region, t.region){

  if (missing(s.region)) s.region <- matrix(c(1,1,0,0,0,1,1,0),ncol=2)
     
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
     
  if (missing(t.region)) {t.region <- c(0,1)}
  
  spp <- rMatClust(lambp, r, mu, win=bsw)
  tpp <- runif(spp$n,t.region[1],t.region[2])
  
  stc1 <- cbind(spp$x,spp$y,tpp)
  stc1 <- as.3dpoints(stc1)

  invisible(return(list(xyt=stc1,s.region=s.region,t.region=t.region)))
}

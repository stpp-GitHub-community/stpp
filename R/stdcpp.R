stdcpp <- function(lambp, a, b, c, mu, s.region, t.region){

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
      
  if (missing(t.region)) t.region <- c(0,1)

  stmc <- NULL
  stPoip <- rpp(lambda=lambp,s.region,t.region)$xyt
  for(i in 1:length(stPoip[,1])){
  PS <- PoiSph(a,b,c,mu=mu,centre=stPoip[i,])
  stmc <- rbind(stmc,PS)}

  ok <- inside.owin(stmc[,1],stmc[,2],w=bsw)
  insw <- data.frame(x=stmc[,1][ok],y=stmc[,2][ok],t=stmc[,3][ok])
  stc3 <- intim(insw,t.region)
  stc3 <- as.3dpoints(stc3)

  invisible(return(list(xyt=stc3,s.region=s.region,t.region=t.region)))
}

PoiSph <- function(a,b,c,mu,centre){
  n <- rpois(1,mu)

  if (a==b & b==c){
    r <- a*((runif(n))^(1/3))

    theta <- 2 * pi * runif(n)
    phi <- acos(2 * runif(n)-1)

    xr <- r * cos(theta)* sin(phi) + centre[1]
    yr <- r * sin(theta)* sin(phi) + centre[2]
    zr <- r * cos(phi) + centre[3]
  }
  else{
  a0 <- a*(runif(n)^(1/3))
  b0 <- b*(runif(n)^(1/3))
  c0 <- c*(runif(n)^(1/3))

  theta <- 2 * pi * runif(n)
  phi <- acos(2 * runif(n)-1)

  x0 <- a0 * cos(theta) * sin(phi)
  y0 <- b0 * sin(theta) * sin(phi)
  z0 <- c0 * cos(phi)

  aa <- 2*pi*runif(1)
  ab <- 2*pi*runif(1)
  ac <- 2*pi*runif(1)

  xr <- (z0*(sin(aa)*sin(ac)+cos(aa)*cos(ac)*sin(ab))-y0*(cos(aa)*sin(ac)-cos(ac)*sin(aa)*sin(ab))+x0*(cos(ac)*cos(ab))) + centre[1]
  yr <- (y0*(cos(aa)*cos(ac)+sin(aa)*sin(ac)*sin(ab))-z0*(cos(ac)*sin(aa)-cos(aa)*sin(ac)*sin(ab))+x0*(cos(ab)*sin(ac))) + centre[2]
  zr <- (z0*(cos(aa)*cos(ab))-x0*(sin(ab))+y0*(cos(ab)*sin(aa))) + centre[3]}

  invisible(return(cbind(xr,yr,zr)))
}

intim <- function(xyt,t.region){
  int <- NULL
  for(i in 1:length(xyt[,1])){
    if (xyt[i,3] > t.region[1] & xyt[i,3] < t.region[2]){
      int <- rbind(int,xyt[i,])}
   }
  
  int <- as.3dpoints(int)

  invisible(return(int))
}

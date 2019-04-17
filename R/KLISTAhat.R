KLISTAhat <-function(xyt, s.region, t.region, dist, times, lambda, correction = "isotropic"){
  
verifyclass(xyt,"stpp")

correc = c("none", "isotropic", "border", "modified.border", "translate")
id <- match(correction, correc, nomatch = NA)
if (any(nbg <- is.na(id))){
  mess <- paste("unrecognised correction method:", paste(dQuote(correction[nbg]), collapse = ","))
stop(mess, call. = FALSE)
}

id = unique(id)
correc2 = rep(0, 5)
correc2[id] = 1

dup <- any(duplicated(data.frame(xyt[,1], xyt[,2], xyt[,3])))
if (dup == TRUE){
  messnbd <- paste("Space-time data contains duplicated points")
  warning(messnbd,call. = FALSE)
}

if (missing(s.region)) 
  s.region <- sbox(xyt[, 1:2], xfrac = 0.01, yfrac = 0.01)
if (missing(t.region)){
  xr = range(xyt[, 3], na.rm = TRUE)
  xw = diff(xr)
  t.region <- c(xr[1] - 0.01 * xw, xr[2] + 0.01 * xw)
}
bsupt <- max(t.region)
binft <- min(t.region)

xp <- s.region[,1]
yp <- s.region[,2]
nedges <- length(xp)
yp <- yp - min(yp) 
nxt <- c(2:nedges, 1)
dx <- xp[nxt] - xp
ym <- (yp + yp[nxt])/2
Areaxy <- -sum(dx * ym)

if (Areaxy > 0){
  bdry <- owin(poly = list(x = s.region[,1], y = s.region[,2]))
}else{
  bdry <- owin(poly = list(x = s.region[,1][length(s.region[,1]):1], y = s.region[,2][length(s.region[,1]):1]))
}

if (missing(dist)){
  rect = as.rectangle(bdry)
  maxd = min(diff(rect$xrange), diff(rect$yrange))/4
  dist = make.even.breaks(maxd, bstep = maxd/512)$r
}
if (missing(times)) {
  maxt = (bsupt - binft)/4
  times = make.even.breaks(maxt, npos = 15)$r
}
dist <- sort(dist)
if (dist[1] == 0) 
  dist = dist[-1]
times <- sort(times)
if (times[1] == 0) 
  times = times[-1]

ok <- inside.owin(xyt[,1],xyt[,2],w=bdry)
xyt.ins <- data.frame(x=xyt[,1][ok],y=xyt[,2][ok],t=xyt[,3][ok])
xyt.in <- .intim(xyt.ins,t.region)

pts <- xyt.in[,1:2]
xytimes <- xyt.in[,3]
ptsx <- pts[, 1]
ptsy <- pts[, 2]
ptst <- xytimes
npt <- length(ptsx)
ndist <- length(dist)
ntimes <- length(times)
area <- areapl(s.region) * (bsupt - binft)
np <- length(s.region[, 1])
polyx <- c(s.region[, 1], s.region[1, 1])
polyy <- c(s.region[, 2], s.region[1, 2])

if (missing(lambda)){
  misl <- 1
  lambda <- rep(npt/area, npt)
}
else misl <- 0
if (length(lambda) == 1) 
  lambda <- rep(lambda, npt)

wbi = array(0, dim = c(npt, ndist, ntimes))
wbimod = array(0, dim = c(npt, ndist, ntimes))
wt = array(0, dim = c(npt, npt))

options(warn = -1)

pppxy = ppp(x = ptsx, y = ptsy, window = bdry)

#  correction=="border" and "modified border"
if (any(correction == "border") | any(correction == "modified.border")) {
  
  bi = bdist.points(pppxy)
  bj = .bdist.times(xytimes, t.region)
  
  for (i in 1 : ndist) {
    for (j in 1 : ntimes) {
      wbi[, i, j] = (bi > dist[i]) * (bj > times[j])/sum((bi > dist[i]) * (bj > times[j])/lambda)
      wbimod[, i, j] = (bi > dist[i]) * (bj > times[j])/(eroded.areas(bdry, dist[i]) * .eroded.areat(t.region, times[j]))
    }
  }
  wbi[is.na(wbi)] <- 0
  wbimod[is.na(wbimod)] <- 0
}

# correction=="translate"
if (any(correction == "translate")) {
  wtt = .overlap.tint(xytimes, t.region)
  wts = edge.Trans(pppxy)
  wt = wtt * wts
  wt = 1/wt
}

options(warn = 0)

all.klista <- sapply(seq(1,npt), function(i) .i.KLISTAhat(i,id,ptsx,ptsy,ptst,npt,polyx,polyy,np,dist,ndist,times,ntimes,bsupt,binft,lambda,wbi,wbimod,wt,correc2,correc,area),simplify="array")

correction = correc[id]
klistatheo <- matrix(0,ncol=length(times),nrow=length(dist))
for(i in 1:length(dist)) klistatheo[i,] <- pi*(dist[i]^2)*times

invisible(return(list(list.KLISTA = all.klista, klistatheo = klistatheo, dist = dist, times = times, correction = correction)))
}

.i.KLISTAhat <- function(i,id,ptsx,ptsy,ptst,npt,polyx,polyy,np,dist,ndist,times,ntimes,bsupt,binft,lambda,wbi,wbimod,wt,correc2,correc,area){
  
  xi <- ptsx[i]
  yi <- ptsy[i]
  ti <- ptst[i]
  klistahat <- array(0, dim = c(ndist, ntimes, 5))
  
  storage.mode(klistahat) <- "double"
  
  kflista <- .Fortran("klistafunction", as.integer(i),
                    as.double(xi), as.double(yi), 
                    as.double(ti), as.double(ptsx), 
                    as.double(ptsy), as.double(ptst), 
                    as.integer(npt), as.double(polyx), 
                    as.double(polyy), as.integer(np), 
                    as.double(dist), as.integer(ndist), 
                    as.double(times), as.integer(ntimes),
                    as.double(bsupt), as.double(binft), 
                    as.double(lambda), (klistahat),
                    as.double(wbi), as.double(wbimod), 
                    as.double(wt), as.integer(correc2))
  
  klistahat <- kflista[[19]]
  
  klistahat[, , c(1, 2, 5)] = klistahat[, , c(1, 2, 5)]/area
  klistahat <- (npt-1)*klistahat
 
   if (length(id) == 1) 
     iKLISTA = as.array(klistahat[, , id])
   else {
     iKLISTA = list()
     for (i in 1:length(id)) iKLISTA[[i]] = klistahat[, , id[i]]
     names(iKLISTA) = correc[id]
   }
   
  return(iKLISTA)
 }
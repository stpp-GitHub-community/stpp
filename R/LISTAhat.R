LISTAhat <-function(xyt, s.region, t.region, dist, times, lambda, ks = "box", hs, kt = "box", ht, correction = "isotropic"){
  
verifyclass(xyt,"stpp")

correc = c("none", "isotropic", "border", "modified.border", "translate")
id <- match(correction, correc, nomatch = NA)
if (any(nbg <- is.na(id))){
  mess <- paste("unrecognised correction method:", paste(dQuote(correction[nbg]), collapse = ", "))
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

frac = 1
if (missing(hs)) {
  d = dist(pts)
  if (ks == "gaussian") 
    hs = dpik(d, kernel = "normal", range.x = c(min(d), max(d)/frac))
  else hs = dpik(d, kernel = ks, range.x = c(min(d), max(d)/frac))
}
if (missing(ht)) {
  d = dist(ptst)
  if (kt == "gaussian") 
    ht = dpik(d, kernel = "normal", range.x = c(min(d), max(d)/frac))
  else ht = dpik(d, kernel = kt, range.x = c(min(d), max(d)/frac))
}
kernel = c(ks = ks, hs = hs, kt = kt, ht = ht)
if (ks == "box") 
  ks = 1
else if (ks == "epanech") 
  ks = 2
else if (ks == "gaussian") 
  ks = 3
else if (ks == "biweight") 
  ks = 4
if (kt == "box") 
  kt = 1
else if (kt == "epanech") 
  kt = 2
else if (kt == "gaussian") 
  kt = 3
else if (kt == "biweight") 
  kt = 4

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

all.lista <- sapply(seq(1,npt), function(i) .i.LISTAhat(i,id,ptsx,ptsy,ptst,npt,polyx,polyy,np,dist,ndist,times,ntimes,bsupt,binft,lambda,ks,kt,hs,ht,wbi,wbimod,wt,correc2,correc,area),simplify="array")

correction = correc[id]
listatheo = array(1+(1/(npt-1)), dim = c(ndist, ntimes))

invisible(return(list(list.LISTA = all.lista, listatheo = listatheo, dist = dist, times = times, kernel = kernel, correction = correction)))
}

.i.LISTAhat <- function(i,id,ptsx,ptsy,ptst,npt,polyx,polyy,np,dist,ndist,times,ntimes,bsupt,binft,lambda,ks,kt,hs,ht,wbi,wbimod,wt,correc2,correc,area){
  
  xi <- ptsx[i]
  yi <- ptsy[i]
  ti <- ptst[i]
  listahat <- array(0, dim = c(ndist, ntimes, 5))
  
  storage.mode(listahat) <- "double"
  
  pcflista <- .Fortran("listafunction", as.integer(i),
                    as.double(xi), as.double(yi), 
                    as.double(ti), as.double(ptsx), 
                    as.double(ptsy), as.double(ptst), 
                    as.integer(npt), as.double(polyx), 
                    as.double(polyy), as.integer(np), 
                    as.double(dist), as.integer(ndist), 
                    as.double(times), as.integer(ntimes),
                    as.double(bsupt), as.double(binft), 
                    as.double(lambda), as.integer(ks), 
                    as.integer(kt), as.double(hs), 
                    as.double(ht), (listahat), 
                    as.double(wbi), as.double(wbimod), 
                    as.double(wt), as.integer(correc2))
  
  listahat <- pcflista[[23]]
  
  listahat[, , c(1, 2, 5)] = listahat[, , c(1, 2, 5)]/area
  listahat <- (npt-1)*listahat/(4 * pi * dist)
  
  if (length(id) == 1) 
    iLISTA = as.array(listahat[, , id])
  else {
    iLISTA = list()
    for (i in 1:length(id)) iLISTA[[i]] = listahat[, , id[i]]
    names(iLISTA) = correc[id]
  }
  
 return(iLISTA)
}
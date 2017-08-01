plot.stpp <- function(xyt, s.region=NULL, t.region=NULL, mark=FALSE, mark.cexmin=0.4, mark.cexmax=1.2, mark.col=1, ...){
  
  verifyclass(xyt,"stpp")
  
  if (missing(s.region)) s.region <- sbox(xyt[,1:2], xfrac=0.01, yfrac=0.01)
  if (missing(t.region)) 
  {
    xr = range(xyt[,3], na.rm = TRUE)
    xw = diff(xr)
    t.region <- c(xr[1]- 0.01 * xw, xr[2] + 0.01 * xw)
  } 
  
  nedges <- length(s.region[,1])
  s.region[,2] <- s.region[,2] - min(s.region[,2])
  nxt <- c(2:nedges, 1)
  dx <- s.region[,1][nxt] - s.region[,1]
  ym <- (s.region[,2] + s.region[,2][nxt])/2
  Areaxy <- -sum(dx * ym)
  
  if (Areaxy > 0){
    bdry = owin(poly = list(x = s.region[,1], y = s.region[,2]))}
  else
    bdry = owin(poly = list(x = s.region[,1][length(s.region[,1]):1], y = s.region[,2][length(s.region[,1]):1]))
  
if (inherits(xyt,"stpp")==TRUE) 
	{ 
	if (mark==FALSE)
	  {
	  par(mfrow=c(1,2),pty="s")
	  if (is.null(s.region))	
	  plot(xyt[,1:2],main="xy-locations",...)
	  else
		{
		  polymap(s.region,xlab="x",ylab="y")
		  points(xyt[,1:2],...)
		  title("xy-locations")	 
		}
	  plot(sort(xyt[,3]),1:length(xyt[,3]),type="l",xlab="t",ylab="",main="cumulative number",las=1,xlim=t.region)
	   }
	if (mark==TRUE)
	 {
	  l=dim(x)[1]
	  CEX=seq(mark.cexmin,mark.cexmax,length=l)
        if (mark.col==0)
	     {
          par(mfrow=c(1,1),pty="s")
		  if (is.null(s.region))	
	  	   plot(x[,1:2],cex=CEX,...)
 		  else
			{
			  polymap(s.region,xlab="x",ylab="y")
			  points(x[,1:2],cex=CEX,...)	 
			}
          }
        else 
         {  
	       if (mark.col=="black" | mark.col==1)	
	           COL=grey((l:1)/l)
           if (mark.col=="red" | mark.col==2)	
	           COL=rgb(l:0, 0, 0, maxColorValue = l)
      	   if (mark.col=="green" | mark.col==3)	
        	   COL=rgb(0, l:0, 0, maxColorValue = l)
	       if (mark.col=="blue" | mark.col==4)	
	           COL=rgb(0, 0, l:0, maxColorValue = l)
	       par(mfrow=c(1,1),pty="s")
    	  if (is.null(s.region))	
	     	plot(x[,1:2],col=COL,cex=CEX,...)
 	      else
		  {
		   polymap(s.region,xlab="x",ylab="y")
		      points(x[,1:2],col=COL,cex=CEX,...)	 
		  }
         }
	 }	
}
  else
  {warning("x should be of class 'stpp'")}
}

getS3method("plot", "stpp", optional = FALSE)


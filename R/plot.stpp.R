plot.stpp <- function(x, s.region=NULL, t.region=NULL, scatter=FALSE, mark=FALSE, mark.cexmin=0.4, mark.cexmax=1.2, mark.col=1, ...){
  
  verifyclass(x,"stpp")
  
if (inherits(x,"stpp")==TRUE) 
	{ 
  
  if (scatter==TRUE){
    par(mfrow=c(1,1),pty="s")
    scatter3D(x[,1],x[,2],x[,3],zlab="\n t",main="xyt-locations",...)
  }
  else{

	if (mark==FALSE)
	  {
	  par(mfrow=c(1,2),pty="s")
	  if (is.null(s.region))	
	  plot(x[,1:2],main="xy-locations",...)
	  else
		{ polymap(s.region,xlab="x",ylab="y")
		  points(x[,1:2],...)
		  title("xy-locations")	 
		}
	  plot(sort(x[,3]),1:length(x[,3]),type="l",xlab="t",ylab="",main="cumulative number",las=1,xlim=t.region)
	}
  
	if (mark==TRUE)
	 {
	  snorm=apply(x[,1:2],MARGIN=1,FUN=norm,type="2")
	  t=x[,3]
	  l=dim(x)[1]
	  CEX=seq(mark.cexmin,mark.cexmax,length=l)
        if (mark.col==0)
	     {
          par(mfrow=c(1,2),pty="s")
		  if (is.null(s.region))	
	  	   plot(x[,1:2],cex=CEX,main="Time mark",...)
 		  else
			{
			  polymap(s.region,xlab="x",ylab="y")
			  points(x[,1:2],cex=CEX,...)
			  title("Time mark")
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
	       par(mfrow=c(1,2),pty="s")
    	  if (is.null(s.region))	
	     	plot(x[,1:2],col=COL,cex=CEX,main="Time mark",...)
 	      else
		  {
		   polymap(s.region,xlab="x",ylab="y")
		   points(x[,1:2],col=COL,cex=CEX,...)	 
		  title("Time mark")
		  }
         }
	  plot(t,snorm,type="h",ylab="||(x,y)||",main="Space mark",...)
	  }	
   }
  }
  else
  {warning("x should be of class 'stpp'")}
}
  

getS3method("plot", "stpp", optional = FALSE)


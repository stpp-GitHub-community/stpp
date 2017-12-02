plot.stpp <- function(x, s.region=NULL, t.region=NULL, style="generic", type="projection", mark=NULL, mark.cexmin=0.4, mark.cexmax=1.2, mark.col=1, ...){
  
  verifyclass(x,"stpp")
  
  if(!is.null(mark)){
    if (mark==TRUE){
      type <- "mark"
      style <- "generic"
      }
    if (mark==FALSE){
      type <- "projection"
      style <- "generic"
    }
  }
  
  if (type=="scatter"){
    if(style=="generic"){
      par(mfrow=c(1,1),pty="s")
      scatter3D(x[,1],x[,2],x[,3],zlab="\n t",main="xyt-locations",...)
    }
    
    if(style=="elegant"){
      par(mfrow=c(1,1),pty="s")
      scatter3D(x[,1],x[,2],x[,3],zlab="\n t",main="xyt-locations",facets=FALSE,curtain=FALSE,col.panel="#EBEBEB",col.grid="F8F8F8",bty="g",border=NA,...)
      }
    }

  if (type=="projection"){
      if(style=="generic"){
        par(mfrow=c(1,2),pty="s")
        if (is.null(s.region)){
          plot(x[,1:2],main="xy-locations",...)
	     }else{
	        polymap(s.region,xlab="x",ylab="y")
		      points(x[,1:2],...)
		      title("xy-locations")}
	        plot(sort(x[,3]),seq(1,length(x[,3])),type="l",xlab="t",ylab="",main="cumulative number",las=1,xlim=t.region)
      }
    
      if(style=="elegant"){
        sd <- data.frame(x=x[,1],y=x[,2])
        td <- data.frame(t=sort(x[,3]),y=seq(1,length(x[,3])))
        if (is.null(s.region)){
          sp <- ggplot(sd,aes_string("x","y"))+geom_point()+theme(text=element_text(...))+ggtitle("xy-locations")
        }else{
          SW <- data.frame(x=c(s.region[,1],s.region[1,1]),y=c(s.region[,2],s.region[1,2]))
          sp <- ggplot(sd,aes_string("x","y"))+geom_path(data=SW,aes_string("x","y"))+geom_point()+theme(text=element_text(...))+ggtitle("xy-locations")
        }
          tp <- ggplot(td,aes_string("t","y"))+geom_line()+labs(y="")+theme(text=element_text(...))+ggtitle("cumulative number")
          grid.arrange(sp,tp,ncol=2,nrow=1)
      }
    }
  
	if (type=="mark"){
	  if(is.null(mark)){
	  sn=apply(x[,1:2],MARGIN=1,FUN=norm,type="2")
	  t=x[,3]}
	  if(style=="generic"){
	  l=dim(x)[1]
	  CEX=seq(mark.cexmin,mark.cexmax,length=l)
    if (mark.col==0){
      if (!is.null(mark)){
      par(mfrow=c(1,1),pty="s")
      }else{
      par(mfrow=c(1,2),pty="s")
      }
      if (is.null(s.region)){	
	  	   plot(x[,1:2],cex=CEX,...)
        if(is.null(mark)){
          title("Time mark")
        }
        }else{
          polymap(s.region,xlab="x",ylab="y")
			    points(x[,1:2],cex=CEX,...)
			    if(is.null(mark)){
			    title("Time mark")
			      }
			    }
      }else{
	      if (mark.col=="black" | mark.col==1){COL=grey((l:1)/l)}
        if (mark.col=="red" | mark.col==2){COL=rgb(l:0, 0, 0, maxColorValue = l)}
      	if (mark.col=="green" | mark.col==3){COL=rgb(0, l:0, 0, maxColorValue = l)}
	      if (mark.col=="blue" | mark.col==4){COL=rgb(0, 0, l:0, maxColorValue = l)}
        if (!is.null(mark)){
          par(mfrow=c(1,1),pty="s")
        }else{
          par(mfrow=c(1,2),pty="s")
        }
        if (is.null(s.region)){
			        plot(x[,1:2],col=COL,cex=CEX,...)
          if(is.null(mark)){
            title("Time mark")
          }
 	      }else{
 	        polymap(s.region,xlab="x",ylab="y")
 	        points(x[,1:2],col=COL,cex=CEX,...)
 	        if(is.null(mark)){
 	          title("Time mark")
 	        }
 	      }
      }
	  if(is.null(mark)){
	        plot(t,sn,type="h",ylab="||(x,y)||",main="Space mark",...)
	        }}
	  if (style=="elegant"){
	    sdm <- data.frame(x=x[,1],y=x[,2],t=x[,3])
	    if (is.null(s.region)){
	    spm <- ggplot(sdm,aes_string("x","y"))+geom_point(aes(size=t),shape=1,show.legend=FALSE)+scale_shape(solid=FALSE)+theme(text=element_text(...))+ggtitle("Time marks")
	    }else{
	      SW <- data.frame(x=c(s.region[,1],s.region[1,1]),y=c(s.region[,2],s.region[1,2]))
	      spm <- ggplot(sdm,aes_string("x","y"))+geom_path(data=SW,aes_string("x","y"))+geom_point(aes(size=t),shape=1,show.legend=FALSE)+theme(text=element_text(...))+ggtitle("Time marks")
	    }
	    tdm <- data.frame(t=t,sn=sn)
	    tpm <- ggplot(tdm,aes_string("t","sn"))+geom_bar(stat="identity",colour="black")+labs(y="||(x,y)||")+theme(text=element_text(...))+ggtitle("Space marks")
	    grid.arrange(spm,tpm,ncol=2,nrow=1)
	  }
	}
}

getS3method("plot", "stpp", optional = FALSE)
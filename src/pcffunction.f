C
C     E. Gabriel, August 2012
C
C     This function provides an edge corrected estimate
C     of the space-time pair correlation function.
C

      subroutine pcffunction(x,y,txy,n,xp,yp,np,s,ns,t,nt,
     +  bsupt,binft,lambda,ks,kt,hs,ht,pcfhat,
     +  wbi,wbimod,wt,correc)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     x,y,txy: coordinates and times of the point process of length n
c     xp,yp: coordinates of the np points defining the polygonal
c            region
c     s: vector of the ns distances at which to calculate the K
c        function,
c     t: vector of the nt times at which to calculate the K function,
c     bint, bsupt: lower and upper boundaries of the time domain,
c     pcfhat: zero matrix of dimension ns x nt.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      integer n,ns,nt,np,iu,iv,ks,kt,correc(5)
      double precision pcfhat(ns,nt,5), two, hs, ht, lambda(n)
      double precision wbi(n,ns,nt), wbimod(n,ns,nt), wt(n,n)
      dimension x(n),y(n),txy(n),xp(np+1),yp(np+1),s(ns),t(nt)
      double precision binf, binft, bsup, bsupt, ti, tij
      double precision vij, wij
      double precision kern, kerns, kernt

      two=2d0
      kerns=0d0
      kernt=0d0

      do iu=1,ns
      do iv=1,nt
        do i=1,n
        xi=x(i)
        yi=y(i)
        ti=txy(i)
        do j=1,n
        if (j.ne.i) then
            hij=dsqrt((xi-x(j))*(xi-x(j)) + (yi-y(j))*(yi-y(j)))
            tij=dabs(ti-txy(j))
            if (ks.eq.1) then
                kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks.eq.2) then
                    kerns=ekernel((s(iu)-hij)/hs,hs)
                    else if (ks.eq.3) then
                          kerns=gausskernel((s(iu)-hij)/hs,hs)
                            else if (ks.eq.4) then
                            kerns=qkernel((s(iu)-hij)/hs,hs)
            end if
            if (kt.eq.1) then
                kernt=boxkernel((t(iv)-tij)/ht,ht)
                else if (kt.eq.2) then
                    kernt=ekernel((t(iv)-tij)/ht,ht)
                    else if (kt.eq.3) then
                        kernt=gausskernel((t(iv)-tij)/ht,ht)
                            else if (kt.eq.4) then
                            kernt=qkernel((t(iv)-tij)/ht,ht)
            end if
            kern=kerns*kernt
            if (kern.ne.0) then
c isotropic
            if(correc(2).eq.1) then
                    bsup=ti+tij
                    binf=ti-tij
                    if ((bsup.le.bsupt).and.(binf.ge.binft)) then
                      vij=1d0
                      else
                      vij=two
                    end if
                    wij=weight(xi,yi,hij,xp,yp,np)
                    wij=kern*vij*wij/(lambda(i)*lambda(j))
                    pcfhat(iu,iv,2)=pcfhat(iu,iv,2)+wij
            end if
c None
            if (correc(1).eq.1) then
                    wij=kern/(lambda(i)*lambda(j))
                    pcfhat(iu,iv,1)=pcfhat(iu,iv,1)+wij
            end if
c border
            if (correc(3).eq.1) then
                  wij=wbi(i,iu,iv)
                  wij=kern*wij/(lambda(i)*lambda(j))
              pcfhat(iu,iv,3)=pcfhat(iu,iv,3)+wij
            end if
c modified border
            if (correc(4).eq.1) then
                    wij=wbimod(i,iu,iv)
                    wij=kern*wij/(lambda(i)*lambda(j))
                    pcfhat(iu,iv,4)=pcfhat(iu,iv,4)+wij
            end if
c translate
            if (correc(5).eq.1) then
              wij=wt(i,j)
              wij=kern*wij/(lambda(i)*lambda(j))
              pcfhat(iu,iv,5)=pcfhat(iu,iv,5)+wij
            end if
            end if
        end if
        end do
        end do
        end do
        end do


      return

      end



c--------------------------------------------------------------------
c
c     boxkernel
c
c--------------------------------------------------------------------

       function boxkernel(x,h)

       implicit double precision (a-h,o-z)

       double precision x, h

       if (dabs(x).le.1) then
           boxkernel=1d0/2d0
       else
           boxkernel=0d0
       end if
       boxkernel=boxkernel/h

       return
       end

c--------------------------------------------------------------------
c
c     Epanechnikov kernel
c
c--------------------------------------------------------------------

       function ekernel(x,h)

       implicit double precision (a-h,o-z)

       double precision x

       if (dabs(x).le.1) then
           ekernel=(3d0/4d0)*(1-x**2)
       else
           ekernel=0d0
       end if
       ekernel=ekernel/h

       return
       end

c--------------------------------------------------------------------
c
c     Gaussian kernel
c
c--------------------------------------------------------------------

       function gausskernel(x,h)

       implicit double precision (a-h,o-z)

       double precision x
       pi=3.14159265d0
       
       gausskernel=exp(-(x**2)/2d0)/sqrt(pi*2d0)
       gausskernel=gausskernel/h

       return
       end


c--------------------------------------------------------------------
c
c     quartic (biweight) kernel
c
c--------------------------------------------------------------------

       function qkernel(x,h)

       implicit double precision (a-h,o-z)

       double precision x, h

       if (dabs(x).le.1) then
           qkernel=(15d0/16d0)*(1-x**2)**2
       else
           qkernel=0d0
       end if
       qkernel=qkernel/h

       return
       end

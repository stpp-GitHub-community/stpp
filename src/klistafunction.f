C
C     Francisco J., Rodriguez-Cortes, March 2019
C
C     This function provides an edge corrected estimate
C     of the space-time K LISTA.
C

      subroutine klistafunction(i,xi,yi,ti,x,y,txy,n,xp,yp,np,s,ns,t,nt,
     +   bsupt,binft,lambda,klistahat,wbi,wbimod,wt,correc)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     x,y,txy: coordinates and times of the point process of length n
c     xp,yp: coordinates of the np points defining the polygonal
c            region
c     s: vector of the ns distances at which to calculate the ith LISTA
c        function,
c     t: vector of the nt times at which to calculate the ith LISTA 
c        function,
c     bint, bsupt: lower and upper boundaries of the time domain,
c     klistahat: zero matrix of dimension ns x nt.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      integer n,ns,nt,np,iu,iv,correc(5),i
      double precision klistahat(ns,nt,5),two,lambda(n),xi,yi,ti
      double precision wbi(n,ns,nt), wbimod(n,ns,nt), wt(n,n)
      dimension x(n),y(n),txy(n),xp(np+1),yp(np+1),s(ns),t(nt)
      double precision binf, binft, bsup, bsupt, tij
      double precision vij, wij

      two=2d0

      do j=1,n
       do iu=1,ns
       do iv=1,nt
        if (j.ne.i) then
            hij=dsqrt((xi-x(j))*(xi-x(j)) + (yi-y(j))*(yi-y(j)))
            tij=dabs(ti-txy(j))
            if ((tij.le.t(iv)).and.(hij.le.s(iu))) then
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
                    wij=vij*wij/(lambda(i)*lambda(j))
                    klistahat(iu,iv,2)=klistahat(iu,iv,2)+wij
             end if
c none
             if (correc(1).eq.1) then
                    vij=1d0
                    wij=vij/(lambda(i)*lambda(j))
                    klistahat(iu,iv,1)=klistahat(iu,iv,1)+wij
             end if
c border
             if (correc(3).eq.1) then
                   wij=wbi(i,iu,iv)
                    wij=wij/(lambda(i)*lambda(j))
                    klistahat(iu,iv,3)=klistahat(iu,iv,3)+wij
             end if
c modified border
             if (correc(4).eq.1) then
                   wij=wbimod(i,iu,iv)
                   wij=wij/(lambda(i)*lambda(j))
                   klistahat(iu,iv,4)=klistahat(iu,iv,4)+wij
             end if
c translate
             if (correc(5).eq.1) then
                 wij=wt(i,j)/(lambda(i)*lambda(j))
                 klistahat(iu,iv,5)=klistahat(iu,iv,5)+wij
             end if
             
        end if     
        end if
        end do
        end do
        end do

      return

      end

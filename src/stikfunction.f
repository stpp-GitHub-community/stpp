C
C     E. Gabriel, september 2005
C
C     This function provides an edge corrected estimate
C     of the space-time inhomogeneous K function.
C

      subroutine stikfunction(x,y,txy,n,xp,yp,np,s,ns,t,nt,
     +     bsupt,binft,lambda,infd,hkhat,wbi,wbimod,wt,correc)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     x,y,txy: coordinates and times of the point process of length n
c     xp,yp: coordinates of the np points defining the polygonal
c            region
c     s: vector of the ns distances at which to calculate the K
c        function,
c     t: vector of the nt times at which to calculate the K function,
c     lambda: vectors of the space-time intensity functions evaluated
c             to (x,y,txy),
c     bint, bsupt: lower and upper boundaries of the time domain,
c     hkhat: zero matrix of dimension ns x nt.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision(a-h,o-z)

      integer n,ns,nt,np,infd,iu,iv,nv,i,j,correc(5), nev(nt)
      double precision hkhat(ns,nt,5), lambda(n), two
      double precision wbi(n,ns,nt), wbimod(n,ns,nt), wt(n,n)
      dimension x(n),y(n),txy(n),xp(np+1),yp(np+1),s(ns),t(nt)
      double precision binf, binft, bsup, bsupt, ti, tij
      double precision vij, wij

      two=2d0

       if (infd.eq.1) then
        do iv=1,nt
            nv=0
            do i=1,n
            if (txy(i).lt.(bsupt-t(iv))) then
                nv=nv+1
            end if
          end do
             nev(iv)=nv
         end do

      do iu=1,ns
      do iv=1,nt
        nv=nev(iv)
        do i=2,nv
        i1=i-1
        xi=x(i)
        yi=y(i)
        ti=txy(i)
        do j=1,i1
            hij=dsqrt((xi-x(j))*(xi-x(j)) + (yi-y(j))*(yi-y(j)))
            tij=dabs(ti-txy(j))
            if ((tij.le.t(iv)).and.(hij.le.s(iu))) then
c isotropic
                wij=weight(xi,yi,hij,xp,yp,np)
                wij=wij/(lambda(i)*lambda(j))
                hkhat(iu,iv,2)=hkhat(iu,iv,2)+wij
         end if
        end do
        end do
        hkhat(iu,iv,2)=hkhat(iu,iv,2)*(n*1d0/nv)
        end do
        end do
       end if

      if (infd.eq.0) then
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
                    hkhat(iu,iv,2)=hkhat(iu,iv,2)+wij
             end if
c none
             if (correc(1).eq.1) then
                    vij=1d0
                    wij=vij/(lambda(i)*lambda(j))
                    hkhat(iu,iv,1)=hkhat(iu,iv,1)+wij
             end if
c border
             if (correc(3).eq.1) then
                   wij=wbi(i,iu,iv)
                    wij=wij/(lambda(i)*lambda(j))
                    hkhat(iu,iv,3)=hkhat(iu,iv,3)+wij
             end if
c modified border
             if (correc(4).eq.1) then
                   wij=wbimod(i,iu,iv)
                   wij=wij/(lambda(i)*lambda(j))
                   hkhat(iu,iv,4)=hkhat(iu,iv,4)+wij
             end if
c translate
             if (correc(5).eq.1) then
                 wij=wt(i,j)/(lambda(i)*lambda(j))
                   hkhat(iu,iv,5)=hkhat(iu,iv,5)+wij
             end if
            end if
        end if
        end do
        end do
        end do
        end do
       end if



      return

      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     functions called by stikfunction:
c     -----------------------------------------
c
c     * iplace
c     * weight
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c--------------------------------------------------------------------
c
c     iplace
c
c--------------------------------------------------------------------

      function iplace(s,ns,t)
c
c which of the variable width bins s is t in?
c
      implicit DOUBLE PRECISION (a-h,o-z)

      dimension s(ns)

      do ib=1,ns
        if(s(ib).ge.t)then
          iplace=ib
          return
        end if
      end do
c
c if it is outside the range of s
c
      iplace=ns+1

      return
      end

c--------------------------------------------------------------------
c
c     weight
c
c--------------------------------------------------------------------

      function weight(x,y,r,xp,yp,np)
c
c find the weight for the point at x,y, radius r
c
      implicit DOUBLE PRECISION (a-h,o-z)

c      include 'bounds.cmn'

      dimension xp(np+1),yp(np+1)

      weight=cncvwt(x,y,r,xp,yp,np)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     function called by weight:
c     --------------------------
c
c     * cncvwt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c--------------------------------------------------------------------
c
c     cncvwt
c
c--------------------------------------------------------------------

      function cncvwt(x,y,r,xp,yp,np)
c
c     compute the weight given to a point at x,y according to how much
c     of a circle of radius r is inside the bounding polygon
c
c
      implicit DOUBLE PRECISION (a-h,o-z)
c      include 'bounds.cmn'
      dimension xp(np+1),yp(np+1)
      parameter(pi=3.141592654d0)
c     store circle/poly intersections here
      parameter(maxcrs=40)
      dimension cross(maxcrs+1)
      parameter(tiny=1.0e-7)
c     set count of crossing points to zero
      ncross = 0
      angtes = 0
c     first loop over the boundary and find the crossing points
      do ic=1,np
c     work with the trial point at origin
         x1=xp(ic)-x
         y1=yp(ic)-y
         x2=xp(ic+1)-x
         y2=yp(ic+1)-y

         cx=x2-x1
         cy=y2-y1

c     these are the coefficients of the quadratic giving the
c     intercept of line and circle.
         a=cx*cx+cy*cy
         b=2*(x1*cx+y1*cy)
         c=x1*x1+y1*y1-r*r

c     find out if real solutions exist...
         b2m4ac=b*b-4*a*c

c     ... and if they do, find them.
         if (b2m4ac.ge.0) then
            t1=(-b+sqrt(b2m4ac))/(2*a)
            t2=(-b-sqrt(b2m4ac))/(2*a)

c     see if the solutions lie in the line segments
            if ((t1.gt.tiny).and.(t1-1.0.le.tiny)) then
               ncross=ncross+1
c     find the angle to this point on the circle
               ctemp=atan2(y1+t1*cy,x1+t1*cx)
               if(ctemp.lt.0)ctemp=2*pi+ctemp
               cross(ncross)=ctemp
c     check crossing of circle with vertex
            else if (abs(t1).le.tiny) then
c     compare this polygon segment's direction with that of the
c     previous one
               nprev = (mod((ic+ (np-2)),np)+1)
               x0 = xp(nprev) - x
               y0 = yp(nprev) - y
               idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
               idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c     see if the polygon passes through the circle here
               if ((idp1-idp2).ne.1 .and.
     +              abs(idp1+idp2).ne.2) then
                  ncross = ncross + 1
                  ctemp = atan2(y1+t1*cy,x1+t1*cx)
                  if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
                  cross(ncross) = ctemp
               end if
            end if

            if ((t2.gt.tiny).and.(t2-1.0.lt.tiny)) then
               ncross=ncross+1
               ctemp=atan2(y1+t2*cy,x1+t2*cx)
               if(ctemp.lt.0)ctemp=2*pi+ctemp
               cross(ncross)=ctemp
c     check crossing of circle with vertex
            else if (abs(t2).le.tiny)then
c     compare this polygon segment's direction with that of the
c     previous one
               nprev = (mod((ic+ (np-2)),np)+1)
               x0 = xp(nprev) - x
               y0 = yp(nprev) - y
               idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
               idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c     see if the polygon passes through the circle here
               if ((idp1-idp2).ne.1 .and.
     +              abs(idp1+idp2).ne.2) then
                  ncross = ncross + 1
                  ctemp = atan2(y1+t2*cy,x1+t2*cx)
                  if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
                  cross(ncross) = ctemp
               end if
            end if
         end if
      end do

c     now we have all the crossing point angles stored in
c     cross(1:ncross)

c     if ncross = 0 then the total angle within the poly is 2*pi
c     unless the circle is large and spans the polygon. this should
c     be checked beforehand so it's okay to assume 2*pi here.

      if (ncross.eq.0) then
         totang=2*pi
      else

c     sort into ascending order
         call sort2(cross,ncross)

c     fix the ncross+1'th element to be the first plus 2pi so that
c     the list is circular...
         cross(ncross+1)=cross(1)+2*pi

c     check that the number of crossings is even - if not then error.
         if (mod(ncross,2).ne.0) then
            cncvwt=-1
            return
         end if
c     now find a nice spot to do the point-in-poly search
         sepmax=0.0
         icm=0

         do ic=1,ncross
            if (cross(ic+1)-cross(ic).gt.sepmax) then
               sepmax=cross(ic+1)-cross(ic)
               icm=ic
            end if
         end do

c     icm is now the index of the crossing with the largest gap
c     between it and the next crossing point

c     test for point in poly of the point on the circle between these
c     points angtes=(cross(icm)+cross(icm+1))/2.

         xtest=x+r*cos(angtes)
         ytest=y+r*sin(angtes)

c     find out if test point is in the polygon boundary
         linpol=ipippa(xtest,ytest,xp,yp,np)

c     find the total angle between (odd-even) crossings
c     (i.e. 1-2 + 3-4 + ...)
        totang = 0.
        do ic=1,ncross-1,2
           totang = totang + (cross(ic+1)-cross(ic))
        end do

c     If the point we tested for p-i-p was on an odd-even
c     section and was in the poly, then totang is the amount of circle
c     inside the polygon. if the point was outside the polygon, then
c     we need to subtract totang from 2*pi radians to get the angle
c     inside the polygon. conversely, if the point tested was between
c     even-odd crossings and outside the polygon, then totang is the
c     angle we want, and if inside the polygon then again we have to
c     do 2*pi-totang

        if ( (((mod(icm,2).eq.1).and.(linpol.eq.0))  .or.
     &       ((mod(icm,2).eq.0).and.(linpol.eq.1)) ) ) then
           totang = 2*pi-totang
        end if

      end if
c     now totang is the angle contained in the polygon

c     weight is proportion of total angle in the poly
      cncvwt = (2*pi)/(totang)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     functions called by cncvwt:
c     ---------------------------
c
c     * isig8
c     * sort2
c     * ipippa
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c--------------------------------------------------------------------
c
c     isig8
c
c--------------------------------------------------------------------

      integer function isig8(value,tiny)

c     return the sign (+1,0,-1) of a value

      DOUBLE PRECISION tiny,value
      if (value.gt.tiny) then
         isig8 = 1
      else if (value.lt.-tiny) then
         isig8 = -1
      else
         isig8 = 0
      end if
      return
      end

c--------------------------------------------------------------------
c
c     sort2
c
c--------------------------------------------------------------------

      subroutine sort2(x,n)
c
c     shellsort algorithm
c     n     : number of elements to be sorted
c     x     : on enter an array of dimension at least n containing
c             real numbers
c             on output first n elements of x are sorted from smallest
c             to largest
c
      implicit DOUBLE PRECISION (a-h,o-z)
      dimension x(n)
      i=1
    1 i=i+1
      if (i.le.n) goto 1
      m=i-1
    2 m=m/2
      if (m.eq.0) return
      k=n-m
      do 4 j=1,k
      kk=j
    3 if (kk.lt.1) goto 4
      if (x(kk+m).ge.x(kk)) goto 4
      w=x(kk+m)
      x(kk+m)=x(kk)
      x(kk)=w
      kk=kk-m
      goto 3
    4 continue
      goto 2
      end


c--------------------------------------------------------------------
c
c     ipippa
c
c--------------------------------------------------------------------

       function ipippa(x,y,xc,yc,nc)
c
c point in polygon routine.
c
c returns 0 if point x,y not in the bound polygon defined by xc,yc
c
c fortran version of C routine by Ken McElvain
c

      implicit DOUBLE PRECISION (a-h,o-z)
c      include 'bounds.cmn'


      dimension xc(nc+1),yc(nc+1)

        iwind = 0
        xlastp = xc(nc)
        ylastp = yc(nc)
        ioldq = iquad(xlastp,ylastp,x,y)
        do i=1,nc
c for each point in the polygon
                xthisp=xc(i)
                ythisp=yc(i)
                inewq = iquad(xthisp,ythisp,x,y)
                if(ioldq.ne.inewq) then
                        if(mod(ioldq+1,4).eq.inewq) then
                          iwind=iwind+1
                        else if(mod(inewq+1,4).eq.ioldq) then
                          iwind = iwind - 1
                        else
                          a = (ylastp-ythisp)*(x-xlastp)
                          b = xlastp-xthisp
                          a = a + ylastp * b
                          b=b*y
                             if (a.gt.b) then
                               iwind=iwind+2
                             else
                               iwind=iwind-2
                             end if
                        end if
                end if
                xlastp=xthisp
                ylastp=ythisp
                ioldq=inewq
      end do
c
c quadrant winding is either -4,0,+4 so divide down and take abs.
c
      ipippa = abs(iwind/4)

      end

      function iquad(xp,yp,xo,yo)
c
c determine which quadrant xp,yp is in relative to xo,yo as origin
c
      implicit DOUBLE PRECISION (a-h,o-z)

        if(xp.lt.xo)then
                if(yp.lt.yo) then
                   iquad=2
                else
                   iquad=1
                end if
        else
                if(yp.lt.yo)then
                   iquad = 3
                else
                   iquad = 0
                end if
        end if

      return
      end

c===========================================================
c
c     Computation of a spatio-temporal covariance function
c
c
c     AUTHOR        :  E. Gabriel
c
c     DATE          :  24/01/2006
c
c     VERSION       :  1
c
c     last modification : 07/03/2007
c
c
c============================================================
c
      subroutine covst(gs,xx,yy,tt,nx,ny,nt,model,param,sigma2,scale,
     &                 aniso,ani)

c       INPUT:
c                xx,yy,tt     : spatial and temporal coordinates
c                model     : model
c                param(6)  : covariance parameters
c                aniso, ani : parameters for anisotropic covariances
c
c       OUTPUT:  gs        : covariance

      implicit none

      integer          nx,ny,nt
      integer          model(3)
      double precision gs(nx,ny,nt)
      double precision xx(nx),yy(ny),tt(nt),param(6)

      integer          i,j,k
      double precision covar,scale(2),sigma2,aniso,ani(2)

      do k=1,nt
        do j=1,ny
         do i=1,nx
            gs(i,j,k)=covar(xx(i),yy(j),tt(k),model,param,sigma2,scale,
     &                      aniso,ani)
         end do
        end do
      end do

      return
      end

c ---------------------------------------------------------------------------
c
c  program
c
      double precision function covar(x,y,t,model,param,sigma2,scale,
     &                                aniso,ani)

c     implicit none

      integer model(3)
      double precision x,y,t,dx,dt,param(6),sigma2,scale(2)
      double precision p1,p2,p3,p4,p5,p6
      double precision mod,mods,modt
      double precision aniso,ani(2), xani,yani
      double precision cauchy,stable,exponential,wave,matern
      double precision gneiting,cesare,theta(3)

      p1 = param(1)
      p2 = param(2)
      p3 = param(3)
      p4 = param(4)
      p5 = param(5)
      p6 = param(6)

      if (aniso.eq.1) then
        xani = x*cos(ani(1)) + y*sin(ani(1))
        yani = -x*sin(ani(1))/ani(2) + y*cos(ani(1))/ani(2)
c        dx=x*cos(ani(1))*x*cos(ani(1)) +
c     &     (x*sin(ani(1))/ani(2))*(x*sin(ani(1))/ani(2))
c     &     + 2*x*y*cos(ani(1))*sin(ani(1))
c     &     - 2*x*y*cos(ani(1))*sin(ani(1))/(ani(2)*ani(2))
c     &     + y*sin(ani(1))*y*sin(ani(1))
c     &     + (y*cos(ani(1))/ani(2))*(y*cos(ani(1))/ani(2))
         dx=dsqrt(xani*xani+yani*yani)/scale(1)
c        dx=dsqrt(dx)/scale(1)
      else
        dx=dsqrt(x*x+y*y)/scale(1)
      end if

      dt=dabs(t)/scale(2)

      mod=0d0
      mods=0d0
      modt=0d0

c
c     model= 0 (none)
c

      if(model(1).eq.0) then
         mods = 1d0
      endif

      if(model(2).eq.0) then
         modt = 1d0
      endif

c
c    model = 1 (Exponential)
c

      if(model(1).eq.1) then
         mods = exponential(dx)
      end if

      if(model(2).eq.1) then
         modt = exponential(dt)
      end if

c
c     model = 2 (Stable)
c

      if(model(1).eq.2) then
         mods = stable(dx,p1)
      end if

      if(model(2).eq.2) then
         modt = stable(dt,p2)
      end if

c
c     model = 3 (Cauchy)
c

      if(model(1).eq.3) then
         mods = cauchy(dx,p1)
      end if

      if(model(2).eq.3) then
         modt = cauchy(dt,p2)
      end if

c
c     model = 4 (Wave)
c

      if(model(1).eq.4) then
         mods = wave(dx)
      endif

      if(model(2).eq.4) then
         modt = wave(dt)
      endif

c
c     model = 7 (Matern)
c

      if(model(1).eq.7) then
         theta(1) = 1d0
         theta(2) = p3
         theta(3) = p1
         mods = matern(theta,dx)
      endif

      if(model(2).eq.7) then
         theta(1) = 1d0
         theta(2) = p4
         theta(3) = p2
         modt = matern(theta,dt)
      endif

c
c     product
c

      mod = mods * modt

c
c     model = 5 (Gneiting)
c

      if(model(3).eq.5) then
         mod = gneiting(dx,dt,param)
      endif

c
c     model = 6 (De Cesare)
c

      if(model(3).eq.6) then
         mod = cesare(dx,dt,p1,p2,p3)
      endif

c
c     result
c

      covar = sigma2*mod

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

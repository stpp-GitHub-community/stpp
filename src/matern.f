        subroutine rkmat(theta,x,y,n)
c
c
c     calculate Matern class covariances
c     If the smoothness is above 50 the Squared
c     Exponential limiting covariance is used.
c     
c     arguments:
c     theta     parameter vector (scale, range, smoothness)
c     x         vector of distances (input)
c     y         vector of Matern covariances (output)
c     n         length of x (and y)
c     
c     for a technical introduction see 
c     
c     Handcock and Stein (1993), 
c     'A Bayesian Analysis of Kriging',  
c     Technometrics, 35, 4, 403-410.
c
c     downloadable from:
c     http://www.stat.ncsu.edu/people/fuentes/st810m/lab/rkmat.f

        implicit double precision (a-h,o-z)
        parameter( zero = 0.0d0, one = 1.0d0 )
        parameter( half = 0.5d0, two = 2.0d0 )
        parameter( fifty = 50.0d0 )
c
        dimension theta(3),x(n),y(n)
        dimension bk(50)
c
        t1  = theta(1)
        t2l = one/theta(2)
        t3  = theta(3)
c
        t2  = two*dsqrt(t3)*t2l
c
        do 10 i = 1, n
          d = x(i) 
          if( d .le. zero ) then
             y(i) = t1
          else
            if(t3 .lt. fifty) then
              cnv = (two**(t3-one))*dgamma(t3)
              p1 = t2 * d
              it3 = dint(t3)
              call rkbesl(p1,t3-it3,it3+1,1,bk,ncalc)
              y(i) = t1*(p1**t3)*bk(it3+1)/cnv
            else
              p1 = t2l * d
              y(i) = t1*dexp(-p1*p1)        
            endif
          endif
10      continue
c
c       finish up
c
        return
        end


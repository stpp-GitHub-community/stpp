c=====================================================================
c
c     Simulation of Spatio-temporal stationary Gaussian random fields:
c     Computation of a the first row of the circulant embedding matrix
c
c
c     AUTHOR        :  E. Gabriel
c
c     DATE          :  09/02/2007
c
c     VERSION       :  2
c                      Modified: 17/03/2017 for anisotropic separable
c                      covariance functions
c
c
c=====================================================================

      subroutine circ(CEM,M,N,xlim,ylim,tlim,model,param,sigma2,scale,
     &                   aniso,ani)

      integer M(3), N(3)
      double precision xlim(2), ylim(2), tlim(2), sigma2, scale(2)
      double precision CEM(M(1)*M(2)*M(3)), MHALF(3)
      double precision COV3, gk
      integer model(3)
      double precision ani(2), aniso
      double precision param(6)

      do 2 I = 1,3
         MHALF(I)=M(I)/DBLE(2)
 2    continue

      do 5 K = 0, M(3) - 1
         do 6 J = 0, M(2) - 1
            do 7 I = 0, M(1) - 1
               NK = I+1+J*M(1)+K*M(1)*M(2)
               if ((I.le.MHALF(1)) .and. (J.le.MHALF(2))
     &              .and. (K.le.MHALF(3))) then
                CEM(NK) = COV3(gk(DBLE(I)/DBLE(N(1)),xlim,dble(N(1))),
     &                 gk(DBLE(J)/DBLE(N(2)),ylim,dble(N(2))),
     &                 gk(DBLE(K)/DBLE(N(3)),tlim,dble(N(3))),
     &                 model,param,sigma2,scale,aniso,ani)
               elseif ((J.le.MHALF(2)) .and. (K.le.MHALF(3))) then
                  CEM(NK) = CEM(M(1)-I+1+J*M(1)+K*M(1)*M(2))
               elseif (K.le.MHALF(3)) then
                  CEM(NK) = CEM(I+1+(M(2)-J)*M(1)+K*M(1)*M(2))
               else
                  CEM(NK) = CEM(I+1+J*M(1)+(M(3)-K)*M(1)*M(2))
               endif

 7          continue
 6       continue
 5    continue

      return
      end


c===========================================================

      double precision function gk(rk,lim,n)

      double precision rk,lim(2),n
      double precision a, b, delta, deltahalf

      a = lim(1)
      b = lim(2)
      delta = (b-a)/n
      deltahalf = delta/dble(2)

      gk = a + n*delta*rk + deltahalf

      return
      end


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
c
c============================================================
c

c
c  program
c
      double precision function COV3(x,y,t,model,param,sigma2,scale,
     &          aniso,ani)

c     implicit none

      integer model(3)
      double precision x,y,t,dx,dt,param(6),sigma2,scale(2)
      double precision xani, yani, ani(2), aniso
      double precision p1,p2,p3,p4,p5,p6
      double precision mod,mods,modt
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
cc      dx=(x*x+y*y)/scale(1)
cc      dt=(t**2)/scale(2)

      mod=0d0
      mods=0d0
      modt=0d0

c
c     model = 0 (none)
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

      COV3 = sigma2*mod

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c     Exponential
c
      double precision function exponential(x)

      double precision x

      exponential=dexp(-x)
      return
      end

c
c     Stable
c
      double precision function stable(x,p)

      double precision x,p

      stable=dexp(-x**p)
      return
      end

c
c     Cauchy
c

      double precision function cauchy(x,p)

      double precision x,p

      cauchy=(1+x**2)**(-p)
      return
      end

c
c     Wave
c
      double precision function wave(x)

      double precision x

      wave=0d0

      if (x.gt.0) then
         wave = dsin(x)/x
      endif
      if (x.eq.0) then
         wave = 1d0
      endif

      return
      end

c
c     Gneiting
c

      double precision function gneiting(x,t,param)

      double precision x,t,param(6)
      double precision psit,res
      double precision cauchy,stable
      double precision p1,p2,p3,p4,p5,p6

      res = 0d0
      p1 = param(1)
      p2 = param(2)
      p3 = param(3)
      p4 = param(4)
      p5 = param(5)
      p6 = param(6)

cc      if(p5.eq.1) psit=dsqrt((t**p3+1)**(p4))
cc      if(p5.eq.2) psit=dsqrt((1+(t**p3)/p4)/(1+(t**p3)))
cc      if(p5.eq.3) psit=dsqrt(-log(t**p3+1/p4)/log(p4))
cc      if(p2.eq.1) res = psit**(p6) * stable(x/psit,p1)
cc      if(p2.eq.2) res = psit**(p6) * cauchy(x/psit,p1)

      if(p5.eq.1) psit=(t**p3+1)**(p4)
      if(p5.eq.2) psit=(1+(t**p3)/p4)/(1+(t**p3))
      if(p5.eq.3) psit=log(t**p3+p4)/log(p4)
      if(p2.eq.1) res = psit**(-p6) * stable(x/psit,p1)
      if(p2.eq.2) res = psit**(-p6) * cauchy(x/psit,p1)

      gneiting = res

      return
      end

c
c     De Cesare
c

      double precision function cesare(x,t,p1,p2,p3)

      double precision x,t,p1,p2,p3

      cesare  = (1 + x**p1 + t**p2)**(-p3)

      return
      end

c
c     Matern
c

      double precision function matern(theta,x)
c
c
c     calculate Matern class covariances
c     If the smoothness is above 50 the Squared
c     Exponential limiting covariance is used.
c
c     arguments:
c     theta     parameter vector (scale, range, smoothness)
c     x         vector of distances (input)
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
        double precision theta(3),x
        dimension bk(50)
c
        t1  = theta(1)
        t2l = one/theta(2)
        t3  = theta(3)
c
        t2  = two*dsqrt(t3)*t2l
c        t2 = t21
c
        d = x
        if( d .le. zero ) then
           matern = t1
        else
           if(t3 .lt. fifty) then
              cnv = (two**(t3-one))*DGAMMAX(t3)
              p1 = t2 * d
              it3 = IDINT(t3)
              call rkbesl(p1,t3-it3,it3+1,1,bk,ncalc)
              matern = t1*(p1**t3)*bk(it3+1)/cnv
           else
              p1 = t2l * d
              matern = t1*dexp(-p1*p1)
           endif
        endif

c
c       finish up
c
        return
        end
      SUBROUTINE RKBESL(X,ALPHA,NB,IZE,BK,NCALC)
C-------------------------------------------------------------------
C
C  This FORTRAN 77 routine calculates modified Bessel functions
C  of the second kind, K SUB(N+ALPHA) (X), for non-negative
C  argument X, and non-negative order N+ALPHA, with or without
C  exponential scaling.
C
C  Explanation of variables in the calling sequence
C
C  Description of output values ..
C
C X     - Working precision non-negative real argument for which
C         K's or exponentially scaled K's (K*EXP(X))
C         are to be calculated.  If K's are to be calculated,
C         X must not be greater than XMAX (see below).
C ALPHA - Working precision fractional part of order for which
C         K's or exponentially scaled K's (K*EXP(X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
C         and 2 if exponentially scaled K's are to be calculated.
C BK    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BK
C         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
C         or the corresponding exponentially scaled functions.
C         If (0 .LT. NCALC .LT. NB), BK(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BK, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C-------------------------------------------------------------------
      INTEGER I,IEND,ITEMP,IZE,J,K,M,MPLUS1,NB,NCALC
      DOUBLE PRECISION
     1    A,ALPHA,BLPHA,BK,BK1,BK2,C,D,DM,D1,D2,D3,ENU,EPS,ESTF,ESTM,
     2    EX,FOUR,F0,F1,F2,HALF,ONE,P,P0,Q,Q0,R,RATIO,S,SQXMIN,T,TINYX,
     3    TWO,TWONU,TWOX,T1,T2,WMINF,X,XINF,XMAX,XMIN,X2BY4,ZERO
      DIMENSION BK(1),P(8),Q(7),R(5),S(4),T(6),ESTM(6),ESTF(7)
C---------------------------------------------------------------------
C  Mathematical constants
C    A = LOG(2.D0) - Euler's constant
C    D = SQRT(2.D0/PI)
C---------------------------------------------------------------------
      DATA HALF,ONE,TWO,ZERO/0.5D0,1.0D0,2.0D0,0.0D0/
      DATA FOUR,TINYX/4.0D0,1.0D-10/
      DATA A/ 0.11593151565841244881D0/,D/0.797884560802865364D0/
C---------------------------------------------------------------------
C  Machine dependent parameters
C---------------------------------------------------------------------
      DATA EPS/2.22D-16/,SQXMIN/1.49D-154/,XINF/1.79D+308/
      DATA XMIN/2.23D-308/,XMAX/705.342D0/
C---------------------------------------------------------------------
C  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
C                                         + Euler's constant
C         Coefficients converted from hex to decimal and modified
C         by W. J. Cody, 2/26/82
C  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
C  T    - Approximation for SINH(Y)/Y
C---------------------------------------------------------------------
      DATA P/ 0.805629875690432845D00,    0.204045500205365151D02,
     1        0.157705605106676174D03,    0.536671116469207504D03,
     2        0.900382759291288778D03,    0.730923886650660393D03,
     3        0.229299301509425145D03,    0.822467033424113231D00/
      DATA Q/ 0.294601986247850434D02,    0.277577868510221208D03,
     1        0.120670325591027438D04,    0.276291444159791519D04,
     2        0.344374050506564618D04,    0.221063190113378647D04,
     3        0.572267338359892221D03/
      DATA R/-0.48672575865218401848D+0,  0.13079485869097804016D+2,
     1       -0.10196490580880537526D+3,  0.34765409106507813131D+3,
     2        0.34958981245219347820D-3/
      DATA S/-0.25579105509976461286D+2,  0.21257260432226544008D+3,
     1       -0.61069018684944109624D+3,  0.42269668805777760407D+3/
      DATA T/ 0.16125990452916363814D-9, 0.25051878502858255354D-7,
     1        0.27557319615147964774D-5, 0.19841269840928373686D-3,
     2        0.83333333333334751799D-2, 0.16666666666666666446D+0/
      DATA ESTM/5.20583D1, 5.7607D0, 2.7782D0, 1.44303D1, 1.853004D2,
     1          9.3715D0/
      DATA ESTF/4.18341D1, 7.1075D0, 6.4306D0, 4.25110D1, 1.35633D0,
     1          8.45096D1, 2.0D1/
C---------------------------------------------------------------------
      ITEMP=0d0

      EX = X
      ENU = ALPHA
      NCALC = MIN(NB,0)-2
      IF ((NB .GT. 0) .AND. ((ENU .GE. ZERO) .AND. (ENU .LT. ONE))
     1     .AND. ((IZE .GE. 1) .AND. (IZE .LE. 2)) .AND.
     2     ((IZE .NE. 1) .OR. (EX .LE. XMAX)) .AND.
     3     (EX .GT. ZERO))  THEN
            K = 0
            IF (ENU .LT. SQXMIN) ENU = ZERO
            IF (ENU .GT. HALF) THEN
                  K = 1
                  ENU = ENU - ONE
            END IF
            TWONU = ENU+ENU
            IEND = NB+K-1
            C = ENU*ENU
            D3 = -C
            IF (EX .LE. ONE) THEN
C---------------------------------------------------------------------
C  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
C                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
C---------------------------------------------------------------------
                  D1 = ZERO
                  D2 = P(1)
                  T1 = ONE
                  T2 = Q(1)
                  DO 10 I = 2,7,2
                     D1 = C*D1+P(I)
                     D2 = C*D2+P(I+1)
                     T1 = C*T1+Q(I)
                     T2 = C*T2+Q(I+1)
   10             CONTINUE
                  D1 = ENU*D1
                  T1 = ENU*T1
                  F1 = LOG(EX)
                  F0 = A+ENU*(P(8)-ENU*(D1+D2)/(T1+T2))-F1
                  Q0 = EXP(-ENU*(A-ENU*(P(8)+ENU*(D1-D2)/(T1-T2))-F1))
                  F1 = ENU*F0
                  P0 = EXP(F1)
C---------------------------------------------------------------------
C  Calculation of F0 =
C---------------------------------------------------------------------
                  D1 = R(5)
                  T1 = ONE
                  DO 20 I = 1,4
                     D1 = C*D1+R(I)
                     T1 = C*T1+S(I)
   20             CONTINUE
                  IF (ABS(F1) .LE. HALF) THEN
                        F1 = F1*F1
                        D2 = ZERO
                        DO 30 I = 1,6
                           D2 = F1*D2+T(I)
   30                   CONTINUE
                        D2 = F0+F0*F1*D2
                     ELSE
                        D2 = SINH(F1)/ENU
                  END IF
                  F0 = D2-ENU*D1/(T1*P0)
                  IF (EX .LE. TINYX) THEN
C--------------------------------------------------------------------
C  X.LE.1.0E-10
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                        BK(1) = F0+EX*F0
                        IF (IZE .EQ. 1) BK(1) = BK(1)-EX*BK(1)
                        RATIO = P0/F0
                        C = EX*XINF
                        IF (K .NE. 0) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
C  ALPHA .GE. 1/2
C--------------------------------------------------------------------
                              NCALC = -1
                              IF (BK(1) .GE. C/RATIO) GO TO 500
                              BK(1) = RATIO*BK(1)/EX
                              TWONU = TWONU+TWO
                              RATIO = TWONU
                        END IF
                        NCALC = 1
                        IF (NB .EQ. 1) GO TO 500
C--------------------------------------------------------------------
C  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
C--------------------------------------------------------------------
                        NCALC = -1
                        DO 80 I = 2,NB
                           IF (RATIO .GE. C) GO TO 500
                           BK(I) = RATIO/EX
                           TWONU = TWONU+TWO
                           RATIO = TWONU
   80                   CONTINUE
                        NCALC = 1
                        GO TO 420
                     ELSE
C--------------------------------------------------------------------
C  1.0E-10 .LT. X .LE. 1.0
C--------------------------------------------------------------------
                        C = ONE
                        X2BY4 = EX*EX/FOUR
                        P0 = HALF*P0
                        Q0 = HALF*Q0
                        D1 = -ONE
                        D2 = ZERO
                        BK1 = ZERO
                        BK2 = ZERO
                        F1 = F0
                        F2 = P0
  100                   D1 = D1+TWO
                        D2 = D2+ONE
                        D3 = D1+D3
                        C = X2BY4*C/D2
                        F0 = (D2*F0+P0+Q0)/D3
                        P0 = P0/(D2-ENU)
                        Q0 = Q0/(D2+ENU)
                        T1 = C*F0
                        T2 = C*(P0-D2*F0)
                        BK1 = BK1+T1
                        BK2 = BK2+T2
                        IF ((ABS(T1/(F1+BK1)) .GT. EPS) .OR.
     1                     (ABS(T2/(F2+BK2)) .GT. EPS))  GO TO 100
                        BK1 = F1+BK1
                        BK2 = TWO*(F2+BK2)/EX
                        IF (IZE .EQ. 2) THEN
                              D1 = EXP(EX)
                              BK1 = BK1*D1
                              BK2 = BK2*D1
                        END IF
                        WMINF = ESTF(1)*EX+ESTF(2)
                  END IF
               ELSE IF (EPS*EX .GT. ONE) THEN
C--------------------------------------------------------------------
C  X .GT. ONE/EPS
C--------------------------------------------------------------------
                  NCALC = NB
                  BK1 = ONE / (D*SQRT(EX))
                  DO 110 I = 1, NB
                     BK(I) = BK1
  110             CONTINUE
                  GO TO 500
               ELSE
C--------------------------------------------------------------------
C  X .GT. 1.0
C--------------------------------------------------------------------
                  TWOX = EX+EX
                  BLPHA = ZERO
                  RATIO = ZERO
                  IF (EX .LE. FOUR) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(1)/EX+ESTM(2))
                        M = INT(D2)
                        D1 = D2+D2
                        D2 = D2-HALF
                        D2 = D2*D2
                        DO 120 I = 2,M
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
  120                   CONTINUE
C--------------------------------------------------------------------
C  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
C    recurrence and K(ALPHA,X) from the wronskian
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(3)*EX+ESTM(4))
                        M = INT(D2)
                        C = ABS(ENU)
                        D3 = C+C
                        D1 = D3-ONE
                        F1 = XMIN
                        F0 = (TWO*(C+D2)/EX+HALF*EX/(C+D2+ONE))*XMIN
                        DO 130 I = 3,M
                           D2 = D2-ONE
                           F2 = (D3+D2+D2)*F0
                           BLPHA = (ONE+D1/D2)*(F2+BLPHA)
                           F2 = F2/EX+F1
                           F1 = F0
                           F0 = F2
  130                   CONTINUE
                        F1 = (D3+TWO)*F0/EX+F1
                        D1 = ZERO
                        T1 = ONE
                        DO 140 I = 1,7
                           D1 = C*D1+P(I)
                           T1 = C*T1+Q(I)
  140                   CONTINUE
                        P0 = EXP(C*(A+C*(P(8)-C*D1/T1)-LOG(EX)))/EX
                        F2 = (C+HALF-RATIO)*F1/EX
                        BK1 = P0+(D3*F0-F2+F0+BLPHA)/(F2+F1+F0)*P0
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(3)*EX+ESTF(4)
                     ELSE
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
C  recurrence, for  X .GT. 4.0
C--------------------------------------------------------------------
                        DM = AINT(ESTM(5)/EX+ESTM(6))
                        M = INT(DM)
                        D2 = DM-HALF
                        D2 = D2*D2
                        D1 = DM+DM
                        DO 160 I = 2,M
                           DM = DM-ONE
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
                           BLPHA = (RATIO+RATIO*BLPHA)/DM
  160                   CONTINUE
                        BK1 = ONE/((D+D*BLPHA)*SQRT(EX))
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(5)*(EX-ABS(EX-ESTF(7)))+ESTF(6)
                  END IF
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
C    K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                  BK2 = BK1+BK1*(ENU+HALF-RATIO)/EX
            END IF
C--------------------------------------------------------------------
C  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
C  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
C--------------------------------------------------------------------
            NCALC = NB
            BK(1) = BK1
            IF (IEND .EQ. 0) GO TO 500
            J = 2-K
            IF (J .GT. 0) BK(J) = BK2
            IF (IEND .EQ. 1) GO TO 500
            M = MIN(INT(WMINF-ENU),IEND)
            DO 190 I = 2,M
               T1 = BK1
               BK1 = BK2
               TWONU = TWONU+TWO
               IF (EX .LT. ONE) THEN
                     IF (BK1 .GE. (XINF/TWONU)*EX) GO TO 195
                     GO TO 187
                  ELSE
                     IF (BK1/EX .GE. XINF/TWONU) GO TO 195
               END IF
  187          CONTINUE
               BK2 = TWONU/EX*BK1+T1
               ITEMP = I
               J = J+1
               IF (J .GT. 0) BK(J) = BK2
  190       CONTINUE
  195       M = ITEMP
            IF (M .EQ. IEND) GO TO 500
            RATIO = BK2/BK1
            MPLUS1 = M+1
            NCALC = -1
            DO 410 I = MPLUS1,IEND
               TWONU = TWONU+TWO
               RATIO = TWONU/EX+ONE/RATIO
               J = J+1
               IF (J .GT. 1) THEN
                     BK(J) = RATIO
                  ELSE
                     IF (BK2 .GE. XINF/RATIO) GO TO 500
                     BK2 = RATIO*BK2
               END IF
  410       CONTINUE
            NCALC = MAX(MPLUS1-K,1)
            IF (NCALC .EQ. 1) BK(1) = BK2
            IF (NB .EQ. 1) GO TO 500
  420       J = NCALC+1
            DO 430 I = J,NB
               IF (BK(NCALC) .GE. XINF/BK(I)) GO TO 500
               BK(I) = BK(NCALC)*BK(I)
               NCALC = I
  430       CONTINUE
      END IF
  500 RETURN
C---------- Last line of RKBESL ----------
      END
      DOUBLE PRECISION FUNCTION DGAMMAX(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      DOUBLE PRECISION
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
  900 DGAMMAX = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END

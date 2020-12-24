C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     spatial r-mark function.
C

       subroutine kmrcore(x,y,txy,n,s,ns,ks,hs,kmr)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision wij,vij,hs,kerns,krm,krn,kmr,x,y,txy
       double precision hij,mij,xi,yi,ti,two
       dimension x(n),y(n),txy(n),s(ns),krm(ns),krn(ns),kmr(ns),ks(3)

       krm=0d0
       krn=0d0

          two=2d0

       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
            mij=ti
              if (ks(1).eq.1) then
               kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks(2).eq.1) then
                 kerns=ekernel((s(iu)-hij)/hs,hs)
                  else if (ks(3).eq.1) then
                   kerns=qkernel((s(iu)-hij)/hs,hs)
              end if
             if (kerns.ne.0d0) then
                    wij=mij*kerns
                    vij=kerns
                    krm(iu)=krm(iu)+wij
                    krn(iu)=krn(iu)+vij
             end if
           end if
          end do
          end do
            kmr(iu)=krm(iu)/krn(iu)
          end do

        return

        end

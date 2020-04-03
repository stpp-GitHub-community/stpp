C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     spatial mark variogram function.
C

       subroutine gspcore(x,y,txy,n,s,ns,ks,hs,gsps)

      implicit double precision (a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision wij,vij,hs,kerns,gspm,gspn,gsps,x,y,txy
       double precision hij,mij,xi,yi,ti,two
       dimension x(n),y(n),txy(n),s(ns),gspm(ns),gspn(ns),gsps(ns),ks(3)

       gspm=0d0
       gspn=0d0

          two=2d0

       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
            mij=((abs(ti-txy(j)))**two)/two
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
                    gspm(iu)=gspm(iu)+wij
                    gspn(iu)=gspn(iu)+vij
             end if
           end if
          end do
          end do
            gsps(iu)=gspm(iu)/gspn(iu)
          end do

        return

        end
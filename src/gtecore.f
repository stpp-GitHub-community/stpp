C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal mark variogram function.
C

       subroutine gtecore(x,y,txy,n,t,nt,kt,ht,gtet)


       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,gtem,gten,gtet,x,y,txy
       double precision tij,mij,ti,xi,yi
       dimension x(n),y(n),txy(n),t(nt),gtem(nt),gten(nt),gtet(nt),kt(3)

       gtem=0d0
       gten=0d0

          two=2d0

       do iv=1,nt
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=abs(ti-txy(j))
            mij=((sqrt(((xi-x(j))**two)+((yi-y(j))**two)))**two)/two
              if (kt(1).eq.1) then
               kernt=boxkernel((t(iv)-tij)/ht,ht)
                else if (kt(2).eq.1) then
                 kernt=ekernel((t(iv)-tij)/ht,ht)
                  else if (kt(3).eq.1) then
                   kernt=qkernel((t(iv)-tij)/ht,ht)
              end if
             if (kernt.ne.0d0) then
                    wij=mij*kernt
                    vij=kernt
                    gtem(iv)=gtem(iv)+wij
                    gten(iv)=gten(iv)+vij
             end if
           end if
          end do
          end do
            gtet(iv)=gtem(iv)/gten(iv)
          end do

        return

        end

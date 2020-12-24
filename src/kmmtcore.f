C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides a non-parametric kernel based estimator of the
C     temporal mark correlation function.
C

       subroutine kmmtcore(snorm,txy,n,t,nt,kt,ht,kmmt)

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,ktm,ktn,kmmt,snorm,txy
       double precision tij,mij,snormi,ti
       dimension snorm(n),txy(n),t(nt),ktm(nt),ktn(nt),kmmt(nt),kt(3)

       ktm=0d0
       ktn=0d0

       do iv=1,nt
        do i=1,n
         snormi=snorm(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            tij=abs(ti-txy(j))
            mij=snormi*snorm(j)
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
                    ktm(iv)=ktm(iv)+wij
                    ktn(iv)=ktn(iv)+vij
             end if
           end if
          end do
          end do
            kmmt(iv)=ktm(iv)/ktn(iv)
          end do

        return

        end

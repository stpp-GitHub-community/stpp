C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides an edge-corrected non-parametric kernel based
C     estimator of the standardized temporal t-mark function. 
C

      subroutine kmtcoreinh(snorm,txy,n,t,nt,tlambda,kt,ht,wrt,wtt,
     +     wbit,wbimodt,wst,edg,kmt)
     
      implicit real*8(a-h,o-z)

      integer i,j,iv,n,nt,kt,edg
      double precision inhwij,inhvij,ht,kernt,kmtminh,kmtninh,kmt
      double precision tij,mij,snormi,ti,wrt,wtt,wbit,txy,snorm
      double precision wbimodt,wst,tlambda
      dimension snorm(n),txy(n),t(nt),kmtminh(nt),kmtninh(nt)
      dimension wrt(n,n),wtt(n,n),wbit(n,nt),wbimodt(n,nt),wst(nt)
      dimension kt(3),edg(6),tlambda(n),kmt(nt)
       
       kmtminh=0d0
       kmtninh=0d0

      do iv=1,nt
      do i=1,n
         snormi=snorm(i)
         ti=txy(i)      
      do j=1,n
      if (j.ne.i) then
      tij=abs(ti-txy(j))
      mij=snormi
      if (kt(1).eq.1) then
      kernt=boxkernel((t(iv)-tij)/ht,ht)
      else if (kt(2).eq.1) then
      kernt=ekernel((t(iv)-tij)/ht,ht)
      else if (kt(3).eq.1) then
      kernt=qkernel((t(iv)-tij)/ht,ht)
      end if
      if (kernt.ne.0d0) then
C    none
      if (edg(1).eq.1) then
       inhwij=(mij*kernt)/(tlambda(i)*tlambda(j))
       inhvij=kernt/(tlambda(i)*tlambda(j))
       kmtminh(iv)=kmtminh(iv)+inhwij
       kmtninh(iv)=kmtninh(iv)+inhvij  
      end if                  
C    isotropic
      if (edg(2).eq.1) then                  
      inhwij=(mij*kernt*wrt(i,j))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wrt(i,j))/(tlambda(i)*tlambda(j))
      kmtminh(iv)=kmtminh(iv)+inhwij
      kmtninh(iv)=kmtninh(iv)+inhvij
      end if
C    border
      if (edg(3).eq.1) then                  
      inhwij=(mij*kernt*wbit(i,iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wbit(i,iv))/(tlambda(i)*tlambda(j))
      kmtminh(iv)=kmtminh(iv)+inhwij
      kmtninh(iv)=kmtninh(iv)+inhvij
      end if
C    modified.border
      if (edg(4).eq.1) then
      inhwij=(mij*kernt*wbimodt(i,iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wbimodt(i,iv))/(tlambda(i)*tlambda(j))
      kmtminh(iv)=kmtminh(iv)+inhwij
      kmtninh(iv)=kmtninh(iv)+inhvij 
      end if                  
C    translate
      if (edg(5).eq.1) then
      inhwij=(mij*kernt*wtt(i,j))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wtt(i,j))/(tlambda(i)*tlambda(j))
      kmtminh(iv)=kmtminh(iv)+inhwij
      kmtninh(iv)=kmtninh(iv)+inhvij
      end if
C    setcovf         
      if (edg(6).eq.1) then
      inhwij=(mij*kernt*wst(iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wst(iv))/(tlambda(i)*tlambda(j))
      kmtminh(iv)=kmtminh(iv)+inhwij
      kmtninh(iv)=kmtninh(iv)+inhvij
      end if
      end if
      end if
       end do
       end do
       kmt(iv)=kmtminh(iv)/kmtninh(iv)
       end do
      
        return
        
        end  

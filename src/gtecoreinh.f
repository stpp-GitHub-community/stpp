C     Francisco J. Rodriguez-Cortes, November 2016
C
C     This code provides an edge-corrected non-parametric kernel based
C     estimator of the standardized temporal mark variogram function. 
C

      subroutine gtecoreinh(x,y,txy,n,t,nt,tlambda,kt,ht,wrt,wtt,
     +     wbit,wbimodt,wst,edg,gtet)
     
      implicit real*8(a-h,o-z)

      integer i,j,iv,n,nt,kt,edg
      double precision inhwij,inhvij,ht,kernt,gteminh,gteninh,gtet
      double precision tij,mij,xi,yi,ti,wrt,wtt,wbit,x,y,txy
      double precision wbimodt,wst,tlambda
      dimension x(n),y(n),txy(n),t(nt),gteminh(nt),gteninh(nt)
      dimension wrt(n,n),wtt(n,n),wbit(n,nt),wbimodt(n,nt),wst(nt)
      dimension kt(3),edg(6),tlambda(n),gtet(nt)
       
       gteminh=0d0
       gteninh=0d0

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
C    none
      if (edg(1).eq.1) then
       inhwij=(mij*kernt)/(tlambda(i)*tlambda(j))
       inhvij=kernt/(tlambda(i)*tlambda(j))
       gteminh(iv)=gteminh(iv)+inhwij
       gteninh(iv)=gteninh(iv)+inhvij  
      end if                  
C    isotropic
      if (edg(2).eq.1) then                  
      inhwij=(mij*kernt*wrt(i,j))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wrt(i,j))/(tlambda(i)*tlambda(j))
      gteminh(iv)=gteminh(iv)+inhwij
      gteninh(iv)=gteninh(iv)+inhvij
      end if
C    border
      if (edg(3).eq.1) then                  
      inhwij=(mij*kernt*wbit(i,iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wbit(i,iv))/(tlambda(i)*tlambda(j))
      gteminh(iv)=gteminh(iv)+inhwij
      gteninh(iv)=gteninh(iv)+inhvij
      end if
C    modified.border
      if (edg(4).eq.1) then
      inhwij=(mij*kernt*wbimodt(i,iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wbimodt(i,iv))/(tlambda(i)*tlambda(j))
      gteminh(iv)=gteminh(iv)+inhwij
      gteninh(iv)=gteninh(iv)+inhvij
      end if                  
C    translate
      if (edg(5).eq.1) then
      inhwij=(mij*kernt*wtt(i,j))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wtt(i,j))/(tlambda(i)*tlambda(j))
      gteminh(iv)=gteminh(iv)+inhwij
      gteninh(iv)=gteninh(iv)+inhvij
      end if
C    setcovf         
      if (edg(6).eq.1) then
      inhwij=(mij*kernt*wst(iv))/(tlambda(i)*tlambda(j))
      inhvij=(kernt*wst(iv))/(tlambda(i)*tlambda(j))
      gteminh(iv)=gteminh(iv)+inhwij
      gteninh(iv)=gteninh(iv)+inhvij
      end if
      end if
      end if
       end do
       end do
       gtet(iv)=gteminh(iv)/gteninh(iv)
       end do
      
        return
        
        end  

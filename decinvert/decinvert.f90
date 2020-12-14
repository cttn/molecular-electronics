module decimod
contains
!%%%%%%%%%%% ELEMENTOS DE LA INVERSA %%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dec_inv_elem(dN,input,el1,el2,Inv12,Inv21) !Inv21 es opcional
 implicit none
 integer,intent(in) :: dN,el1,el2 
 real(kind(0d0)),intent(in) :: input(dN,dN)
 real(kind(0d0)),intent(out) :: Inv12
 real(kind(0d0)),optional,intent(out) :: Inv21
 integer :: i,j,d
 real(kind(0d0)) :: A(dN,dN)
 logical :: decimado(dN)
 if ((el1 .lt. 1 .or. el1 .gt. dN) .or. (el2 .lt. 1 .or. el2 .gt. dN)) &
	stop 'Error: el elemento 1 o el 2 es incorrecto'
 decimado=.false.					
 A=input						
 do d=1,dN						
   if (d .eq. el1 .or. d .eq. el2 ) cycle		
   do i=1,dN
     do j=1,dN
       decimado(d)=.true.				
       if((.not. decimado(i)).and.(.not. decimado(j))) then
	 a(i,j)=a(i,j)-a(i,d)*a(d,j)/a(d,d)		
       endif
     enddo
   enddo
 enddo
 if (el1 .ne. el2) then 				
   Inv12=-a(el1,el2)/(a(el1,el1)*a(el2,el2)-a(el1,el2)*a(el2,el1)) 
   if(present(Inv21)) Inv21=-a(el2,el1)/(a(el2,el2)*a(el1,el1)-a(el2,el1)*a(el1,el2))
 else
   Inv12=1.d0/a(el1,el1)
 endif
end subroutine dec_inv_elem

!%%%%%%%%%%% INVERSA TOTAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dec_inv_total(dN,input,output)
 implicit none
 integer,intent(in) :: dN !Dimension
 real(kind(0d0)),intent(in) :: input(dN,dN)
 real(kind(0d0)),intent(out) :: output(dN,dN)
 integer :: l1,l2,i,j,d
 real(kind(0d0)) :: inv(dN,dN),A(dN,dN)
 logical :: decimado(dN)
 do l1=1,dN						
   do l2=l1,dN						
     decimado=.false.					
     A=input						
     do d=1,dN						
       if (d .eq. l1 .or. d .eq. l2 ) cycle		
       do i=1,dN
	 do j=1,dN
	   decimado(d)=.true.				
           if((.not. decimado(i)).and.(.not. decimado(j))) then
	     a(i,j)=a(i,j)-a(i,d)*a(d,j)/a(d,d)		
	   endif
         enddo
       enddo
     enddo
     if (l1 .ne. l2) then 				
       output(l1,l2)=-a(l1,l2)/(a(l1,l1)*a(l2,l2)-a(l1,l2)*a(l2,l1)) 
       output(l2,l1)=-a(l2,l1)/(a(l2,l2)*a(l1,l1)-a(l2,l1)*a(l1,l2))
     else
       output(l1,l1)=1.d0/a(l1,l1)
     endif
   enddo
 enddo
end subroutine dec_inv_total
end module decimod

program decinvert
 use decimod
 implicit none
 integer,parameter :: N=10
 integer :: i,j,d,l1,l2,idum=1131
 real(kind(0d0)) :: ran2,checkij,checkji
 real(kind(0d0)) :: B(N,N),Inv(N,N),check(N,N)

 do i=1,N
   do j=1,N
   B(i,j)=ran2(idum)*100.
   enddo
 enddo

 call dec_inv_total(N,B,inv)
 check=matmul(B,inv)

 call dec_inv_elem(N,B,1,3,checkij)
   print *, checkij,inv(1,3)

 do i=1,N
   write(88,111) (check(i,j),j=1,N)
 enddo
 
 111 format(200E12.4)
end program decinvert




function ran2(idum)
 integer,parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1,IA1=40014,	&
	IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB
 integer :: idum,idum2,j,k,iv(NTAB),iy
 real(kind(0d0)),parameter :: AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS
 real(kind(0d0)) :: ran2
 save iv,iy,idum2
 data idum2/123456789/, iv/NTAB*0/, iy/0/ 
  if (idum.le.0) then
    idum=max(-idum,1)
    idum2=idum
      do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
      enddo
    iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
end function ran2

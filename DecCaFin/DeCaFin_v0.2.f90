module global_DeCaFin
 implicit none
  integer,parameter :: puntos_ener=1,max_sitios=5000
  integer,save :: num_sitios,semilla
  integer :: s
  complex(8),save :: ener(max_sitios,puntos_ener),hop(max_sitios-1,puntos_ener),	&
	delta_lateral(max_sitios,puntos_ener),epsi(puntos_ener)
  real(8), parameter :: Pi=dacos(-1.0d0),ener_lead=0.0d0,hop_lead=3.0d0
  real(8),save :: Emin, Emax, eta,dos(puntos_ener),Trans(puntos_ener)
  data dos /puntos_ener*0.0d0/

  contains
  function sigma(e)
     complex(8) :: sigma
     integer :: e
     if    ((real(epsi(e))-ener_lead) .gt. 2.0d0*hop_lead) then
	sigma=dcmplx((real(epsi(e))-ener_lead)/2.0d0                  &
        -sqrt(((real(epsi(e))-ener_lead)/2.0d0)**2-hop_lead**2),0.0d0)   
     elseif(dabs(real(epsi(e))-ener_lead) .le. dabs(2.0d0*hop_lead)) then
	sigma=dcmplx((real(epsi(e))-ener_lead)/2.0d0,-sqrt(hop_lead**2 &
        -((real(epsi(e))-ener_lead)/2.0d0)**2))
     elseif((real(epsi(e))-ener_lead) .lt. -2.0d0*dabs(hop_lead)) then
	sigma=dcmplx((real(epsi(e))-ener_lead)/2.0d0 +                &
        sqrt(((real(epsi(e))-ener_lead)/2.0d0)**2-hop_lead**2),0.0d0)
     endif
  end function sigma
end module global_DeCaFin



program DeCaFin
 use global_DeCaFin
!*******************************************************************************!
!  Caldula transmitancia "coherente" de una cadena unidimensional desordenada   !
! en la que cada sitio está conectado a uuna cadena finita arbitraria		!
!*******************************************************************************!
 implicit none


!-- Parametros Generales 
  eta=0.0000000001d0
  Emin=0.00d0
  Emax=0.00d0
  semilla=1776776

!-- Parametros Generales -----------------------
!  num_sitios=1000 !Numero de sitios en el sistema (cadena horizontal)


  do s=1,puntos_ener
   epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s)/real(puntos_ener),eta)
  enddo

  open(unit=55,file='T_vs_L.dat')

   do num_sitios=10,5000,10

    call asignaciones
 
    call calculos

  write(55,*)num_sitios,trans(1)

  enddo


end program DeCaFin



!_______________________________________________________________________________!
subroutine asignaciones
 use global_DeCaFin
 implicit none
  integer,parameter :: max_lateral=5000
  integer :: largo_lateral(max_sitios),i,j
  complex(8) :: delta_fin(0:max_lateral,puntos_ener),ener_sitio(max_sitios),	&
	hop_sitio(max_sitios),ener_lateral(max_lateral),hop_lateral(0:max_lateral)
  real(8) :: W,ran3
 

  W=0.2d0 !Parametro de desorden para energias del sistema

!*** CADENA LATERAL ****************************

  do i=1,1!num_sitios

!-- Largo de cada cadena lateral 
   largo_lateral(i)=5000

!-- Tunneling a la cadena lateral
   hop_lateral(0)=0.0d0

!-- Asignacion de parametros laterales
  do j=1,largo_lateral(i)
   ener_lateral(j)=0.0d0
   if (j .lt. largo_lateral(i)) then
    hop_lateral(j)=1.0d0
   endif
  enddo


!-- Delta lateral
   do s=1,puntos_ener

    delta_fin(largo_lateral(i),s)=dcmplx(0.0d0,0.0d0)
    do j=largo_lateral(i)-1,0,-1
     delta_fin(j,s)=( hop_lateral(j)*hop_lateral(j)/(epsi(s)-ener_lateral(j+1)-delta_fin(j+1,s)) )
    enddo!j
    delta_lateral(i,s)=delta_fin(0,s)

   enddo!s

  enddo !i



!*** CADENA SISTEMA (HORIZONTAL) *****************
 do i=1,num_sitios
  ener_sitio(i)=(ran3(semilla)-0.5)*W
  if (i .lt. num_sitios) then
   hop_sitio(i)=1.0d0
  endif
 enddo



!*** ASIGNACION FINAL ****************************
 do i=1,num_sitios
  do s=1,puntos_ener
   ener(i,s)=ener_sitio(i)+delta_lateral(1,s)!!!!!!!!
   if (i .lt. max_sitios) then
    hop(i,s)=hop_sitio(i)
   endif
  enddo!s
 enddo !i

end subroutine asignaciones



!_______________________________________________________________________________!
subroutine calculos
 use global_DeCaFin
 implicit none

  integer :: n,h,i,j !sólo para la info de avanze de la subrutina
  complex(8) :: Green(max_sitios,max_sitios),ener_renor(max_sitios,puntos_ener)
  type self
  complex(8) :: mas(0:max_sitios+1,puntos_ener),men(0:max_sitios+1,puntos_ener)
  end type self
  type (self) delta


!-- Correccion a energia de sitio por renormalizacion
  do s=1,puntos_ener
      delta%men(0,s)=sigma(s)  !Sigma es la Self-energy debida al lead.
      delta%men(1,s)=sigma(s)
      delta%mas(num_sitios,s)=sigma(s)
      delta%mas(num_sitios+1,s)=sigma(s)

  do j=2,num_sitios+1
      delta%men(j,s)= (hop(j-1,s)*hop(j-1,s))                  &
      /(epsi(s)-ener(j-1,s)-delta%men(j-1,s))

      delta%mas(num_sitios+1-j,s)=(hop(num_sitios+1-j,s)             &
      *hop(num_sitios +1-j,s)/(epsi(s)-ener(num_sitios+2-j,s)  &
      -delta%mas(num_sitios+2-j,s)))

  enddo !j


!-- Energias renormalizadas
  do j=1,num_sitios
      ener_renor(j,s)=ener(j,s)+delta%mas(j,s)+delta%men(j,s)
      dos(s)=dos(s)+(-1.0d0/pi)*aimag(1.0d0/(epsi(s)-ener_renor(j,s)))
  enddo !j
  enddo !s

!-- Funciones de Green sitio a sitio
  do s=1,puntos_ener

      do i=1,num_sitios
       Green(i,i)=1.0d0/(epsi(s)-ener_renor(i,s))
      enddo

      do i=2,num_sitios+1
       do j=1,i-1
        Green(j,i)=Green(j,i-1)*delta%mas(i-1,s)/hop(i-1,s)
       enddo
      enddo

   trans(s)=4.0d0*(cdabs(green(1,num_sitios))**2)*aimag(sigma(s))**2


! Info de avanze de esta subrutina
   do n=20,100,20
   if (n .eq. int(real(s)*100./real(puntos_ener))) then
        if (h .eq. n) then
		goto 118
        else
		print *, n,'% comletado'
        	h=n
        endif
   endif
   enddo
 118 continue

  enddo !s

!---
!  open(unit=34,file='DeCaFin.dat')
!   do s=1,puntos_ener
!    write(34,*)real(epsi(s)),trans(s)!,dos(s)
!   enddo
!  close(43)

end subroutine calculos





!*********************************************************
!*      Generador de numeros aleatorios                  *
!*********************************************************
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
!     PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
 11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12        continue
 13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END FUNCTION
!  (C) Copr. 1986-92 Numerical Recipes Software.

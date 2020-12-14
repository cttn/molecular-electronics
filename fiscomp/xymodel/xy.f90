!##########################################################!
!###		FÍSICA COMPUTACIONAL: FINAL		###!
!### 	  Simulación Monte Carlo para el modelo XY	###!
!##########################################################!
module Global_XY
use mtmod
implicit none
 integer,save :: seed=15791,L,paso,paso_acep,N_MCS,interv,	&
	N_muestreo,instan_cont=0
 real(kind(0d0)),parameter :: pi=dacos(-1.d0)
 real(kind(0d0)),save :: E,T,sum_E,sum_E2,sum_MH,sum_MH2,sum_EMH
 real(kind(0d0)),allocatable,save :: ang(:,:)
 logical,save :: Cond_Inicial_ordenado
 character(len=7),save :: Cond_Contorno
 !# Interfaces que dependen de las condiciones de contorno #!
 pointer :: deltE,helicidad
 interface
   function deltE(x,y,na)
    real(kind(0d0)) :: deltE
    integer,intent(in) :: x,y
    real(kind(0d0)),intent(in),optional :: na
   end function deltE
   subroutine helicidad(HT)
    real(kind(0d0)),intent(in) :: HT
   end subroutine helicidad
 end interface
contains

!###----------- Subrutina de Inicialización ------------###!
 subroutine inicial
  implicit none
  integer :: x,y

 !-- Configuración inicial
  do y=1,L
    do x=1,L
      if (Cond_Inicial_ordenado) then
        ang(x,y) = 0.d0 		
      else
        ang(x,y) = 2.d0*pi*grnd()
      endif
    enddo
  enddo

 !-- Energía inicial
  E=0.0d0
  do y=1,L
    do x=1,L
      E = E + deltE(x,y)
    enddo
  enddo
  E = E*0.5d0 !Eliminando conteo Múltiple

 !paso_acep=0
 N_muestreo=0
 sum_E  = 0.0d0
 sum_E2 = 0.0d0
 sum_MH = 0.0d0
 sum_MH2= 0.d0         !!
 sum_EMH= 0.d0         !!
 end subroutine inicial

!###------ CC Periodica: Subrutina Delta Energía ------###!
 function deltE_CC_Periodica(x,y,nuevo_ang)
  implicit none
  real(kind(0d0)) :: deltE_CC_Periodica
  integer :: xp1,xm1,yp1,ym1
  integer,intent(in) :: x,y
  real(kind(0d0)),intent(in),optional :: nuevo_ang
 
  xp1=1+mod(x,L)
  xm1=1+mod(L+x-2,L)
  yp1=1+mod(y,L)
  ym1=1+mod(L+y-2,L)
  if ( present(nuevo_ang) ) then
    deltE_CC_Periodica = 			 &
	    - dcos( nuevo_ang - ang(xp1,y) )	 &
	    - dcos( nuevo_ang - ang(xm1,y) ) 	 &
	    - dcos( nuevo_ang - ang(x,yp1) )	 &
	    - dcos( nuevo_ang - ang(x,ym1) )
  else
    deltE_CC_Periodica = 			 &
	    - dcos( ang(x,y) - ang(xp1,y) )	 &
	    - dcos( ang(x,y) - ang(xm1,y) ) 	 &
	    - dcos( ang(x,y) - ang(x,yp1) )	 &
	    - dcos( ang(x,y) - ang(x,ym1) )
  endif
 end function deltE_CC_Periodica

!###------ CC Abierta: Subrutina Delta Energía ---------###!
 function deltE_CC_abierta(x,y,nuevo_ang)
  implicit none
  real(kind(0d0)),target :: deltE_CC_Abierta
  integer,intent(in) :: x,y
  real(kind(0d0)),intent(in),optional :: nuevo_ang
 
  deltE_CC_Abierta = 0.0d0
  if ( present(nuevo_ang) ) then
    if(x .ne. L) deltE_CC_Abierta = deltE_CC_Abierta - dcos(nuevo_ang-ang(x+1,y))
    if(x .ne. 1) deltE_CC_Abierta = deltE_CC_Abierta - dcos(nuevo_ang-ang(x-1,y))
    if(y .ne. L) deltE_CC_Abierta = deltE_CC_Abierta - dcos(nuevo_ang-ang(x,y+1))
    if(y .ne. 1) deltE_CC_Abierta = deltE_CC_Abierta - dcos(nuevo_ang-ang(x,y-1))
  else
    if(x .ne. L) deltE_CC_Abierta = deltE_CC_Abierta - dcos(ang(x,y)-ang(x+1,y))
    if(x .ne. 1) deltE_CC_Abierta = deltE_CC_Abierta - dcos(ang(x,y)-ang(x-1,y))
    if(y .ne. L) deltE_CC_Abierta = deltE_CC_Abierta - dcos(ang(x,y)-ang(x,y+1))
    if(y .ne. 1) deltE_CC_Abierta = deltE_CC_Abierta - dcos(ang(x,y)-ang(x,y-1))
  endif
 end function deltE_CC_Abierta

!###------- CC Periodica: Subrutina de Helicidad --------###!
 subroutine helicidad_CC_Periodica(H_term)
  implicit none
  real(kind(0d0)) :: Ex,H_term !Direccion arbitraria dada
  integer :: x,y,xp1,xm1,yp1,ym1
 
  Ex=0.6d0
  H_term=0.0d0
  do x=1,L
    xp1=1+mod(x,L)
    xm1=1+mod(L+x-2,L)
    do y=1,L
      yp1=1+mod(y,L)
      ym1=1+mod(L+y-2,L)
      H_term = H_term 						&
	+ dsin( ang(x,y)-ang(xp1,y) )*(+1.d0)*Ex		&
	+ dsin( ang(x,y)-ang(xm1,y) )*(-1.d0)*Ex 	 	&
	+ dsin( ang(x,y)-ang(x,yp1) )*(+1.d0)*(1.d0-Ex*Ex)	&
	+ dsin( ang(x,y)-ang(x,ym1) )*(-1.d0)*(1.d0-Ex*Ex)	
    enddo
  enddo
  H_term = H_term*0.5d0
 end subroutine helicidad_CC_Periodica

!###-------- CC Abierta: Subrutina de Helicidad --------###!
 subroutine helicidad_CC_Abierta(H_term)
  implicit none
  real(kind(0d0)) :: Ex,H_term !Direccion arbitraria dada
  integer :: x,y
 
  Ex=0.6d0
  H_term=0.0d0
  do x=1,L
    do y=1,L
      if(x .ne. L) H_term = H_term + dsin(ang(x,y)-ang(x+1,y))*(+1.d0)*Ex
      if(x .ne. 1) H_term = H_term + dsin(ang(x,y)-ang(x-1,y))*(-1.d0)*Ex
      if(y .ne. L) H_term = H_term + dsin(ang(x,y)-ang(x,y+1))*(+1.d0)*(1.d0-Ex*Ex)
      if(y .ne. 1) H_term = H_term + dsin(ang(x,y)-ang(x,y-1))*(-1.d0)*(1.d0-Ex*Ex)
    enddo
  enddo
  H_term = H_term*0.5d0
 end subroutine helicidad_CC_Abierta

!###------------- Subrutina de Metropolis --------------###!
 subroutine metropolis
  implicit none
  integer :: i,x,y
  real(kind(0d0)) :: nang

  do x=1,L
    do y=1,L
      nang = 2.d0*pi*grnd()
      if (grnd() .le. dexp( (deltE(x,y)-deltE(x,y,nang))/T) ) then   
        !paso_acep = paso_acep + 1
        E = E - deltE(x,y) + deltE(x,y,nang)
        ang(x,y) = nang
      endif
    enddo
  enddo
 end subroutine metropolis

!###------------- Subrutina de Muestreo ----------------###!
 subroutine muestreo
  implicit none
  real(kind(0d0)) :: H_term
  integer :: i  

  N_muestreo = N_muestreo + 1

 !-- Promedio E y E^2
  sum_E  = sum_E  + E
  sum_E2 = sum_E2 + E*E

 !-- Promedio MH
  call helicidad(H_term)
  sum_MH  = sum_MH + H_term*H_term
  sum_MH2 = sum_MH2 + H_term**4       !!
  sum_EMH = sum_EMH + E*H_term*H_term !!
 end subroutine muestreo

!###------------- Subrutina de Datos -------------------###!
 subroutine datos
  implicit none
  real(kind(0d0)) :: CE,MH,MH2
  integer :: i

   sum_E  = sum_E /dble(N_muestreo)
   sum_E2 = sum_E2/dble(N_muestreo)
   sum_MH = sum_MH/dble(N_muestreo)

   sum_MH2 = sum_MH2/dble(N_muestreo) !!
   sum_EMH = sum_EMH/dble(N_muestreo) !!

   CE = (sum_E2-sum_E**2)/(T**2)
   MH = sum_E*(-0.5d0) - sum_MH/T

   MH2 = sum_E2*0.25 + sum_MH2/(T*T) + sum_EMH/T !!
   write(34,1133)T,sum_E/dble(L*L),dsqrt(sum_E2-sum_E**2)/dble(L*L),CE/dble(L*L),MH/dble(L*L),dsqrt(MH2-MH**2)/dble(L*L)

  1133 format(10E16.8)
 end subroutine datos

!###- Subrutina q toma instantaneas de las posiciones -###!
 subroutine instantaneas
  implicit none
  integer :: x,y

   instan_cont=instan_cont+1
   do x=1,L
     do y=1,L
       write(instan_cont*10,1139)dble(x),dble(y),ang(x,y),dcos(ang(x,y)),dsin(ang(x,y)),T
     enddo
   enddo

  1139 format(10E16.8)
 end subroutine instantaneas
end module Global_XY

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!
!@@@			Programa			@@@!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!
program XY
use Global_XY
implicit none
 integer :: N_equ,N_temp,i_temp
 real(kind(0d0)) :: t_inicial,t_final

 Cond_Inicial_ordenado=.true.
 Cond_Contorno='periodica'! 'abierta' 'semiperiodica'

 if (Cond_Contorno .eq. 'abierta') then
   deltE     => deltE_CC_Abierta
   helicidad => helicidad_CC_Abierta
 else
   deltE     => deltE_CC_Periodica
   helicidad => helicidad_CC_Periodica
 endif

 L=32
 allocate(ang(L,L)) 

 N_equ=1000
 N_mcs=10000
 interv=70

 t_inicial = 1.d-10
 t_final   = 2.d0  
 N_temp	   = 200

 !---------------- Imprime Info Util ------------------------------!
 write(*,*)
 write(*,'(A50)')'     #############################################'
 write(*,'(A10,A12,A10,A12,A6)')'     #####',' ','MODELO XY',' ',' #####'
 write(*,'(A50)')'     #############################################'
 write(*,'(A10,A8,A31,A)'	     	  )'     # CC: ',Cond_Contorno,' ','#'
 if (Cond_Inicial_ordenado) then
   write(*,'(A10,A9,A30,A)'	     	  )'     # CI: ',' Ordenado',' ','#'
 else
   write(*,'(A10,A12,A27,A)'	     	  )'     # CI: ',' Desordenado',' ','#'
 endif
 write(*,'(A16,I4,A30,A)'	     	  )'     # Tamaño: ',L,'','#'
 write(*,'(A12,F6.2,A7,F6.2,A3,I4,A8,A4,A)')'     # T: de',t_inicial,' hasta ',t_final,' (',N_temp, 'puntos)','#'
 write(*,'(A50)') '     # ----------------------------------------- #'
 write(*,'(A19,A30,A)'		     	  )'     # Monte Carlo:',' ','#'
 write(*,'(A6,A4,A13,A12,A,I6,A7,A)' 	  )'     #',' ','Equilibracion',' ',':',N_equ,' ','#'
 write(*,'(A6,A4,A15,A11,A,I6,A7,A)' 	  )'     #',' ','Nº de pasos MC',' ',':',N_MCS,' ','#'
 write(*,'(A6,A4,A21,A4,A,I6,A7,A)'  	  )'     #',' ','Intervalo de muestreo',' ',':',interv,' ','#'
 write(*,'(A50)')'     #############################################'
 write(*,*)


 do i_temp=1,N_temp
   T=t_inicial+(t_final-t_inicial)*dble(i_temp)/dble(N_temp)

   call inicial
   do paso=1,N_equ
     call metropolis
   enddo
   do paso=1,N_mcs
     call metropolis
     if (mod(paso,interv) .eq. 0 ) then
	call muestreo
	!call instantaneas
     endif
   enddo
   call datos
 enddo

end program XY

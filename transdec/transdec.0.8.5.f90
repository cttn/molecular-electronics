module transdec_global								!matrices 2x2 => que: (ej: 201)
  integer, parameter :: numanillos=401,& !numsitios=numanillos*5+1,	& 	!El numero de anillos -1 debe ser divisible por 4!!!!!!!!!!!
	dimensiones=2,indice=(numanillos-1)*7/2+4,num_conf=1,puntos_ener=1000	!de lo contrario no quedan armadas las matrices (caso contrario puede
  integer,save :: semilla,dist_qui,num_blocks						! terminar en un pC)
  integer,save :: registro(numanillos,num_conf),tolerancia
  complex(8),save :: hop_block(indice,dimensiones,dimensiones),		&
	ener_block(indice,dimensiones,dimensiones),epsi(puntos_ener),	&
	inversa(dimensiones,dimensiones)
  real(8),save :: porc_max, porc_min,eta,porcentaje,dos(puntos_ener)
  real(8),parameter :: ener_lead=0.0d0,hop_lead=5.0d0
  data dos /puntos_ener*0.0d0/
contains
!************************************************************************!
!*	  		Funcion sigma			       		*!
!************************************************************************!
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
end module transdec_global



program transdec
!***********************************************************************!
! Calcula la transmitancia efectiva de macromoleculas bidimensionales.  !
! Versión adaptada para el modelo de bipolarones de PAni-HCl.		    !
!-------------------------------					                    !
! Version 0.0.2 (alfa 2)						                        !
! 20 sep 2009								                            !
! Entendiendo mejor el metodo MCF					                    !
!***********************************************************************!
! v0.5.1 : de ahora en más no uso la formula de thouless		        !
! v0.5.2 : utilizo una funcion para invertir matrices			        !
! v0.5.3 : Deltas y Green diag evaluados				                !
! v0.6.0 : Chequeadas las partes fundamentales, aún no da.		        !
! v0.7.0 : Vuelvo a utilizar la formula de Thouless para Green no diag  !
!		Toma como base la v0.6.0			                        	!
! v0.8.0 : Buscando errores en el Hamiltoniano, ya que las versiones 	!
!		0.6.0 y 0.7.0 dan practicamente igual.			                !
! v0.8.1 : Hamiltoniano corregido y cheq hasta los primeros 21 bloques	!
! v0.8.2 : Contiene chuequeos de las diferentes subrutinas		        !
! v0.8.3 : Probando cadenas simples.					                !
! v0.8.5 : Agrega DoS							                        !
!************************************************************************

use transdec_global
implicit none
  integer :: i_conf,s
  real(8) :: Emax,Emin,Pi


!--- Opciones generales ************************************************!
  semilla=1539111

  porc_max=0.255d0 !Topes para el criterio de aceptacion de configuraciones
  porc_min=0.245d0
  tolerancia= 10000 !Cantidad de pruebas buscando configuraciones con porcentajes adecuados
  dist_qui=1 !Distamcia minima entre quinoides en unidades de bencenos

  Emax= 10.0d0  !Limite superior del barrido en energia
  Emin=-10.0d0  !  "	inferior  " 	"     "    "
  eta=0.00000001d0 !Parte imaginaria para eliminar singularidades


  do s=1,puntos_ener
  epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s)/real(puntos_ener),eta)
  enddo

!--- Asignacion de datos ***********************************************!
  do i_conf=1,num_conf
 	call asignaciones(i_conf)
!	call asignacion_simple(i_conf)

  print *

!--- Seccion de calculos ************************************************!
	call calculos(i_conf)

!--- Datos en Pantalla **************************************************!
  print *
  print *, '****************************************************'
  print *, 'Puntos de energia: 	',puntos_ener
  print *, 'Numero de anillos:  	',numanillos
  print *, 'Porcentaje de Quinoides: 	',porcentaje*100
  print *, 'Numero de matrices: 	',num_blocks
  print *, 'Numero de sitios:   	',num_blocks*2
  print *, '****************************************************'
  print *

  enddo !i_conf

  print *, 'Terminado correctamente'
end program transdec

!******* PRUEBA ****************************************************************!
subroutine asignacion_simple(i_conf)
 use transdec_global
 implicit none
  integer :: i_conf,i,jo

  jo=0
  do i=1,indice
	!-- eners
	ener_block(i,1,1)=dcmplx(1.0d0,0.0d0)
	ener_block(i,2,2)=dcmplx(-1.0d0,0.0d0)
	ener_block(i,1,2)=dcmplx(0.0d0,0.0d0)!!
	ener_block(i,2,1)=dcmplx(0.0d0,0.0d0)!!
	!---hops
	hop_block(i,1,1)=dcmplx(0.0d0,0.0d0)
	hop_block(i,2,2)=dcmplx(0.0d0,0.0d0)
	hop_block(i,1,2)=dcmplx(0.0d0,0.0d0)
	hop_block(i,2,1)=dcmplx(0.0d0,0.0d0)!!
        jo=jo+1
 enddo

 num_blocks=jo

end subroutine asignacion_simple
!******* FIN PRUEBA *********************************************************!



subroutine asignaciones(i_conf)
!--- Asignacion de valores y sorteo
use transdec_global
implicit none
  integer :: i_conf,m,jo,i,j,h,cont,contbenceno,contquinoide,unidad,i_toleran
  real(8) :: ran3
  complex(8) :: hop_Npara,hop_Npara_q,ener_N,ener_Nmas,enerpara,enerorto,	&
	vop_qui,voo_qui,vop_ben,voo_ben


!--- Energias de sitio y hoppings de todos los sitios en BL de PAni-HCl
  hop_Npara=dcmplx(-3.4d0,0.0d0)
  hop_Npara_q=dcmplx(-4.2d0,0.0d0)
  ener_N=dcmplx(-4.2d0,0.0d0)
  ener_Nmas=dcmplx(-4.5d0,0.0d0)

  enerpara=dcmplx(1.1d0,0.0d0)
  enerorto=dcmplx(0.0d0,0.0d0)
  vop_qui=dcmplx(-3.4d0,0.0d0)
  voo_qui=dcmplx(-3.8d0,0.0d0)
  vop_ben=dcmplx(-3.8d0,0.0d0)
  voo_ben=dcmplx(-3.35d0,0.0d0)


!--- Sorteo del orden de los anillos *******************

  i_toleran=0
 8900 continue !si el sorteo no fue satisfactorio

  cont=0
  contbenceno=0
  contquinoide=0

  h=dist_qui
  do j=1,numanillos
        if (h .lt. dist_qui) then !Si h < dist_qui, pongo un benceno
	        registro(j,i_conf)=0 !benceno
        	contbenceno = contbenceno + 1
       	 	cont=cont+1
		h=h+1
	else !Si h >= dist_qui, sorteo
	   unidad = 1 + int(3*ran3(semilla))
           if (unidad .eq. 1) then
		registro(j,i_conf)=1 !quinoide
		contquinoide=contquinoide + 1
		cont=cont+1
                h=0
           else
		registro(j,i_conf)=0 !benceno
	        contbenceno = contbenceno + 1
	        cont=cont+1
		h=h+1
           endif
  	endif
  enddo

! *** Evaluacion del sorteo ******************************
  i_toleran=i_toleran+1
  if (i_toleran .ge. tolerancia) then
    print *
    print *, 'No se encontraron configuraciones que satisfagan ',   &
	'los porcentajes requridos dentro del limite de tolerancia ',&
	'elegido: ',tolerancia
    stop 'bluff'
  endif

  porcentaje=real(contquinoide)/real(cont)

  if (porcentaje .gt. porc_max .or. porcentaje .lt. porc_min) then
	goto 8900
  endif

!--- Asignacion de lugares en base al sorteo *************
   jo=0
   do m=1,numanillos
      if (registro(m,i_conf) .eq. 0) then !benceno
	if (real(m)/2.0d0 .eq. real(m/2) ) then !jo par

	!-- 0 --
	          if ( (m .gt. 1) .and. (registro(m-1,i_conf) .eq. 1) ) then
		    ener_block(jo,2,2)=ener_Nmas
		  else
		    ener_block(jo,2,2)=ener_N
		  endif
		ener_block(jo,2,1)=hop_Npara

		hop_block(jo,2,1)=hop_Npara
		hop_block(jo,2,2)=dcmplx(0.0d0,0.0d0)

	!-- 1 --
		ener_block(jo+1,1,1)=enerpara
		ener_block(jo+1,1,2)=vop_ben
		ener_block(jo+1,2,1)=vop_ben
		ener_block(jo+1,2,2)=enerorto

		hop_block(jo+1,1,1)=vop_ben
		hop_block(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,2)=voo_ben

	!-- 2 --
		ener_block(jo+2,1,1)=enerorto
		ener_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,2)=enerorto

		hop_block(jo+2,1,1)=voo_ben
		hop_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,2)=vop_ben

	!-- 3 --
		ener_block(jo+3,1,1)=enerorto
		ener_block(jo+3,1,2)=vop_ben
		ener_block(jo+3,2,1)=vop_ben
		ener_block(jo+3,2,2)=enerpara

		hop_block(jo+3,1,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,2,1)=hop_Npara
		hop_block(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

	jo = jo+3
	else !jo impar

	!-- 1 --
	          if ( (m .gt. 1) .and. (registro(m-1,i_conf) .eq. 1)) then
		    ener_block(jo+1,1,1)=ener_Nmas
		  else
		    ener_block(jo+1,1,1)=ener_N
		  endif
		ener_block(jo+1,1,2)=hop_Npara
		ener_block(jo+1,2,1)=hop_Npara
		ener_block(jo+1,2,2)=enerpara

		hop_block(jo+1,1,1)=dcmplx(0.0d0,0.0d0)	!siempre los del lado derecho
		hop_block(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,1)=vop_ben
		hop_block(jo+1,2,2)=vop_ben

	!-- 2 --
		ener_block(jo+2,1,1)=enerorto
		ener_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,2)=enerorto

		hop_block(jo+2,1,1)=voo_ben
		hop_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,2)=voo_ben !vop_ben

	!-- 3 --
		ener_block(jo+3,1,1)=enerorto
		ener_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+3,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+3,2,2)=enerorto

		hop_block(jo+3,1,1)=vop_ben
		hop_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,2,1)=vop_ben
		hop_block(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

	!-- 4 --
		ener_block(jo+4,1,1)=enerpara
		ener_block(jo+4,1,2)=hop_Npara

		hop_block(jo+4,1,1)=dcmplx(0.0d0,0.0d0) !vop_ben
		hop_block(jo+4,1,2)=dcmplx(0.0d0,0.0d0)

		if (m .eq. numanillos) then
  	          if ( registro(m-1,i_conf) .eq. 1 ) then
		    ener_block(jo+4,2,2)=ener_Nmas
		  else
		    ener_block(jo+4,2,2)=ener_N
		  endif
		ener_block(jo+4,2,1)=hop_Npara

		hop_block(jo+4,2,1)=hop_Npara
		hop_block(jo+4,2,2)=dcmplx(0.0d0,0.0d0)
                endif

	jo = jo+4
	endif

	contbenceno = contbenceno + 1

  elseif (registro(m,i_conf) .eq. 1) then !quinoide
	if (real(m)/2.0d0 .eq. real(m/2)) then !jo par

	!-- 0 --
	        ener_block(jo,2,2)=ener_Nmas
		ener_block(jo,2,1)=hop_Npara_q

		hop_block(jo,2,1)=hop_Npara_q
		hop_block(jo,2,2)=dcmplx(0.0d0,0.0d0)

	!-- 1 --
		ener_block(jo+1,1,1)=enerpara
		ener_block(jo+1,1,2)=vop_qui
		ener_block(jo+1,2,1)=vop_qui
		ener_block(jo+1,2,2)=enerorto

		hop_block(jo+1,1,1)=vop_qui
		hop_block(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,2)=voo_qui

	!-- 2 --
		ener_block(jo+2,1,1)=enerorto
		ener_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,2)=enerorto

		hop_block(jo+2,1,1)=voo_qui
		hop_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,2)=vop_qui

	!-- 3 --
		ener_block(jo+3,1,1)=enerorto
		ener_block(jo+3,1,2)=vop_qui
		ener_block(jo+3,2,1)=vop_qui
		ener_block(jo+3,2,2)=enerpara

		hop_block(jo+3,1,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,2,1)=hop_Npara_q
		hop_block(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

	jo = jo+3
	else !jo impar

	!-- 1 --
	    	ener_block(jo+1,1,1)=ener_Nmas
		ener_block(jo+1,1,2)=hop_Npara_q
		ener_block(jo+1,2,1)=hop_Npara_q
		ener_block(jo+1,2,2)=enerpara

		hop_block(jo+1,1,1)=dcmplx(0.0d0,0.0d0)	!siempre los del lado derecho
		hop_block(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+1,2,1)=vop_qui
		hop_block(jo+1,2,2)=vop_qui

	!-- 2 --
		ener_block(jo+2,1,1)=enerorto
		ener_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+2,2,2)=enerorto

		hop_block(jo+2,1,1)=voo_qui
		hop_block(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+2,2,2)=vop_qui

	!-- 3 --
		ener_block(jo+3,1,1)=enerorto
		ener_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+3,2,1)=dcmplx(0.0d0,0.0d0)
		ener_block(jo+3,2,2)=enerorto

		hop_block(jo+3,1,1)=vop_qui
		hop_block(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
		hop_block(jo+3,2,1)=vop_qui
		hop_block(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

	!-- 4 --
		ener_block(jo+4,1,1)=enerpara
		ener_block(jo+4,1,2)=hop_Npara_q

		hop_block(jo+4,1,1)=vop_qui
		hop_block(jo+4,1,2)=dcmplx(0.0d0,0.0d0)

	jo = jo+4
	endif

	contquinoide=contquinoide + 1

  else
      print *, 'Existe un error en el registro del sorteo'
      print *, 'Solo debe haber unos y ceros'
      print *, 'Asi no puedo seguir.'
	stop 'bluff'

  endif

  enddo !m

  num_blocks=jo

end subroutine asignaciones



subroutine calculos(i_conf)
use transdec_global
implicit none
  integer :: i_conf,i_dim,j_dim,s,j,i,m,l,n,h
  complex(8) :: matv(dimensiones,dimensiones),mat_hop(dimensiones,dimensiones),	&
	matd(dimensiones,dimensiones),inver(dimensiones,dimensiones),		&
	Green(indice,indice,dimensiones,dimensiones),				&
	mata(dimensiones,dimensiones),mate(dimensiones,dimensiones),		&
	hop_ef(indice,indice,dimensiones,dimensiones),mat_hop_fin(dimensiones,dimensiones)
  real(8) :: ID(dimensiones,dimensiones),TCoher(puntos_ener),pi
  type self
  complex(8) :: mas(0:indice+1,dimensiones,dimensiones),men(0:indice+1,dimensiones,dimensiones)
  end type self
  type (self) delta,sim_delta

  Pi=dacos(-1.0d0)

!***** Matriz identidad **************************
  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
	if (i_dim .eq. j_dim) then
	 ID(i_dim,j_dim)=1.0d0
	else
	 ID(i_dim,j_dim)=0.0d0
  	endif
  enddo
  enddo

!******** DELTAS ***************************************************

  do s=1, puntos_ener   !Loop s. Hay que quebrarlo para quitar del loop partes que no dependen de s*****************************


  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
     if (i_dim .eq. j_dim ) then
       if (i_dim .eq. 1) then
         delta%men(0,1,1)=sigma(s)
         delta%men(1,1,1)=sigma(s)
         delta%mas(num_blocks,1,1)=dcmplx(0.0d0,0.0d0)
         delta%mas(num_blocks,1,1)=dcmplx(0.0d0,0.0d0)
	 sim_delta%mas(num_blocks,1,1)=dcmplx(0.0d0,0.0d0)
       elseif(i_dim .eq. 2) then
         delta%men(0,2,2)=dcmplx(0.0d0,0.0d0)
         delta%men(1,2,2)=dcmplx(0.0d0,0.0d0)
         delta%mas(num_blocks,2,2)=sigma(s)
         delta%mas(num_blocks+1,2,2)=sigma(s)
	 sim_delta%mas(num_blocks,2,2)=sigma(s)
       endif
     else
      delta%men(1,i_dim,j_dim)=dcmplx(0.0d0,0.0d0)
      delta%mas(num_blocks,i_dim,j_dim)=dcmplx(0.0d0,0.0d0)
      sim_delta%mas(num_blocks,i_dim,j_dim)=dcmplx(0.0d0,0.0d0)
     endif
  enddo !i_dim
  enddo !j_dim

!*** Delta L (menos) ************************************************
  do l=2,num_blocks+1

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
        mat_hop(i_dim,j_dim)=hop_block(l-1,i_dim,j_dim)
	inver(i_dim,j_dim)=epsi(s)*ID(i_dim,j_dim)-ener_block(l-1,i_dim,j_dim)-delta%men(l-1,i_dim,j_dim)
  enddo !i_dim
  enddo !j_dim

  call inversora(inver)
  matv=matmul(transpose(mat_hop),inversa)
  matd=matmul(matv,mat_hop)

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
    delta%men(l,i_dim,j_dim)=matd(i_dim,j_dim)
  enddo !i_dim
  enddo !j_dim


  enddo !l

!******** Delta R (mas) **************************************
  do l=2,num_blocks+1

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
        mat_hop(i_dim,j_dim)=hop_block(num_blocks+1-l,i_dim,j_dim)
        inver(i_dim,j_dim)=epsi(s)*ID(i_dim,j_dim)-ener_block(num_blocks+2-l,i_dim,j_dim)	&
		-delta%mas(num_blocks+2-l,i_dim,j_dim)
  enddo !i_dim
  enddo !j_dim

  call inversora(inver)
  matv=matmul(mat_hop,inversa)
  matd=matmul(matv,transpose(mat_hop))

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
    delta%mas(num_blocks+1-l,i_dim,j_dim)=matd(i_dim,j_dim)
    sim_delta%mas(num_blocks+1-l,i_dim,j_dim)=matv(i_dim,j_dim)
  enddo !i_dim
  enddo !j_dim

  enddo !l


!** Funciones de Green ******************************************
!------------ Elementos diagonales -----------------------------!

  do i=1,num_blocks

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
	inver(i_dim,j_dim)=epsi(s)*ID(i_dim,j_dim)-ener_block(i,i_dim,j_dim)			&
		-delta%men(i,i_dim,j_dim)-delta%mas(i,i_dim,j_dim)
  enddo !i_dim
  enddo !j_dim

  call inversora(inver)

  do i_dim=1,dimensiones
  do j_dim=1,dimensiones
     Green(i,i,i_dim,j_dim)=inversa(i_dim,j_dim)
  enddo
  enddo

  enddo !i

!***************************************************************!
!******************** DoS **************************************!
  do l=1,num_blocks						!
  do i_dim=1,dimensiones					!
  dos(s)=dos(s)-(1.0d0/pi)*aimag(Green(l,l,i_dim,i_dim))	!
!  write(201,*)real(epsi(s)),(-1.0d0/pi)*aimag(Green(1,1,1,1))
!  write(250,*)real(epsi(s)),(-1.0d0/pi)*aimag(Green(50,50,1,1))
  enddo								!
  enddo								!
!****************** Fin DoS ************************************!
!***************************************************************!

!------ Elem no diagonales de la func de Green ---------------
 !--
  do i=2,num_blocks
  do l=1,i-1

    do i_dim=1,dimensiones
    do j_dim=1,dimensiones
	mat_hop(i_dim,j_dim)=Green(l,i-1,i_dim,j_dim)
	mat_hop_fin(i_dim,j_dim)=sim_delta%mas(i-1,i_dim,j_dim)
    enddo
    enddo

    matv=matmul(mat_hop,mat_hop_fin)

    do i_dim=1,dimensiones
    do j_dim=1,dimensiones
       Green(l,i,i_dim,j_dim)=matv(i_dim,j_dim)
    enddo
    enddo

  enddo !l
  enddo !i


! Info de avanze de esta subrutina
   do n=20,100,20
   if (n .eq. int(real(s)*100./real(puntos_ener))) then
        if (h .eq. n) then
		goto 118
        else
		print *, n,'% completado'
        	h=n
!******* Lo que sigue es totalmente innecesario ******************************************!*
!		if (n .eq. 20) then							  !*
!		call system ('notify-send -u critical Transdec "20 % completado"')   	  !*
!		elseif (n .eq. 40) then							  !*
!		call system ('notify-send -u critical Transdec "40 % completado"')   	  !*
!		elseif (n .eq. 60) then							  !*
!		call system ('notify-send -u critical Transdec "60 % completado"')   	  !*
!		elseif (n .eq. 80) then							  !*
!		call system ('notify-send -u critical Transdec "80 % completado"')   	  !*
!		endif									  !*
!******* Lo anterior es totalmente innecesario *******************************************!*
        endif
   endif
   enddo
 118 continue

!***** EVALUACIONES! *************************************
!********** Test Trans Coher *******************
  TCoher(s)=4.0d0*(dabs(aimag(sigma(s)))**2)*(cdabs(Green(1,num_blocks,1,2)))**2
  write(113,*)4.0d0*(dabs(aimag(sigma(s)))**2)*(cdabs(Green(1,num_blocks-4,1,2)))**2
!***********************************************
  enddo !s

!********** Test trans coher CONT ****************
   print *
   print *, 'Test Transmitancia Coherente'
   open(unit=189,file='transcoher.tst')
   do s=1,puntos_ener
   write(189,*)real(epsi(s)),TCoher(s),dos(s)
   enddo
   close(189)
!**********************************************
end subroutine calculos

!=========================================================!
! 		FUNCIONES UTILES			  !
!=========================================================!

!*********************************************************
!*			Inversora			 *
!*********************************************************
subroutine inversora(matriz)
 use transdec_global
  integer :: i,j
  complex(8) :: matriz(dimensiones,dimensiones),deter

    do i=1,dimensiones
    do j=1,dimensiones
	if (i .ne. j)then
	  inversa(i,j)=-matriz(j,i)
	else
	  if (i .eq. 1) then
	   inversa(i,i)=matriz(2,2)
	  elseif (i .eq. 2) then
           inversa(i,i)=matriz(1,1)
	  endif
	endif
    enddo
    enddo
  deter=inversa(1,1)*inversa(2,2)-inversa(1,2)*inversa(2,1)
  inversa=inversa/deter

end subroutine inversora

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

!************************************************************************
!   Idem que icrdm-ngX.f90 pero se le agrega el calculo de la 		!
!                 transmitancia incoherente				!
!-----------------------------------------------------------------------!
!        	  It works! - 11-nov-2008				!				
!       Arreglé sigma, ahora S=D-iG (en lugar de S=D+iG)		!
!    Corregí las T's, ahora dependen del abs(aimag(sigma))		!
!-----------------------------------------------------------------------!
!     A diferencia de icrdm-ng6, aquí se promedian diferentes		! 
!                      configuraciones					!
!-----------------------------------------------------------------------!
!  A diferencia de icrdm-ng7, aquí corren varios Gama automáticamente   !
!       Genera un archivo por cada Gama. En cada Archivo está la 	!
!  			transmitancia promedio.				!
!-----------------------------------------------------------------------!
! Corregido el calculo del nivel de Fermi.				!
!-----------------------------------------------------------------------!
! Corregido el traumático error de indices en el calculo de la T 	!
!  efectiva								!
!************************************************************************
! Version ng19:								!
! En principio he corregido todos los errores que tenia TEFLIN.FOR,	!
! el prgrama de D'Amato.						!
!-------------------------						!
! Agrega la impresion de resultados de T Coherente			!
!************************************************************************

module globalesrdm
  use diag_icrdm
  use decbipol
  integer,parameter :: indice=numanillos*3+1
  integer :: s,j,i,m,porcentaje,puntos_decbipol
  integer, save :: numsitios,contbenceno,contquinoide,	 &
	hhh,registro(numanillos,num_conf),i_conf
  complex(8),save :: e1q(puntos_ener),v12q(puntos_ener),e1b(puntos_ener),&
	v12b(puntos_ener),ener_sitio(0:indice+1,puntos_ener),		 &
	dos(puntos_ener),epsi(puntos_ener),ener_N,ener_Nmas,hop_Npara_q, &
	hop_Npara,ener_renor(indice,puntos_ener),Gama_dec,		 &
	hop_sitio(0:indice,puntos_ener)!,sigma(puntos_ener)
  real(8), parameter :: Pi=dacos(-1.0d0),ener_lead=0.0d0,hop_lead=5.0d0, &
      tun_left=5.0d0,tun_right=5.0d0 
  real(8),save :: eta,tef(puntos_ener),alpha,Fermi(num_conf),tcoher(puntos_ener)
  character(len=13) :: entrada
  character(len=30) :: resultado,configuraciones
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
  function moder(complejo)
     complex(8) :: complejo
     real(8) :: moder
     moder = sqrt((real(complejo))**2 + (aimag(complejo))**2)  
  end function moder
end module globalesrdm

module gamma_par
!---Parametros de Gamma ----!
  integer,parameter :: N_alpha=1
  real(8),parameter :: kT=0.025d0, alpha_min= 2.3d0*kT,alpha_max= 2.3d0*kT
!---Fin parametros Gamma ---!
end module gamma_par



program icrdm_ng19
  use diag_icrdm
  use globalesrdm
  use gamma_par
  implicit none
  integer :: mm,hh,nn,i_alpha
  real(8) :: tef_final(puntos_ener,num_conf),varA,tcoher_final(puntos_ener,num_conf)


!-- Variables manuales de MAIN__
  entrada='decbipol.inp'
  configuraciones='config-ng19.inp'

!-- Leyendo archivo de decimaciones
  call lectora(entrada)

  do i_conf=1,num_conf

! Comienza el loop en los Gama's

  varA= (alpha_max - alpha_min + 1.0d0)**(1.0d0/real(N_alpha))

!MOD****************************************
  do i_alpha=1,N_alpha 
  alpha=varA**i_alpha + alpha_min-1.0d0


  !OJO: iconf modificado para que empieze en inc!!!!!!!!!!!!!!!!!!!!

  resultado='icrdm-ng19_'// char(48+int((i_conf+inc)/100))//		 &
	char(48+int((i_conf+inc)/10)-10*int(int((i_conf+inc)/10)/10))//	 &
	char(48+i_conf+inc-10*int((i_conf+inc)/10))//'p_'// 		 &
	char(48+int((i_alpha)/1000))//					 &
        char(48+int((i_alpha)/100)-10*int((i_alpha)/1000))//	 	 &
	char(48+int((i_alpha)/10)-10*int(int((i_alpha)/10)/10))//	 &
	char(48+i_alpha-10*int((i_alpha)/10))//'G.dat'

  Gama_dec=dcmplx(0.0d0,-alpha)   


!-- Empezando a calcular configuraciones
  print *,'-------------------------------------------------------------------'
  print *, 'Comenzando a calcular configuracion: ',i_conf+inc
  print *, 'Gama: ', aimag(Gama_dec)


!-- Sorteando sitios
  call mezcladora
 

!-- Nucleo de calculos
  call calculadora


!-- Sumando nuevos resultados
  do s=1,puntos_ener
  tef_final(s,i_conf)=tef(s)
  tcoher_final(s,i_conf)=tcoher(s)

  if (dabs(tef_final(s,i_conf)) .lt. 1.0d-99) then
    tef_final(s,i_conf)=0.0d0
  endif
  if (dabs(tcoher_final(s,i_conf)) .lt. 1.0d-99) then
    tcoher_final(s,i_conf)=0.0d0
  endif

  enddo !s


!-- Impresion de resultados
  print *
  print *, '************ Concluído ********************'
  print *, 'Conductancia Incoherente - Bipolarones: '
  print *
  print *, 'Puntos de energia:        ',puntos_ener
  print *, 'Anillos:                  ',numanillos
  print *, 'Total de sitios:          ',numsitios
  print *, '# de anillos bencenicos : ', contbenceno
  print *, '# de anillos quinoides:   ', contquinoide
        porcentaje = (real(contquinoide)/real(numanillos))*100
  print *, 'Hay un ', porcentaje , ' % de anillos quinoides'
  print *, 'Archivo de entrada:       ',entrada
  print *, 'Archivo de salida:        ',resultado
  print *, '*******************************************'
  print *  

  open(unit=44, file=resultado)

  do s=1,puntos_ener
  write(44,1011) real(epsi(s)),tef_final(s,i_conf)!,tcoher_final(s,i_conf)
  enddo

  close(44)

 1011 format (6E16.8)

  enddo ! loop en los Gama's

  enddo ! i_conf
  print *,'Terminado'
end program icrdm_ng19



subroutine lectora(archivo)
  use diag_icrdm
  use globalesrdm
  implicit none
  real(8) :: real_hop_Npara,imag_hop_Npara,                              &
     real_ener_N,imag_ener_N,real_ener_Nmas,imag_ener_Nmas,              &
     real_hop_Npara_q,imag_hop_Npara_q,imag_v12b,       	         &
     epsi_inp,real_e1q,imag_e1q,real_v12q,imag_v12q,                     &
     real_e1b,imag_e1b,real_v12b
  character (len=30) archivo, archivo2

  archivo=entrada
  archivo2=configuraciones
  open(unit=75,file=archivo, status='old')

  read(75,1077)puntos_decbipol
     if (puntos_decbipol .ne. puntos_ener) then
		print *, 'Los rangos de energía no coinciden'
  		print *, 'Es necesario correr nuevamente decbipol.f90'
		stop
     endif
  read(75,1076)eta,real_hop_Npara,imag_hop_Npara,                        &
     real_ener_N,imag_ener_N,real_ener_Nmas,imag_ener_Nmas,              &
     real_hop_Npara_q,imag_hop_Npara_q

!**************************************************************
!  print * 
!  print *, 'MOD: Partes imaginarias forzadas a cero'
!  print *
!******* MOD .. Partes imaginarias forzadas a cero ************
!  ener_N=dcmplx(real_ener_N,0.0d0)
!  hop_Npara=dcmplx(real_hop_Npara,0.0d0)
!  hop_Npara_q=dcmplx(real_hop_Npara_q,0.0d0)
!  ener_Nmas=dcmplx(real_ener_Nmas,0.0d0)
  ener_N=dcmplx(real_ener_N,imag_ener_N)
  hop_Npara=dcmplx(real_hop_Npara,imag_hop_Npara)
  hop_Npara_q=dcmplx(real_hop_Npara_q,imag_hop_Npara_q)
  ener_Nmas=dcmplx(real_ener_Nmas,imag_ener_Nmas)


  do s=1,puntos_ener
  read(75,1075)epsi_inp,real_e1q,imag_e1q,real_v12q,imag_v12q,           &
      real_e1b,imag_e1b,real_v12b,imag_v12b

!**********************
!	eta=0.0d0
!**********************

!**************************************************************
!******* MOD .. Partes imaginarias forzadas a cero ************
!  epsi(s)=dcmplx(epsi_inp,0.0d0)
!  e1q(s)=dcmplx(real_e1q,0.0d0)
!  v12q(s)=dcmplx(real_v12q,0.0d0)
!  e1b(s)=dcmplx(real_e1b,0.0d0)
!  v12b(s)=dcmplx(real_v12b,0.0d0)
  epsi(s)=dcmplx(epsi_inp,-eta)
  e1q(s)=dcmplx(real_e1q,imag_e1q)
  v12q(s)=dcmplx(real_v12q,imag_v12q)
  e1b(s)=dcmplx(real_e1b,imag_e1b)
  v12b(s)=dcmplx(real_v12b,imag_v12b)


  enddo
 1075   format(11E16.8)
 1076   format(9E16.8)
 1077   format(1I5)
  close(75)

  print *
  print *, 'Archivo de datos leido'


!*** Archivo de configuraciones ****************
  open(57, file=archivo2, status='old')
  read(57,1056) (Fermi(i), i=1,num_conf)
  do m=1,numanillos
  read(57,1055)(registro(m,i),i=1,num_conf)
  enddo
  close(57)
 1056 format(100E16.8)
 1055 format(100I4)
 
end subroutine lectora



subroutine mezcladora
  use diag_icrdm
  use globalesrdm
  implicit none
  integer ::  jo,ho,io
  real(8) :: re_ener,im_ener,re_hop,re_sigma,im_sigma,im_gama


  contbenceno=0
  contquinoide=0
  jo=1   !j corre sobre los sitios 

  do m=1,numanillos

    if (registro(m,i_conf) .eq. 0) then !benceno:
  	do s=1,puntos_ener
		  if (registro(m-1,i_conf) .eq. 1) then
		    ener_sitio(jo,s)=ener_Nmas+Gama_dec
		  else
		    ener_sitio(jo,s)=ener_N +Gama_dec
		  endif
		hop_sitio(jo,s)=hop_Npara
		ener_sitio(jo+1,s)=e1b(s)+Gama_dec
		hop_sitio(jo+1,s)=v12b(s)
		ener_sitio(jo+2,s)=e1b(s)+Gama_dec
		hop_sitio(jo+2,s)=hop_Npara
	enddo
	jo=jo+3
        contbenceno = contbenceno + 1
   else 
	do s=1,puntos_ener
		ener_sitio(jo,s)=ener_Nmas +Gama_dec
		hop_sitio(jo,s)=hop_Npara_q
		ener_sitio(jo+1,s)=e1q(s)+Gama_dec
		hop_sitio(jo+1,s)=v12q(s)
		ener_sitio(jo+2,s)=e1q(s)+Gama_dec
		hop_sitio(jo+2,s)=hop_Npara_q
	enddo
       	jo=jo+3
	contquinoide=contquinoide + 1
  endif
  enddo !m


  do s=1,puntos_ener
  if (real(jo/2) .eq. real(jo)/2.0d0 ) then  !anda bien chequeado
  ener_sitio(jo,s)=ener_Nmas+Gama_dec
  else
  ener_sitio(jo,s)=ener_N+Gama_dec
  end if
  end do

  numsitios=jo

! Agregando los extremos: ------------->!En principio no necesito esto
  do s=1,puntos_ener
  hop_sitio(0,s)= -tun_left
  hop_sitio(numsitios,s)= -tun_right
  ener_sitio(0,s)= -ener_lead
  ener_sitio(numsitios+1,s)= -ener_lead
  enddo

end subroutine mezcladora



subroutine calculadora
  use diag_icrdm
  use globalesrdm
  implicit none
  integer :: n,h !sólo para la info de avanze de la subrutina
  complex(8) :: Green(0:indice+1,0:indice+1)
  real(8) :: trans(0:indice+1,0:indice+1)
  type self
  complex(8) :: mas(0:indice+1,puntos_ener),men(0:indice+1,puntos_ener)
  end type self
  type (self) delta


!  print *, 'Calculando transmitancias...'

!-- Correccion a energia de sitio por renormalizacion
  do s=1,puntos_ener
      delta%men(0,s)=sigma(s)  !Sigma es la Self-energy debida al lead.
      delta%men(1,s)=sigma(s)
      delta%mas(numsitios,s)=sigma(s)
      delta%mas(numsitios+1,s)=sigma(s)

  do j=2,numsitios+1
      delta%men(j,s)= (hop_sitio(j-1,s)*hop_sitio(j-1,s))                &
      /(epsi(s)-ener_sitio(j-1,s)-delta%men(j-1,s))

      delta%mas(numsitios+1-j,s)=(hop_sitio(numsitios+1-j,s)             &
      *hop_sitio(numsitios +1-j,s)/(epsi(s)-ener_sitio(numsitios+2-j,s)  &
      -delta%mas(numsitios+2-j,s)))

  enddo !j


!-- Energias renormalizadas
  do j=1,numsitios
      ener_renor(j,s)=ener_sitio(j,s)+delta%mas(j,s)+delta%men(j,s)
!      dos(s)=dos(s)+(-1.0d0/pi)*aimag(1.0d0/(epsi(s)-ener_renor(j,s)))
  enddo !j
  enddo !s

!-- Funciones de Green sitio a sitio
  do s=1,puntos_ener

      do i=1,numsitios
      Green(i,i)=1.0d0/(epsi(s)-ener_renor(i,s))
      enddo

      do i=2,numsitios+1 !Correccion al D'Amato
      do j=1,i-1
      Green(j,i)=Green(j,i-1)*delta%mas(i-1,s)/hop_sitio(i-1,s)
      enddo
      enddo

!*********** Mi correción al programa de D'Amato ******************
      Green(0,0)=Green(1,1)
      Green(0,numsitios+1)=Green(1,numsitios)
      Green(numsitios+1,numsitios+1)=Green(numsitios,numsitios)
      do i=1,numsitios
      Green(0,i)=Green(1,i)
      Green(numsitios+1,i)=Green(numsitios,i)
      enddo


!Transmitancias

! Transmitancia coherente:
  trans(0,numsitios+1)=4.0d0*(cdabs(green(0,numsitios+1))**2)*aimag(sigma(s))**2 !No la necesito

!********* Si G=0 saltear D'amato-Pastawski ********************
   if (aimag(Gama_dec) .eq. 0.0d0 ) then
           goto 1245
   endif 
!***************************************************************

  trans(0,0)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(sigma(s)))*green(0,0)-1.0d0))**2 

  trans(numsitios+1,numsitios+1)=cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*        &  ! No la necesito
          dabs(aimag(sigma(s)))*green(numsitios+1,numsitios+1)-1.0d0)**2    !

  do i=1,numsitios
  trans(0,i)=4.0d0*dabs(aimag(sigma(s)))*(dabs(aimag(Gama_dec)))*((cdabs(green(0,i)))**2)
  trans(i,i)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(Gama_dec))*green(i,i)-1.0d0))**2 
  trans(i,numsitios+1)=4.0d0*dabs(aimag(sigma(s)))*dabs(aimag(Gama_dec))*(cdabs(green(i,numsitios+1)))**2 !No las necesito
  enddo
  

  do i=1,numsitios-1
  do j=i+1,numsitios
  trans(i,j)=4.0d0*((aimag(Gama_dec))**2)*(cdabs(green(i,j))**2)
  enddo
  enddo


  do i=numsitios,1,-1
  do j=0,i-1
  do m= j, i-1
    
    trans(j,m)=trans(j,m)+trans(j,i)*trans(m,i)/(1.0d0-trans(i,i))

  enddo
  enddo
  enddo

! Tansmitancia total efectiva y transmitancia coherente:

 1245 continue ! Si Gama = 0 saltear D'Amato-Pastawski

  if (aimag(Gama_dec) .eq. 0.0d0) then
	  tef(s)=trans(0,numsitios+1) 
  else
  	tef(s)=1.0d0-trans(0,0) 
  	tcoher(s)=trans(0,numsitios+1)
  endif

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
 
end subroutine calculadora

program cor_icrdm01
!***********************************************************************!
!   Pertenece al grupo de programas para calcular las propiedades 	!
! electricas de cadenas tight-binding.					!
!-----------------------------------------------------------------------!
!   Integra la transmitancia efectiva multiplicada por la diferencia	!
! entre las funciones de fermi (de los electrodos) Left y Right		!
! para obtener la corriente a un potencial dado.			!
!-----------------------------------------------------------------------!
!   Necesita los datos de Transmitancia efectiva (y los puntos de 	!
! energia) que se obtienen mediante el programa icrdm-ngXX.f90		!
!***********************************************************************!
implicit none
  integer, parameter :: max_conf=50,puntos=1024,max_volt_pun=50
  integer :: unid,cant_conf,puntos_ener,j,m,punto_inicial,punto_final,	&
	volt_pun,k,inc,i_alpha
  real(8) :: energia(max_conf,puntos),transmitancia(max_conf,puntos), 	&
	voltaje(max_volt_pun),corriente(max_conf,max_volt_pun),		&
	const_ini,diflr,prom_corr(max_volt_pun),volt_max,volt_min
  real :: alpha_cor
  character :: archivo*30,resultado_prom*18,resultados*30
  data corriente /2500*0.0/ !esta porqueria... ver como ponerlo!!!!!!!!!!!!!!!
  data prom_corr /max_volt_pun*0.0/


!-- Variables locales **********************************
  cant_conf=50  !Numero de configuraciones a examinar.
!  alpha_cor=0.901 !OJO: esta variable determina el nombre de los archivos a leer
  puntos_ener=200  !OJO: este dato viene de icrdm CUIDADO!!!
  const_ini=-1.0d0  !Constante inicial 2e/h
!  const_ini=
  punto_inicial= 1 !=(-2)  !indice del punto en el q comienza la ventana de integracion
  punto_final= 200 !=(+2)   !indice del punto en el q termina la ventana de integracion
  resultado_prom='cor-icrdm-ng19.dat'

  inc=0
  i_alpha=1

!-- Leyendo archivos de transmitancias ******************
  do j=1,cant_conf

    unid=j+10
    archivo='icrdm-ng19_'// char(48+int((j+inc)/100))//		 &
	char(48+int((j+inc)/10)-10*int(int((j+inc)/10)/10))//	 &
	char(48+j+inc-10*int((j+inc)/10))//'p_'// 		 &
	char(48+int((i_alpha)/1000))//					 &
        char(48+int((i_alpha)/100)-10*int((i_alpha)/1000))//	 	 &
	char(48+int((i_alpha)/10)-10*int(int((i_alpha)/10)/10))//	 &
	char(48+i_alpha-10*int((i_alpha)/10))//'G.dat'

    open (unit=unid, file=archivo, status='old')

      do m=1,puntos_ener
        read(unid,1000) energia(j,m), transmitancia(j,m) !j da el numero de conf y m la funcion.

      enddo

    close(unid)

!    print *, energia(j,1),transmitancia(j,1)
  enddo


!-- Definiendo el voltaje  ***********************
  volt_max= 1.5d0
  volt_min=-1.5d0
  volt_pun=50  !Numero de puntos de medicion de voltaje
  do k=1,volt_pun
  voltaje(k)=(real(k)*(volt_max-volt_min))/real(volt_pun)+volt_min
  enddo


!-- Integrando... ********************************

  do k=1,volt_pun
  do j=1,cant_conf
  do m=punto_inicial,punto_final-2
  
  corriente(j,k) = corriente(j,k) + transmitancia(j,m+1)*		 & 
	diflr(energia(j,m+1),voltaje(k))*(energia(j,m+2)-energia(j,m))


  enddo
  enddo
  enddo

  corriente=corriente*const_ini

!-- Promedio de las corrientes  *******************
  do k=1,volt_pun
  do j=1,cant_conf
  prom_corr(k)=prom_corr(k)+corriente(j,k)
  enddo
  enddo
  prom_corr=prom_corr/real(cant_conf)


!-- Escritura de Resultados  **********************

  open(unit=121,file=resultado_prom)

  do k=1,volt_pun
  write(121,1121)voltaje(k),prom_corr(k)
  enddo

  close(121)


  do j=1,cant_conf
    resultados='cor-icrdm-ng19_'// char(48+int( j/100))//		 &
	char(48+int( j/10)-10*int(int( j/10)/10))//			 &
	char(48+ j-10*int( j/10)) // '.dat'

    unid=j*10

    open(unit=unid, file=resultados)

    do k=1,volt_pun
    write(unid,1111) voltaje(k),corriente(j,k)
    enddo

    close(unid)
  enddo


  open (unit=99999, file='dif_fermi.dat')
  do j=punto_inicial,punto_final-1
  write(99999,1119) energia(1,j),(diflr(energia(1,j),voltaje(k)),k=1,volt_pun)
  enddo
  close(99999)
  
  open (unit=2999, file='voltajes.dat')
  do k=1,volt_pun
  write(2999,*)k,voltaje(k)
  enddo
  close(2999)


  print *
  print *,'******************************************************'
  print *,'Resultados de corriente (cor-icrdm):'
  print *
  print *,'Ventana de energia entre ',energia(1,punto_inicial),' y ', energia(1,punto_final)
  print *,'Ventana de voltaje entre ',volt_min, ' y ',volt_max
  print *,'Promedio final en: ',resultado_prom
  print *,'******************************************************'
  print *

 1119 format(201E16.8)
 1111 format(2E16.8)
 1121 format(2E16.8)
 1000 format(1X,2E16.8)
end program cor_icrdm01


!***********************************************************************!
!		FUNCIONES Y SUBRUTINAS					!
!***********************************************************************!

function diflr(epsi,volt)
!Diferencia entre las funciones de fermi de los leads left y right
 implicit none
 real(8), parameter :: kt=0.025d0 !eV
 real(8) :: diflr,epsi,volt,asim

	asim=0.5d0
  diflr = 1.0d0/(exp( (epsi+(1.0d0- asim)*volt)/kt ) +1.0d0)-1.0d0/(exp((epsi-asim*volt)/kt) +1.0d0)

end function diflr


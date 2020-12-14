program histotef_ng19
!*************************************************************************
! V01:									 *
! Calcula histogramas a partir de diagbipol-ng14 y icrdm-ng14		 *
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*
! V02:									 *
! Ademas de los histogramas realiza calculos de probabilidad acumulativa *
! para los valores de transmitancia					 *
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*
! NG16:									 *
! Versión adaptada para toda la suite de programas ng16 para calculos de *
! trasmitancia efectiva.
!*************************************************************************
  use decbipol
  use diag_icrdm
  use histo_mod_G
  use gamma_par
  implicit none
  integer, parameter :: max_barras=1000!,num_conf=100
  character :: fermi_input*15,archivo*30,resultado*17,acumulativo*15,	 &
	datos*20,cabec_nom*11,Largo_nom,efectivas*21
  integer :: ios,j,unid,m,dummy1,maximo,minimo,num_barras,i,ind,	 &
	hist_num_conf  !,inc,num_conf
  real(8) :: Efermi(num_conf),ener(puntos_ener),trans(puntos_ener),	 &
	dif(puntos_ener),trans_ef(num_conf),sum_trans(0:max_barras),	 &
	x(0:max_barras),ancho,alpha,promedio,desv,desv_cuad,prev,post
  data sum_trans /max_barras*0.0d0,0.0d0/

!-- Constantes utiles ****************************
  fermi_input='config-ng19.inp'
  resultado='histotef-ng19.dat'
  acumulativo='acumul-ng19.dat'
  datos='datos_hist-02.dat'
  efectivas='histotef-ng19_tef.dat'
  Largo_nom='C' !  Si los archivos a leer estan numerados con 3 cifras, elegir C;
                ! pero si estan numerados con 4 cifras, elegir M.
  cabec_nom='icrdm-ng19_' !  Cabecera del nombre de los archivos a leer, 
			  ! indica el programa del cual vienen. OJO: ajustar el tamaño
			  ! del nombre en la declaracion de variables!!!!
!  inc=0
  hist_num_conf=num_conf-inc
!  alpha=7.001
  ancho=0.05d0
  minimo=0.0d0
  maximo=1.0d0

!-- Leyendo energías de Fermi *********************
  open (unit=197, file=fermi_input, iostat=ios, status='old')
    read(197,1197)(Efermi(j), j=1,hist_num_conf)
  close(197)

!-- Leyendo archivos de transmitancias ************
  do j=1,hist_num_conf

    unid=j+10

    if (Largo_nom .eq. 'C') then
      archivo=cabec_nom// char(48+int(inc+j/100))//			 &
	char(48+int(inc+j/10)-10*int(int( j/10)/10))//			 &
	char(48+inc+j-10*int(inc+j/10))//'p_'//		 		 &
	char(48+int((i_alpha)/1000))//					 &
        char(48+int((i_alpha)/100)-10*int((i_alpha)/1000))//	 	 &
	char(48+int((i_alpha)/10)-10*int(int((i_alpha)/10)/10))//	 &
	char(48+i_alpha-10*int((i_alpha)/10))//'G.dat'

!******* 3 cifras para alpha *********************************************
!	char(48+int((i_alpha)/100))//		 			 &
!	char(48+int((i_alpha)/10)-10*int(int((i_alpha)/10)/10))//	 &
!	char(48+i_alpha-10*int((i_alpha)/10))//'G.dat'

    elseif (Largo_nom .eq. 'M') then

      archivo=cabec_nom// char(48+int((j+inc)/1000))//			 &
        char(48+int((j+inc)/100)-10*int((j+inc)/1000))//	 	 &
	char(48+int((j+inc)/10)-10*int(int((j+inc)/10)/10))//	 	 &
	char(48+j+inc-10*int((j+inc)/10))//'p_'//		 	 &
	char(48+int((i_alpha)/1000))//					 &
        char(48+int((i_alpha)/100)-10*int((i_alpha)/1000))//	 	 &
	char(48+int((i_alpha)/10)-10*int(int((i_alpha)/10)/10))//	 &
	char(48+i_alpha-10*int((i_alpha)/10))//'G.dat'
     endif
!    print *, archivo

    open (unit=unid, file=archivo, status='old')

      do m=1,puntos_ener
        read(unid,1000) ener(m), trans(m)
        dif(m)=dabs(Efermi(j)-ener(m))
      enddo

      trans_ef(j)= trans(minloc(dif,1))

    close(unid)

!  print *,trans_ef(j),j
    
  enddo
  
 
!-- Histograma con promedio ***********************
  num_barras=int((maximo-minimo)/ancho)
  do i=0,num_barras
     do j=1,hist_num_conf
        if((trans_ef(j).ge.(minimo+ancho*(real(i)-0.5d0)).and.		 &
  	   (trans_ef(j) .lt. (minimo+ancho*(real(i)+0.5d0))))) then
           sum_trans(i)=sum_trans(i)+1 
        endif
     enddo !fin loop j
  enddo !fin loop i

  sum_trans=100.*sum_trans/real(hist_num_conf)
!  print*,sum(sum_trans),num_barras


  if (num_barras .gt. max_barras) then
     print *, '**** Atención POR FAVOR!!! *******'
     print *, 'El # de beams es muy grande, aumentar el parametro max_pics'
     print *, 'Estimado: ',num_barras+1
  endif

  do i=0,num_barras
      x(i)=minimo + ancho*(real(i))
  enddo


!-- Promedio y desviacion cuadratica  ************
  promedio=sum(trans_ef)/real(hist_num_conf)

  desv=0.0d0
  do i=1,hist_num_conf
  desv=desv+(trans_ef(i)-promedio)**2
!  print*, desv,trans_ef(i),promedio
  enddo

  desv_cuad=sqrt(desv/real(hist_num_conf))


!-- Probabilidad acumulativa  *********************
  do i=1,hist_num_conf

  prev=maxval(trans_ef(i:hist_num_conf),1)
  ind=maxloc(trans_ef(i:hist_num_conf),1)
  post=trans_ef(i)

  trans_ef(i)=prev
  trans_ef(ind+i-1)=post

!  print*, trans_ef(i)
  enddo

 
  print*
  print*, '***************  Promedios  *******************'  
  print*, 'Promedio = ',promedio
  print*, 'Desviacion Cuadratica Media= ',desv_cuad
  print*, 'Numero de configuraciones= ',hist_num_conf
  print*, 'Histograma en: ',resultado
  print*, 'Acumulativo en: ',acumulativo
  print*, '***********************************************'
  print*

!-- Escritura de resultados  **********************
  open(unit=200, file = resultado)
      do i=0,num_barras
         write(200,1200)x(i),sum_trans(i)
      enddo
  close(200)

  open(unit=202, file=acumulativo)
      do i=1,hist_num_conf
	write(202,1202)trans_ef(i), real(i)/real(hist_num_conf)
      enddo
  close(202)

  open(unit=204, file=datos)
      do i=1,hist_num_conf
        write(204,1204)trans_ef(i)
      enddo
  close(204)

  open(unit=206, file=efectivas)
      write(206,1206)promedio,((alpha_max-alpha_min+1.0d0)**	&
	(1.0d0/real(N_alpha)))**i_alpha+alpha_min-1.0d0
  close(206)

 1206 format(E12.5,E16.8)
 1204 format(E12.5)
 1202 format(F5.3,F7.3)
 1200 format(F4.2,E12.5)
 1197 format(100E16.8)
 1000 format(1X,2E16.8)
end program histotef_ng19

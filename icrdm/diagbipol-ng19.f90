program diagbipol_ng19
! *******************************************************************************
! Nota para ingresar en forma correcta los valores de la matriz a diagonalizar	*
!										*
! Si KA es el numero de superdiagonales superiores, entonces KA+1 es el valor	*
! maximo del primer indice del vector AB, y corresponde a la diagonal. La	*
! diagonal se ingresa con el vector AB(KA+1,j) La siguiente superdiagonal se 	*
! construye como AB(KA,j). La última superdiagonal, sera ingresada mediante 	*
! el vector AB(1,j)								*
! ------------------------------------------------------------------------------*
! Nota para compilar bien el programa.						*
! 										*
! @ Hasta el momento solo es compilable teniendo las librerias:			*
! liblapack-doc, libblas, lapack3-dev, lapack3					*
!										*
! @ Al compilar, hay que linkear a lapack y blas, por ej mediante:		*
!    [directorio]$ g77 diagon.f -llapack -lblas -o diagon			* 
!-------------------------------------------------------------------------------*
!  Corregido el nivel de Fermi 13-12-2008					*
!-------------------------------------------------------------------------------*
!  Creación de Modulo comun para diagbipol e icrdm			*
!-------------------------------------------------------------------------------*
!  Version 14: corrige desviación en concentracion de quinoides			*
! sub-version a: altera la cantidad de electrones pra cambiar el nivel d fermi	*
! *******************************************************************************
  use diag_icrdm
  implicit none
  external DSBGV
  integer,parameter :: KA=2,N=(numanillos*7+1),KB=0,LDAB=N,LDBB=N,LDZ=N
  integer :: INFO,j,i,m,i2,h,unidad,cont,contbenceno,contquinoide,	 &
	seed(0:num_conf),registro(numanillos,num_conf),i_conf,		 &
	dist_qui,num_electrones,dif_electrones,cont_sort,cont_tol
  real(8) :: AB(LDAB,N),BB(LDBB,N),W(N),ran3,WORK(3*N),Z(LDZ,N), 	 &
	Fermi(num_conf),porcentaje,nivelfermi
  real :: max_qui,min_qui
  character (len=20) :: autov,configuraciones
  character :: JOBZ,UPLO,diagonq


  seed(0)=173811		! ej: 113317131
  dist_qui=1 	!Distancia mínima entre quinoides, en unidades de anillos bencénicos
  dif_electrones=0  !Agrega o quita electrones
  diagonq='N'   !Diagonaliza los hamiltonianos? N=no, S=si
  nivelfermi=-0.48d0  !valido solamente si se evitan las diagonalizaciones

  cont_tol=100000   !Toleranciade una semilla para encontrar configuraciones adecuadas
  min_qui=0.24   ! Fraccion de quinoides minima aceptable
  max_qui=0.26   !Fraccion de quinoides maxima aceptable

  configuraciones='config-ng19.inp'
  JOBZ='N'	!Opciones:  V => con autovectores, N => Sin autovectores
  UPLO='U'	!Como UPLO

! **** Informe de Variables *****
   print*
   print *, '************ Set de Variables ***************'
   print *, 'Seed= ', seed(0)
   print *, 'Cantidad de anillos=',numanillos
   print *, 'Numero de configuraciones=',num_conf
   print *, 'Configuracion inicial=', inc+1
   print *, 'Distancia entre quinoides=',dist_qui,' bencenos'
   print *
   if (diagonq .eq. 'N') then
      print *, 'No se incluye calculo exacto del nivel de Fermi'
      print *, 'Nivel de Fermi elegido:', nivelfermi
   elseif (diagonq .eq. 'S') then
      print *, 'Niveles de Fermi calculados exacatamente'
      print *, 'Alteracion en el # de electrones =', dif_electrones
   endif
   print *, '*********************************************'
   print*


  do i_conf=1,num_conf

! *****  Opciones y otras  *******************************
  autov='autobip15-'// char(48+int((i_conf+inc)/1000))//		 &
        char(48+int((i_conf+inc)/100)-10*int((i_conf+inc)/1000))//	 &
	char(48+int((i_conf+inc)/10)-10*int(int((i_conf+inc)/10)/10))//	 &
	char(48+i_conf+inc-10*int((i_conf+inc)/10))//'.dat'


! *** Sorteo del orden de los anillos *******************
  cont_sort=0
 8900 continue !Vuelve hasta aca si el sorteo no fue satisfactorio
  if (cont_sort .gt. cont_tol) then
	print *, 'La semilla actual impide encontrar configuraciones'
	print *, 'Con los margenes indicados de Quinoides'
	print *, '--------------------------------------------------'
	print *, 'Modificar alguno de estos parametros:'
	print *, 'seed(0)=', seed(0)
	print *, 'max_qui=', max_qui
	print *, 'min_qui=', min_qui
	print *, 'cont_tol=', cont_tol
	print *
	stop 
  endif

  seed(i_conf)=111+int(999999998*ran3(seed(i_conf-1))) !enteros entre 111 y 999999997
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
	   unidad = 1 + int(3*ran3(seed(i_conf)))
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
  porcentaje=real(contquinoide)/real(cont)
!  print *, porcentaje
  if (porcentaje .gt. max_qui .or. porcentaje .lt. min_qui) then
	cont_sort=cont_sort+1
	goto 8900
  endif


 11445  continue
  if (diagonq .eq. 'N' ) then
     Fermi(i_conf)=nivelfermi
     goto 11225
  elseif (diagonq .eq. 'S') then
     goto 11335
  else
     print*, 'Obtener las energias de Fermi exactas?(S/N)'
     read*, diagonq
     goto 11445
  endif

 11335 continue !Hacer diagonalizaciones


! *****  Arrays de diagonales  ***************************

! Diagonal (energias de sitio)  ---------------->
	j=0
	do m=1,numanillos
	if (registro(m,i_conf) .eq. 0) then !benceno:

	   if (m .ge. 1 ) then
	   if (registro(m-1,i_conf) .eq. 1) then
		AB(3,1+j)=-4.5d0	!n+
	   else
		AB(3,1+j)=-4.2d0  	!N
	   endif
	   else
		AB(3,1+j)=-4.2d0  	!N		
	   endif
		AB(3,2+j)=1.1d0 	!pC
		AB(3,3+j)=0.0d0		!oC
		AB(3,4+j)=0.0d0		!oC
		AB(3,5+j)=0.0d0		!oC
		AB(3,6+j)=0.0d0		!oC
		AB(3,7+j)=1.1d0		!pC

	else !quinoide:

		AB(3,1+j)=-4.5d0  	!N
		AB(3,2+j)=1.1d0 	!pC
		AB(3,3+j)=0.0d0		!oC
		AB(3,4+j)=0.0d0		!oC
		AB(3,5+j)=0.0d0		!oC
		AB(3,6+j)=0.0d0		!oC
		AB(3,7+j)=1.1d0		!pC
	endif
	j=j+7
	enddo
	AB(3,N)=-4.2d0		!N (termina en un nitrogeno)


! Primer superdiagonal  ----------------->
	j=1
	do m=1,numanillos
	if (registro(m,i_conf) .eq. 0) then !benceno:

		AB(2,1+j)=-3.4d0  	!N-pC
		AB(2,2+j)=-3.8d0 	!pC-oC
		AB(2,3+j)=0.0d0		!
		AB(2,4+j)=0.0d0		!
		AB(2,5+j)=0.0d0		!
		AB(2,6+j)=-3.8d0	!pC-oC 
		AB(2,7+j)=-3.4d0	!pC-N+

	else !quinoide:

		AB(2,1+j)=-4.2d0  	!N+-pC
		AB(2,2+j)=-3.4d0 	!pC-oC
		AB(2,3+j)=0.0d0		!
		AB(2,4+j)=0.0d0		!
		AB(2,5+j)=0.0d0		!
		AB(2,6+j)=-3.4d0	!pC-oC
		AB(2,7+j)=-4.2d0	!pC-N+
	endif
	j=j+7
	enddo


! Segunda superdiagonal  ---------------->
	j=2
	do m=1,numanillos
	if (registro(m,i_conf) .eq. 0) then !benceno:

		AB(1,1+j)=0.0d0  	!
		AB(1,2+j)=-3.8d0 	!pC-oC
		AB(1,3+j)=-3.35d0	!oC-oC
		AB(1,4+j)=-3.35d0	!oC-oC
		AB(1,5+j)=-3.8d0	!pC-oC
		AB(1,6+j)=0.0d0		!
	           if ((6+j) .eq. N ) then	!!!!!!!!!!!!!
  		   go to 333 !se termia de cargar los datos !!!!
	           endif   		!!!!!!!!!!!!!!!!!!!!!
		AB(1,7+j)=0.0d0	!

	else !quinoide:

		AB(1,1+j)=0.0d0  	!
		AB(1,2+j)=-3.4d0 	!pC-oC
		AB(1,3+j)=-3.8d0	!oC-oC
		AB(1,4+j)=-3.8d0	!oC-oC
		AB(1,5+j)=-3.4d0	!pC-oC
		AB(1,6+j)=0.0d0		!
	           if ((6+j) .eq. N ) then	!!!!!!!!!!!!!
  		   go to 333 !se termia de cargar los datos !!!!
	           endif   		!!!!!!!!!!!!!!!!!!!!!
		AB(1,7+j)=0.0d0	!	
	endif	
	j=j+7
	enddo

 333	continue             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Matriz identidad   ---------------------------------------------->
	do j=1,N
	BB(1,j)=1.0d0 !Como es la identidad, solo tiene una diagonal
	enddo


! ***** Llamada a la subrutina principal  ******************
        call DSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB,		 &
           LDBB, W, Z, LDZ, WORK, INFO )


! ***** ESeccion de control de errores  *************************** 
	print *
	print *,'**********************************************'

	if (INFO .LT. 0) then
		print *, 'El argumento # ', -INFO, ' es erroneo'
	elseif (INFO .gt. 0) then
	   if (INFO .le. N) then
		print *, 'ATENCION'
		print *
		print *,'El algoritmo no pudo converger. Los elementos ',&
      'diagonales ', INFO, ' no pudieron converger' 
		print *
           else
		print *, 'ATENCION'
		print *
		print *, 'La factorizacion de B no pudo completarse ',	 &
      'probablemente no sea definida positiva'
		print *,'La variable INFO arrojo el valor: ', info
		print *
	    endif
!	elseif (INFO .eq. 0) then
!		print *, 'Diagonalizacion exitosa!'
!		print *, 'Autovalores en orden ascendente en: ',autov
      
	endif

	print *, '**********************************************'
	print *


! *** Calculo de la Energia de Fermi ***********************
  num_electrones=8*numanillos-2*contquinoide+2 +dif_electrones   !Sólo para bipolarones
  Fermi(i_conf)=w(num_electrones/2)  


! ***** Escritura de resultados autoval  *************************** 
	open (unit=10, file=autov)
	do i=1,N
	write (10,1010) w(i),16.0d0
	enddo
	close(10)

 1010   format(2D16.8)

 11225 continue !Se saltea la diagonalizacion

  enddo ! i_conf


!***** Conf.inp  *******************************************
  open(15, file=configuraciones)
      write (15,1015) (Fermi(i_conf), i_conf=1,num_conf)
      do j = 1,numanillos
        write (15,1016) (registro(j,i_conf), i_conf=1,num_conf)
      enddo
  close(15)

 1014 format(I5)
 1015 format(100E16.8)
 1016 format(100I4)

end program diagbipol_ng19

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

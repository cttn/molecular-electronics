module MOD_nmrRAN
    !*** Propósito: Modulo general en el que se encuentran todos los métodos y
    !variables globales necesitados por el programa principal (depende, a su
    !vez, de MOD_MT)
    use mtmod, only: grnd
    implicit none
    integer, save :: xrango = 20
    contains
        subroutine leerdatos(N,arch,x,y)
            !*** Propósito: Leer archivo de datos de acuerdo a los parametros
            !esperados
            implicit none
            integer, intent(in) :: N
            integer :: i
            real(kind(0d0)), intent(out) :: x(N), y(N)
            character(len=30), intent(in) :: arch
            open(50, file = arch, status='old')
              do i = 1,N; read(50,*) x(i), y(i); enddo
            close(50)
        end subroutine leerdatos
        subroutine ransel_ppm_from_amp(i,N,y)            
            !*** Propósito: selecciona un valor ppm al azar, con probabilidad
            !directamente proporcional a la amplitud de ese ppm.
            implicit none
            integer, intent(out) :: i
            integer, intent(in)  :: N
            real(kind(0d0)), intent(in) :: y(N)
            real(kind(0d0)) :: ran, intervalo, suma
           
            intervalo = 0.d0
            ran       = grnd()
            suma      = sum(dabs(y))
           
            busqueda: do i = 1, N
              intervalo = intervalo + dabs(y(i))/suma
              if (ran .lt. intervalo) exit busqueda
            enddo busqueda
        end subroutine ransel_ppm_from_amp
        subroutine ransel_compen_cuat(i,N,x,y)
            !*** Propósito: selecciona un valor ppm al azar, con probabilidad
            !directamente proporcional a la amplitud de ese ppm. Compensa a los
            !carbonos cuaternarios aumentando su amplitud artificialmente.
            implicit none
            integer, intent(out) :: i
            integer, intent(in) :: N
            real(kind(0d0)), intent(in) :: x(N), y(N)
            integer, parameter :: Ncuat = 3 !# de ppm's de cuaternarios
            real(kind(0d0)) :: cuat(Ncuat), ymod(N), factor
            integer :: icuat(Ncuat),j,k, nrango

            nrango = 10 !rango de ppm dentro del cual compenso a las amplitudes
            factor = 2.d0 !factor por el cual lo compenso

            !-- ppm de carbonos cuaternarios, ordenados en forma creciente
            cuat(1) = 141.6d0
            cuat(2) = 145.7d0
            cuat(3) = 157.3d0

            !-- posicion de los ppm cuaternarios y compensacion
            do k=1,Ncuat
              icuat(k) = minloc(dabs(x(:)-cuat(k)),1)
              do j= icuat(k)-nrango/2, icuat(k)+nrango/2
                ymod(icuat(k)) = factor*ymod(icuat(k))
              enddo
            enddo

            !-- Calculo de probabilidad con amplitues compensadas
            call ransel_ppm_from_amp(i,N,ymod) 
        end subroutine ransel_compen_cuat
        subroutine redist_amp(i,N,y)
            !*** Propósito: Dado el punto i, redistribuir algo de la amplitud
            !asociada al ppm(i) en los ppm's cercanos (dentro del rango "xrango")
            implicit none
            integer, intent(in) :: i, N
            real(kind(0d0)), intent(inout) :: y(N)
            integer :: xsel
            real(kind(0d0)) :: yfrac, extracto

            yfrac  = 1.d-3

            ! -- Disminuir amplitud del punto indicado
            extracto  = y(i)*yfrac*grnd()
            y(i)      = y(i) - extracto

            ! -- Distribuir l amplitud
            xsel = i + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
            do while ( (xsel .lt. 0) .or. (xsel .gt. N) )
              xsel = i + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
            enddo
            y(xsel) = y(xsel) +  extracto
        end subroutine redist_amp
        subroutine escribir_res(uni,iter,N,x,y)
            !*** Propósito: Escribir resultados en un formato leíble por gnuplot
            !para mapas de colores (pm3d map)
            implicit none
            integer, intent(in) :: uni, iter, N
            real(kind(0d0)), intent(in) :: x(N), y(N)
            integer :: i
            do i = 1, N
                write(uni, 1000) real(iter,kind(0d0)), x(i), y(i)
            enddo
            write(uni,*)
            1000 format(3E16.8)
        end subroutine escribir_res
        subroutine nombrar(i, archivo_inicial,archivo)
            !*** Propósito: Generar un nombre para el archivo de resultados que
            !tome como base el nombre del archivo de datos inicial y algun
            !parámetro de interés (como "xrango")
            implicit none
            integer, intent(in) :: i
            character,intent(in) :: archivo_inicial*30
            character,intent(out) :: archivo*60
            write(archivo, *) i 
            archivo  = trim(adjustl(archivo_inicial)) // '_' // trim(adjustl(archivo)) // '.dat'
        end subroutine nombrar
        subroutine info_avance(i,maxi)
            !*** Propósito: Escribir en standar output informacion del avance de
            !las iteraciones
            implicit none
            integer, intent(in) :: i, maxi
            real :: veces, avance
            real, save :: curren = 1
            veces = 5.0
            avance = real(i)*100.0/real(maxi)
            if ( ( int(avance) .eq. int(100.0*curren/veces)) ) then
                write(*,*) int(avance), ' %'
                curren = curren + 1.0
            endif
        end subroutine info_avance         
end module MOD_nmrRAN

program nmrRAN
    !*** Propósito: Simular desorden en un espectro de 13C CPMAS NMR, en forma
    !progresiva
    use MOD_nmrRAN
    implicit none
    integer, parameter :: Npuntos = 4096
    integer :: isel, int_max, int_min, int_step, iter, Niter
    real(kind(0d0)) :: ppm(Npuntos), amp(Npuntos) !!Eje x y eje y del espectro
    character :: archivo_datos*30, archivo_res*60

    ! -- Parametros iniciales
    int_max  = 500000
    int_min  = 5000
    int_step = 5000
    archivo_datos = 'PAniES09.dat'

    call leerdatos(Npuntos,archivo_datos,ppm,amp)
    call nombrar(xrango,archivo_datos,archivo_res)

    open(70, file=archivo_res)

    do Niter = int_min, int_max, int_step
      do iter = 1, Niter
        !call ranselec_ppm_from_amp(isel,Npuntos,amp)
        call ransel_compen_cuat(isel,Npuntos,ppm,amp)
        call redist_amp(isel,Npuntos,amp)
      enddo
      call escribir_res(70,Niter,Npuntos,ppm,amp)
      call info_avance(Niter,int_max)
    enddo

    close(70)
end program nmrRAN

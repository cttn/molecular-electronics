program nmrRAN
    use mtmod, only: grnd
    implicit none
    integer, parameter :: Npuntos = 4096
    integer :: i, iter, xsel, thirdmax, thirdmin, thirdstep
    integer, target :: Niter, xrango
    integer, pointer :: thirdaxis
    real(kind(0d0)) :: xdat(Npuntos), ydat(Npuntos), suma, ran, intervalo,  &
        yfrac, extracto
    character :: inputfile*30, outputfile*60

    ! Parametros basicos
    xrango      = 100
    yfrac       = 1.d-3
    Niter       = 500000
    inputfile   = 'PAniES09.dat'

    ! Parametros 3er eje
    thirdaxis => Niter
    thirdmax  = 100000
    thirdmin  = 5000
    thirdstep = 5000

    open(50, file = inputfile, status='old')
    do i = 1 , Npuntos; read(50,*) xdat(i), ydat(i); enddo
    close(50); suma = sum(dabs(ydat))

    write(outputfile, *) xrango              !! Escribiendo resultados
    outputfile  = trim(adjustl(inputfile)) // '_' // trim(adjustl(outputfile)) //  '.dat'
    open(80, file=outputfile)

    extraparam: do thirdaxis = thirdmin, thirdmax, thirdstep
      do iter = 1, Niter
        ! Seleccionar un punto con probabilidad proporcional a la amplitud del punto
        intervalo = 0.d0
        i = 0
        ran = grnd()
        busqueda: do while( ran .gt. intervalo )
          i = i + 1
          intervalo = intervalo + dabs(ydat(i))/suma
          if (i .eq. Npuntos - 1) exit busqueda
        enddo busqueda

        ! Al punto seleccionado, bajarle la amplitud con probabilidad proporcional a
        ! la misma
        extracto  = ydat(i+1)*yfrac*grnd()
        ydat(i+1) = ydat(i+1) - extracto

        ! Sumarle esa amplitud a un punto x seleccionado al azar dentro de un rango
        ! xrango de sitios
        xsel = i + 1 + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
        do while ( (xsel .lt. 0) .or. (xsel .gt. Npuntos) )
          xsel = i + 1 + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
        enddo
        ydat(xsel) = ydat(xsel) + extracto
      enddo 

      do i = 1, Npuntos
        write(80, 1080) real(thirdaxis,kind(0d0)), xdat(i), ydat(i)
      enddo
      write(80,*)
    enddo extraparam

    1080 format(3E16.8)
    close(80)

end program nmrRAN

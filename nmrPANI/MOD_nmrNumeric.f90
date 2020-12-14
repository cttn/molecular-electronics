module mod_nmrNumeric
    use mtmod, only: grnd
    implicit none
    real(kind(0d0)), parameter :: ppm_max = 447.7517089843750000000d0 , & !ppm final
        ppm_min = -214.784764661552458165d0     !ppm inicial
    contains
    subroutine initdat(N,dat)
      !*** Propósito: Inicializar los "datos" es decir, las frecuencias a la que
      !responden los diferentes 13C de la muestra.
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(out) :: dat(N)
      dat(1) = 145.7d0
      dat(2) = 123.1d0 +0.5d0 !+2.d0
      dat(3) = 123.1d0 +0.5d0 !+2.d0
      dat(4) = 141.6d0
      dat(5) = 141.6d0
      dat(6) = 115.0d0 +4.5d0
      dat(7) = 157.3d0
      dat(8) = 135.7d0                           
    end subroutine initdat
    
    subroutine gen_ppm(ppm_max,ppm_min,N,ppm)
      !*** Propósito: Generar ppms similares a los del experimento
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(in) :: ppm_min, ppm_max
      real(kind(0d0)),intent(out) :: ppm(N)
      integer :: i
      do i=1,N
        ppm(i) = ppm_min+(ppm_max-ppm_min)*real(i-1,kind(0d0))/real(N-1,kind(0d0))
      enddo
    end subroutine gen_ppm

    subroutine gen_factorQ(N,factorQ)
      !*** Propósito: Generar un factor para las "poblaciones" de 13C 
      !dependiente de la posición. En principio, para distinguir entre 13C
      !cuaternarios y no cuaternarios.
      !Controla e tiempo de contacto del CP en la simulación.
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(out) :: factorQ(N)
      real(kind(0d0)) :: val
      val = 0.4d0
      factorQ(:) = 1.0d0
      factorQ(1) = val
      factorQ(4) = val  
      factorQ(5) = val
      factorQ(7) = val
    end subroutine gen_factorQ

    subroutine gen_factorQ_1param(N,factorQ,val)
      !*** Propósito: Generar un factor para las "poblaciones" de 13C 
      !dependiente de la posición. En principio, para distinguir entre 13C
      !cuaternarios y no cuaternarios.
      !Controla e tiempo de contacto del CP en la simulación.
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(out) :: factorQ(N)
      real(kind(0d0)),intent(in) :: val
      !val = 0.4d0
      factorQ(:) = 1.0d0
      factorQ(1) = val
      factorQ(4) = val  
      factorQ(5) = val
      factorQ(7) = val
    end subroutine gen_factorQ_1param

    subroutine gen_cantD_PAni(N,cantD)
      !*** Propósito: inicializa un vector cantD que contiene la cantidad de
      !veces que aparece cada 13C indicado en data
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)) :: cantD(N)
      cantD = 4.d0
      cantD(1) = 2.d0
      cantD(4) = 2.d0
      cantD(5) = 2.d0
      cantD(7) = 2.d0
      cantD(:) = cantD(:)*1.1d-2 !! Probabilidad de encontrar un 13C
    end subroutine gen_cantD_PAni

    subroutine read_cantD_PAni(infile,N,D)
      !*** Propósito: Lee de un archivo la cantidad de cada 13C
      implicit none
      integer,intent(in) :: N
      character,intent(in) :: infile*65
      real(kind(0d0)),intent(out) :: D(N)
      integer :: i
      open(45,file=infile,status='old')
      do i=1,N
        read(45,*)D(i)
      enddo
      close(45)
    end subroutine read_cantD_Pani

    subroutine gen_cpmas_spectra(eta,Ndat,dat,factorQ,cantD,Npts,ppm,amp)
      !*** Propósito: Genera un espectro de RMN mediante una suma de
      !lorentzianas compensadas. Cada lorenziana está centrada en la frecuencia
      !de cada uno de los 13C del sistema (datos).
      implicit none
      integer,intent(in) :: Ndat,Npts
      real(kind(0d0)),intent(in) :: eta, dat(Ndat), factorQ(Ndat), ppm(Npts), cantD(Ndat)
      real(kind(0d0)),intent(out) :: amp(Npts)
      integer :: i,j
      amp = 0.d0
      do j=1,Ndat
        do i=1,Npts
          amp(i) = amp(i) + lorentz(ppm(i),dat(j),eta)*factorQ(j)*cantD(j)
        enddo
      enddo
    end subroutine gen_cpmas_spectra

    subroutine read_frecPar(infile,N,P)
      !*** Propósito: Leer de un archivo las probabilidades de cruces
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(out) :: P(N,N)
      character,intent(in) :: infile*65
      integer :: i,j
      open(44,file=infile,status='old')
      do i=1,N
        read(44,*) (P(i,j),j=1,N)
      enddo
      close(44)
    end subroutine read_frecPar

    subroutine gen_collapseNMR_PAni(eta,rate,Ndat,dat,factorQ,Npts,ppm,amp,Col_amp,P)
      !*** Propósito: Colapsa las lineas de PAni buscando pares de 13C por
      !anillo
      implicit none
      integer,intent(in) :: Ndat, Npts
      integer :: i,j,k
      real(kind(0d0)),intent(in) :: dat(Ndat), ppm(Npts), amp(Npts),    &
          eta, factorQ(Ndat), rate, P(Ndat,Ndat)
      real(kind(0d0)) :: col_amp(Npts)
      col_amp(:) = amp(:)
      !-- Restando componentes originales del espectro
      do k=1,Ndat
        do j=1,Ndat
          do i=1,Npts
            col_amp(i)=col_amp(i)-lorentz(ppm(i),dat(j),eta)*P(k,j)*factorQ(j)
          enddo
        enddo
      enddo
      !-- Agregando colapso de las lineas restadas previamente
      do k=1,Ndat
        do j=1,Ndat
          do i=1,Npts
           col_amp(i) = col_amp(i) + &
            collapse_NMRlines_Qcomp(ppm(i),dat(k),dat(j),rate,eta,factorQ(k),factorQ(j))*P(k,j)
          enddo
        enddo
      enddo
    end subroutine gen_collapseNMR_PAni

    subroutine result_write1D(resfile,N,x,y)
      !*** Propósito: Escribe resultados cuando lo que se busca es un resultado
      !1D, es decir ppm vs amp.
      implicit none
      integer,intent(in) :: N
      real(kind(0d0)),intent(in) :: x(N), y(N)
      character(len=65) :: resfile
      integer :: i
      open(55, file=resfile)
        do i=1,N
          write(55,*) x(i), y(i)
        enddo
      close(55)
    end subroutine result_write1D

    function lorentz(epsil,centro,etaf)
      !*** Propósito: Construye la distribución de Cauchy o Lorentziana
      implicit none
      real(kind(0d0)), parameter :: Pi=dacos(-1.0d0)
      real(kind(0d0)) :: lorentz,centro,epsil,etaf
        lorentz=(1.0d0/pi)*(etaf/((epsil-centro)**2+etaf**2))
    end function      

    function cLorentz(epsil,centro,etaf)
      !*** Propósito: Construye una función imaginaria cuyo valor absoluto es la
      ! distribución de Cauchy o Lorentziana
      implicit none
      real(kind(0d0)), parameter :: Pi=dacos(-1.d0)
      real(kind(0d0)), intent(in) :: etaf, epsil
      complex(kind(0d0)),intent(in) :: centro
      complex(kind(0d0)) :: cLorentz, ci
      ci = dcmplx(0.d0,1.d0)
      cLorentz = (1.d0/pi)*(1.d0/(etaf + ci*(epsil-centro)))
    end function cLorentz

    function collapse_NMRlines(epsil,w1,w2,rate,etaf)
      !*** Propósito: Colapsa dos líneas de RMN a partir de sus frecuencias
      !iniciales, ensanchamiento inicial y tasa del intercambio 
      implicit none
      real(kind(0d0)),intent(in) :: w1, w2, epsil, rate, etaf
      real(kind(0d0)) :: collapse_NMRlines, wdif, wsum
      complex(kind(0d0)) :: R, cSignal, ci
      ci = dcmplx(0.d0,1.d0)
      wdif = w1 - w2
      wsum = (w1 + w2)*0.5d0
      R = cdsqrt(dcmplx((Wdif*0.5d0)**2,0d0)-dcmplx(rate**2,0d0)) 
      cSignal = 0.5d0*(1.d0 - ci*rate/R)*cLorentz(epsil,wsum + R, etaf + rate)      &
          + 0.5d0*(1.d0 + ci*rate/R)*cLorentz(epsil,wsum - R, etaf + rate)
      collapse_NMRlines = real(cSignal,kind(0d0))
    end function collapse_NMRlines

    function collapse_NMRlines_Qcomp(epsil,w1,w2,rate,etaf,factor1,factor2)
      !*** Propósito: Colapsa dos líneas de RMN a partir de sus frecuencias
      !iniciales, ensanchamiento inicial y tasa del intercambio 
      implicit none
      real(kind(0d0)),intent(in) :: w1, w2, epsil, rate, etaf, factor1, factor2
      real(kind(0d0)) :: collapse_NMRlines_Qcomp, wdif, wsum
      complex(kind(0d0)) :: R, cSignal, ci
      ci = dcmplx(0.d0,1.d0)
      wdif = w1 - w2
      wsum = (w1 + w2)*0.5d0
      R = cdsqrt(dcmplx((Wdif*0.5d0)**2,0d0)-dcmplx(rate**2,0d0)) 
      cSignal = 0.5d0*(1.d0 - ci*rate/R)*cLorentz(epsil,wsum + R, etaf + rate)*factor1      &
          + 0.5d0*(1.d0 + ci*rate/R)*cLorentz(epsil,wsum - R, etaf + rate)*factor2
      collapse_NMRlines_Qcomp = real(cSignal,kind(0d0))
    end function collapse_NMRlines_Qcomp

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
    
    subroutine redist_amp(i,xrango,N,y)
      !*** Propósito: Dado el punto i, redistribuir algo de la amplitud
      !asociada al ppm(i) en los ppm's cercanos (dentro del rango "xrango")
      implicit none
      integer, intent(in) :: i, N, xrango
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

    subroutine redist_amp_selec_param(i,xrango,fac,xsel,N,y)
      implicit none
      integer, intent(in) :: i, N, xrango
      real(kind(0d0)),intent(out) :: fac
      real(kind(0d0)), intent(inout) :: y(N)
      integer,intent(out) :: xsel
      real(kind(0d0)) :: yfrac, extracto
      yfrac = 1.d-3
      fac   = yfrac*grnd()
      ! -- Disminuir amplitud del punto indicado
      extracto  = y(i)*fac
      y(i)      = y(i) - extracto
      ! -- Distribuir l amplitud
      xsel = i + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
      do while ( (xsel .lt. 1) .or. (xsel .gt. N) )
        xsel = i + int(real(xrango,kind(0d0))*(grnd() - 0.5d0))
      enddo
      y(xsel) = y(xsel) +  extracto
    end subroutine redist_amp_selec_param

    subroutine redist_amp_take_param(i,xsel,fac,N,y)
      implicit none
      integer, intent(in) :: i, N, xsel
      real(kind(0d0)),intent(in) :: fac
      real(kind(0d0)), intent(inout) :: y(N)
      real(kind(0d0)) :: extracto
      extracto = y(i)*fac
      y(i)     = y(i) - extracto
      y(xsel)  = y(xsel) + extracto
    end subroutine redist_amp_take_param

    subroutine nombrar(i, archivo, prefix, sufix)
      !*** Propósito: Generar un nombre para el archivo de resultados que
      !tome como base el nombre del archivo de datos inicial y algun
      !parámetro de interés (como "xrango")
      implicit none
      integer, intent(in)   :: i
      character,intent(in)  :: prefix*65
      character,intent(out) :: archivo*65
      character,optional,intent(in) :: sufix*65
      character :: default_sufix*65 = '.dat'
      write(archivo, *) i
      if (present(sufix)) default_sufix = sufix
      archivo = trim(adjustl(prefix)) // trim(adjustl(archivo)) // default_sufix
      archivo = trim(adjustl(archivo))
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
      if (i .eq. 1) print *, 'Avance:'
      if ( ( int(avance) .eq. int(100.0*curren/veces)) ) then
         write(*,*) int(avance), ' %'
         curren = curren + 1.0
      endif
    end subroutine info_avance  
end module mod_nmrNumeric



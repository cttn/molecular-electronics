program numNMRPAni
  use modNMRPAni
  use modstatPAni
  use modIOUtil
  implicit none
  integer, parameter :: Nanil = 400, NtipAnil = 5, NC_anil = 6, N13C = 8 !! stat PAni
  integer, parameter :: Npts = 4096  !!para linea RMN
  integer, parameter :: Niter = 100 !!Numero de iteraciones a promediar
  integer :: cadena(Nanil), anillo(NtipAnil,NC_anil), Rango, i, j, iter
  real(kind(0d0)) :: fracInd(N13C), fracPar(N13C,N13C),prat,pran
  real(kind(0d0)) :: ppm(Npts), amp(Npts), col_amp(Npts), chemShift(N13C), &
    eta, factor(N13C), tc, rate, AV_amp(Npts)=0.d0

  eta = 8.d0                                !Ensanchamiento original
  tc  = 0.8d0                               !Parámetro tiempo contacto
  pRat = 0.35d0
  pRan = 2.50d0

!------------
  rate  = pRat*eta
  Rango = nint(pRan*eta)

  !-- Linea RMN
  call initdat(N13C,chemShift)
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)
  call gen_factorQ_1param(N13C,factor,tc)

  call genPAniAnil(NtipAnil,NC_anil,anillo)

  do iter=1,Niter
    !-- Calculando parámetros estáticos de PAni
    call sortQui(Nanil,cadena)
    call genCad(Nanil,cadena)
    call getFracInd(NtipAnil,NC_anil,anillo,Nanil,cadena,N13C,fracInd)

    !-- Info rapida
    call printFracAnil(NtipAnil,Nanil,cadena)

    !-- Linea RMN sin colapsos
    call cpmas_spectra(eta,N13C,chemShift,factor,fracInd,Npts,ppm,amp)

    !-- Linea RMN con colapsos
    call getFracPar(NtipAnil,NC_anil,anillo,Nanil,cadena,N13C,fracPar,Rango)
    call collapseNMR_PAni(eta,rate,N13C,chemShift,factor,Npts,ppm,amp,col_amp,fracPar)  

    AV_amp(:) = AV_amp(:)+col_amp(:)
  enddo; AV_amp(:)=AV_amp(:)/real(Niter,kind(0d0))


  do i=1,Npts
    write(77,*)ppm(i),amp(i),AV_amp(i)
  enddo

  !-- Info rapida
  !call printFracAnil(NtipAnil,Nanil,cadena)

end program numNMRPAni

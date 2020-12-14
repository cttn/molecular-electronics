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
  real(kind(0d0)) :: ppm(Npts), amp(Npts), col_amp(Npts), chemShift(N13C),  &
    eta, factor(N13C), tc, rate, AV_amp(Npts)=0.d0, AV_col_amp(Npts)=0.d0,  &
    normalize_AV, width(N13C),etaQ, etaNQ, normalize_colAv,                 &
    fracAnil(NtipAnil),AV_fracAnil(Ntipanil)
  logical :: collapse
  character :: namefile*20

  AV_fracAnil(:) = 0.0d0 

  !-- Parámetros fijos
  collapse = .true.
  namefile = 'etaRange.dat'
  etaQ  = 3.0d0  !+12.d0   !2.8d0
  etaNQ = 3.7d0  !+12.d0   !3.83d0
  tc    = 0.7d0                                  !Parámetro tiempo contacto

  !-- Parámetros de colapso
  eta = 4.50d0                                       !Ensanchamiento extra
  pRat = 0.35d0
  pRan = 2.5d0

  !-- Rate and Range
  rate  = pRat*eta ; Rango = int(pRan*eta)

  !-- Linea RMN
  if (.not. collapse) eta = 0
  call initdat(N13C,chemShift)
  call initWidth(N13C,width,eta,etaQ,etaNQ)
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)
  call gen_factorQ_1param(N13C,factor,tc)

  call genPAniAnil(NtipAnil,NC_anil,anillo)

  do iter=1,Niter
    !-- Calculando parámetros estáticos de PAni
    call sortQui(Nanil,cadena)
    call genCad(Nanil,cadena)
    call getFracInd(NtipAnil,NC_anil,anillo,Nanil,cadena,N13C,fracInd)

    !-- Estadistica de Quinoides
    call getFracAnil(NtipAnil,Nanil,cadena,fracAnil)
    AV_fracAnil(:) = Av_fracAnil(:) + fracAnil(:)

    !-- Linea RMN sin colapsos
    call cpmas_spectra(N13C,chemShift,width,factor,fracInd,Npts,ppm,amp)
    AV_amp(:) = AV_amp(:) + amp(:)

    if(collapse) then
      !-- Linea RMN con colapsos
      call getFracPar(NtipAnil,NC_anil,anillo,Nanil,cadena,N13C,fracPar,Rango)
      call collapseNMR_PAni(rate,N13C,chemShift,width,factor,Npts,ppm,amp,col_amp,fracPar)  
      AV_col_amp(:) = AV_col_amp(:)+col_amp(:)
    endif

    call info_avance(iter,Niter)
  enddo 
  
  AV_amp(:) = Av_amp(:)/real(Niter,kind(0d0))
  normalize_AV   = maxval(AV_amp)

  open(77,file=namefile)
  if (collapse) then
    AV_col_amp(:)=AV_col_amp(:)/real(Niter,kind(0d0))
    normalize_colAV = maxval(AV_col_amp)
    do i=1,Npts
      write(77,'(10E16.8)')ppm(i),AV_amp(i),AV_col_amp(i),AV_amp(i)/normalize_AV,AV_col_amp(i)/normalize_colAV
    enddo
  else
    do i=1,Npts
      write(77,'(10E16.8)')ppm(i),AV_amp(i),AV_amp(i)/normalize_AV
    enddo
  endif
  close(77)


  print *, 'AV % Quinoides: ', AV_fracAnil(3)*100.0d0/real(Niter)
  print *, 'Archivo: ', namefile

  !-- Info rapida
  !print *, 'Eta: ', eta
  !print *, 'Rate: ', rate
  !print *, 'Rango: ', rango
  !print *
  !call printFracAnil(NtipAnil,Nanil,cadena)

end program numNMRPAni

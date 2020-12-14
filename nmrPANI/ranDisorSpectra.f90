program nmrRanDisor
  use mod_nmrNumeric
  implicit none
  integer, parameter :: Npts = 4096, Ndat = 8
  integer :: isel, iter, Niter, xrango, xsel
  real(kind(0d0)) :: ppm(Npts),amp(Npts),dat(Ndat),eta,factor(Ndat),fac,fakeAmp(Npts),cantD(Ndat)
  character(len=65) :: outfile, prefix
  logical :: compensar_13C_cuaternarios

  Niter  = 5*10**5  ; print *, '# iteraciones: ', Niter
  xrango = 100      ; print *, 'Rango ppm: ', xrango
  eta    = 4.d0     ; print *, 'Ensanchamiento previo: ', eta
  compensar_13C_cuaternarios = .true.
  if (compensar_13C_cuaternarios) then
    prefix = 'c_DisSpec'
  else  
    prefix = 'nc_DisSpec'
  endif
  
  call nombrar(xrango,outfile,prefix)                                 !Generando nombre de archivo para resultados
  call initdat(Ndat,dat)                                              !Inicializando datos
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)                              !Generando valores ppm
  call gen_factorQ(Ndat,factor)                                        !Generando factor de compensación para 13C cuaternarios
  call gen_cantD_PAni(Ndat,cantD)
  call gen_cpmas_spectra(eta,Ndat,dat,factor,cantD,Npts,ppm,amp)            !Suma de Lorentzianas para generar el espectro de RMN.

  if (compensar_13C_cuaternarios) then
    factor(:) = 1.d0       !-- Generando espectro sin compensar
    call gen_cpmas_spectra(eta,Ndat,dat,factor,cantd,Npts,ppm,fakeAmp)   
    do iter = 1, Niter
      call ransel_ppm_from_amp(isel,Npts,fakeamp)                     !Distribución sin compensancion a los 13C no cuaternarios
      call redist_amp_selec_param(isel,xrango,fac,xsel,Npts,fakeamp)  !Re-distribución sin compensación
      call redist_amp_take_param(isel,xsel,fac,Npts,amp)              !Redistribución del espectro real
      call info_avance(iter,Niter)
    enddo
  else
    print *,'Aviso: No se compensa la probabilidad para 13C Cuaternarios'
    do iter = 1, Niter
      call ransel_ppm_from_amp(isel,Npts,amp)                         !Distribución sin compensancion a los 13C no cuaternarios
      call redist_amp(isel,xrango,Npts,amp)                           !Redistribución del espectro real
      call info_avance(iter,Niter)
    enddo
  endif

  call result_write1D(outfile,Npts,ppm,amp)

end program nmrRanDisor



program paniCollapse
  use mod_nmrNumeric
  implicit none
  integer, parameter :: Npts = 4096, Ndat = 8
  integer :: i,j
  real(kind(0d0)) :: ppm(Npts), amp(Npts), dat(Ndat), eta, factor(Ndat),    &
      val, cantD(Ndat), col_amp(Npts), rate, P(Ndat,Ndat)
  character(len=65) :: outfile,prefix,sufix,cantdfile,frecParFile

  prefix = 'xPAniLorentz_tc80%_Q'              !Archivo de resultados
  sufix  = '_R10.dat'
  cantdfile   = 'cantDPAni_N400.dat'
  frecParFile = 'frecParPAni_N400_R10.dat'
  eta    = 4.0d0                               !Ensanchamiento a usar
  rate   = 8.d0
  val    = 0.8d0

  call initdat(Ndat,dat)                       !Inicializando frecuencias
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)       !Generando valores ppm
  call gen_factorQ_1param(Ndat,factor,val)     !Compensación para 13C cuatern
 call gen_cantD_PAni(Ndat,cantD)              !Compensación por cant de 13C
!  call read_cantD_Pani(cantDfile,Ndat,cantD)

  !Suma de Lorenzianas para generar el espectro de RMN.
  call gen_cpmas_spectra(eta,Ndat,dat,factor,cantD,Npts,ppm,amp)

  call read_frecPar(frecParFile,Ndat,P)        !Leyendo probab. pares 13C
  call gen_collapseNMR_PAni(eta,rate,Ndat,dat,factor,Npts,ppm,amp,col_amp,P)

! call nombrar(int(val*100.),outfile,prefix,sufix) !Generando nombre de archivo
  call nombrar(int(eta),outfile,prefix,sufix)      !Generando nombre de archivo
  call result_write1D(outfile,Npts,ppm,col_amp)    !Guardando resultados.


end program paniCollapse

program pani_eta_map
  use mod_nmrNumeric
  implicit none
  integer, parameter :: Npts = 4096, Ndat = 8
  integer :: i,j
  real(kind(0d0)) :: ppm(Npts), amp(Npts), dat(Ndat), eta, 	&
    eta_min, eta_max, eta_pts, factor(Ndat)
  character(len=65) :: outfile

  outfile='map_cpmas.dat'
  eta_min = 4.d0
  eta_max = 8.d0
  eta_pts = 100
  
  call initdat(Ndat,dat)                               !Inicializando datos
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)               !Generando valores ppm
  call gen_factor(Ndat,factor)                         !Generando factor de compensación para 13C cuaternarios
  
  open(45,file=outfile)
    do i=1,100
      eta = eta_min+(eta_max-eta_min)*real(i-1)/real(eta_pts-1)
      call gen_cpmas_spectra(eta,Ndat,dat,factor,Npts,ppm,amp)
      do j=1600,2600
        write(45,*) ppm(j), eta, amp(j)
      enddo
      write(45,*)
    enddo
  close(45)

  print *, 'MAPA: Resultados guardados en: ', outfile

end program pani_eta_map

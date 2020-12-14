program nmrRanDisor
  use mod_nmrNumeric
  implicit none
  integer, parameter :: Npts = 4096, Ndat = 8
  integer :: isel, int_max, int_min, int_step, iter, Niter, xrango
  real(kind(0d0)) :: ppm(Npts), amp(Npts), dat(Ndat), eta, factor(Ndat), &
       fakeAmp(Npts)
  character(len=65) :: outfile, prefix

  prefix   = 'ranDisMap'
  int_max  = 140000
  int_min  = 140000
  int_step = 140000
  xrango   = 100
  eta 	   = 4.d0
  
  call nombrar(xrango,outfile,prefix)                   !Generando nombre de archivo para resultados

  call initdat(Ndat,dat)                                      !Inicializando datos
  call gen_ppm(ppm_max,ppm_min,Npts,ppm)                      !Generando valores ppm

  !-- Generando espectro sin compensar
  factor(:)  = 1.d0
  call gen_cpmas_spectra(eta,Ndat,dat,factor,Npts,ppm,fakeAmp) 
  
  !-- Generando espectro compensado
  call gen_factor(Ndat,factor)                                !Generando factor de compensación para 13C cuaternarios
  call gen_cpmas_spectra(eta,Ndat,dat,factor,Npts,ppm,amp)    !Suma de Lorenzianas para generar el espectro de RMN.

  open(70, file=outfile)

    do Niter = int_min, int_max, int_step
      do iter = 1, Niter
        call ransel_ppm_from_amp(isel,Npts,fakeAmp)
        call redist_amp(isel,xrango,Npts,amp)
      enddo
      call escribir_res(70,Niter,Npts,ppm,amp)
      call info_avance(Niter,int_max)
    enddo

  close(70)

end program nmrRanDisor

subroutine escribir_res(uni,iter,N,x,y)
   !*** Propósito: Escribir resultados en un formato leíble por gnuplot
   !para mapas de colores (pm3d map)
   implicit none
   integer, intent(in) :: uni, iter, N
   real(kind(0d0)), intent(in) :: x(N), y(N)
   integer :: i
   do i = 1600,2600 !N
     write(uni, 1000) real(iter,kind(0d0)), x(i), y(i)
   enddo
   write(uni,*)
   1000 format(3E16.8)
end subroutine escribir_res

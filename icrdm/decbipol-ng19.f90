!************************************************************************
! Decima anillos bencenoides. Aquí se usa para decimar bencenos y	*
! quinoides, que serán pasados a los programas de transporte incoherente	*
!************************************************************************
! NOTA1: El rango [Emin,Emax] debe coincidir con el ancho que se quiera
! ver en los programas de transporte incoherente, [-2V,2V] donde V
! es el parametro de hopping para los leads
module decbipol
  integer, parameter :: puntos_ener = 200
end module decbipol


program decbipol_ng19
  use decbipol
  implicit none
  integer :: s
  complex(8) :: enerorto,enerpara,vop_qui,vop_ben,voo_qui,voo_ben,	 &
	epsi(puntos_ener),hop_Npara,hop_Npara_q,ener_N,ener_Nmas,	 &
	Ep_ben(puntos_ener),Ep_qui(puntos_ener),vpp_ben(puntos_ener),	 &
	vpp_qui(puntos_ener)
  real(8) :: eta,Emax,Emin,pi
  character (len=13) archivo

! Variables manuales
  archivo='decbipol.inp'
  Emax=2.0d0
  Emin=-2.0d0
  eta=0.000000001d0
  pi=dacos(-1.0d0)

  hop_Npara=dcmplx(-3.4d0,0.0d0)
  hop_Npara_q=dcmplx(-4.2d0,0.0d0)
  ener_N=dcmplx(-4.2d0,0.0d0)
  ener_Nmas=dcmplx(-4.5d0,0.0d0)

  enerpara=dcmplx(1.1d0,0.0d0)	
  enerorto=dcmplx(0.0d0,0.0d0)
  vop_qui=dcmplx(-3.4d0,0.0d0)
  voo_qui=dcmplx(-3.8d0,0.0d0)
  vop_ben=dcmplx(-3.8d0,0.0d0)
  voo_ben=dcmplx(-3.35d0,0.0d0)


  open(unit=75,file=archivo)
	
  write (75,1077) puntos_ener
  write(75,1076) eta,real(hop_Npara),aimag(hop_Npara),		 	 &
  real(ener_N),aimag(ener_N),real(ener_Nmas),aimag(ener_Nmas),		 &
  real(hop_Npara_q),aimag(hop_Npara_q)


  do s=1,puntos_ener
  epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s)/real(puntos_ener),-eta)


! Benceno:
  Ep_ben(s)=enerpara+2.*cdabs(vop_ben)**2/(epsi(s)-enerorto-		 &
	 cdabs(voo_ben)**2/(epsi(s)-enerorto))

  Vpp_ben(s)=2.*voo_ben*cdabs(vop_ben)**2/((epsi(s)-enerorto)**2-		 &
	cdabs(voo_ben)**2)

! Quinoide:
  Ep_qui(s)=enerpara+2.*cdabs(vop_qui)**2/(epsi(s)-enerorto-		 &
	 cdabs(voo_qui)**2/(epsi(s)-enerorto))

  Vpp_qui(s)=2.*voo_qui*cdabs(vop_qui)**2/((epsi(s)-enerorto)**2-		 &
	cdabs(voo_qui)**2)


! Resultados:  
  write(75,1075) real(epsi(s)),real(Ep_qui(s)),aimag(Ep_qui(s)), 	 &
	real(vpp_qui(s)),aimag(vpp_qui(s)),real(Ep_ben(s)),		 &
	aimag(Ep_ben(s)),real(vpp_ben(s)),aimag(vpp_ben(s))


  enddo

  close(75)
 1075   format(15E16.8)
 1076   format(9E16.8)
 1077	format(1I5)

end program decbipol_ng19

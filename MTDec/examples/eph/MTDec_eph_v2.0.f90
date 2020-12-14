!########################################################################!
!#  FALTANTES A OPTIMIZAR y/o AGREGAR :					#!
!#     * Calculos cuando T(i,j) .ne. T(j,i)				#!
!#     * Calculos con matrices matv,math,matd a eliminar		#!
!#     * Hay un Lead por canal!						#!
!#     * Revertir Locacion de memoria, especialmente subrutinas.	#!
!#     * Diferenciar el caso matriz T simétrica o no y elegir subrut    #!
!#        eb funcion de eso.						#! 
!########################################################################!
module MultiDP_global
  implicit none
  integer,save :: BLOCK_DIMENSION, BLOCK_NUMBER, RESOLUTION, CHANNEL_NUMBER,		&
	DEC_NUMBER
  integer :: i,j,m,s,k,l,H,N,M_comp=0,M_real=0
  integer,allocatable,save :: chan_pos(:),dec_pos(:)
  complex(kind(0d0)),allocatable,save :: ener(:,:,:),hop(:,:,:),			&
	Sigma_D(:),epsi(:)						
  real(kind(0d0)),allocatable,save :: ener_lead(:),hop_lead(:),tun_lead(:)
  real(kind(0d0)),parameter :: pi=dacos(-1.d0)
  character,save :: decoherence,arch_trans*30
  contains
!%%% Sigma LEAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  function Sigma_L(var,e)
   complex(kind(0d0)) :: Sigma_L
   integer :: e,var
     if    ((real(epsi(e))-ener_lead(var)) .gt. 2.0d0*dabs(hop_lead(var))) then 
	Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0                   		&
        -dsqrt(((real(epsi(e))-ener_lead(var))/2.0d0)**2-hop_lead(var)**2),0.0d0)   
     elseif(dabs(real(epsi(e))-ener_lead(var)) .le. dabs(2.0d0*hop_lead(var))) then
	Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0,-dsqrt(hop_lead(var)**2 	&
        -((real(epsi(e))-ener_lead(var))/2.0d0)**2))
     elseif((real(epsi(e))-ener_lead(var)) .lt. -2.0d0*dabs(hop_lead(var))) then 
	Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0                 		&
        +dsqrt(((real(epsi(e))-ener_lead(var))/2.0d0)**2-hop_lead(var)**2),0.0d0)
     endif
!-- Tunneling hacia la muestra -----
	Sigma_L=Sigma_L*(tun_lead(var)/hop_lead(var))**2
   end function
!%%% COMPLEX SQUARE MATRIX INVERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine invert(N,A)
    external ZGETRF
    external ZGETRI
       integer :: INFO, LDA, LWORK, N
       integer :: IPIV(N)
       complex(kind(0d0)),intent(inout) :: A(N,N)
       complex(kind(0d0)),allocatable :: WORK(:),Asim(:,:)
      if (m_comp .eq. 0) then
        LWORK=-1
        LDA=N
        allocate(work(N))
        allocate(ASim(N,N))
      else
	LDA=N
	LWORK=m_comp
        allocate(work(m_comp))
      endif
 223 continue
       if((lwork .eq. -1) .and. (m_comp .eq. 0)) then
        CALL ZGETRF(N,N,Asim,LDA,IPIV,INFO)
        CALL ZGETRI(N,Asim,LDA,IPIV,WORK,LWORK,INFO)
        M_comp=work(1)
	deallocate(work)
	lwork=m_comp
        allocate(work(m_comp))
        CALL ZGETRF(N,N,A,LDA,IPIV,INFO)
        CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
	goto 223
       else
        CALL ZGETRF(N,N,A,LDA,IPIV,INFO)
        CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
       endif       
  end subroutine
!%%% REAL SQUARE MATRIX INVERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine invert_real(d,A)
   external DGETRF
   external DGETRI
     integer :: INFO, LWORK, d
     integer :: IPIV(d)
     real(kind(0d0)),intent(inout) :: A(d,d)
     real(kind(0d0)),allocatable :: WORK(:),Asim(:,:)
       select case (m_real)
	 case(0)
           LWORK=-1
           allocate(work(d))
  	   allocate(ASim(d,d))
           CALL DGETRF(d,d,Asim,d,IPIV,INFO)
           CALL DGETRI(d,Asim,d,IPIV,WORK,LWORK,INFO)
           M_real=work(1)
           deallocate(work)
	   lwork=m_real
           allocate(work(m_real))
           CALL DGETRF(d,d,A,d,IPIV,INFO)
           CALL DGETRI(d,A,d,IPIV,WORK,LWORK,INFO)
         case default
	   LWORK=m_real
           allocate(work(m_real))
           CALL DGETRF(d,d,A,d,IPIV,INFO)
           CALL DGETRI(d,A,d,IPIV,WORK,LWORK,INFO)
       end select
  end subroutine
end module MultiDP_global

program MultiDP
  use MultiDP_global
  implicit none
  character(len=60) :: arch_all,arch_tot
!--- Interfase ------------------------------------------------------------------!
  decoherence='Y' !Yes/NO
  arch_trans='mtd_e-ph'

  RESOLUTION=10000   	!Number of energy points
  allocate(epsi(RESOLUTION))


  BLOCK_DIMENSION=8 	!Dimension of blocks
  BLOCK_NUMBER=1  	!Number of blocks
  allocate(ener(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  if (block_number .gt. 1) then
   allocate(hop(BLOCK_NUMBER-1,BLOCK_DIMENSION,BLOCK_DIMENSION))
  endif


  CHANNEL_NUMBER=2*BLOCK_DIMENSION  	!Numbers of real channels
  allocate(chan_pos(CHANNEL_NUMBER))
  allocate(ener_lead(CHANNEL_NUMBER))
  allocate(hop_lead(CHANNEL_NUMBER))
  allocate(tun_lead(CHANNEL_NUMBER))
   do i=1,CHANNEL_NUMBER
    chan_pos(i)=1+int((i-1)/2)   !Impar: Left, Par: Right
   enddo


  if (decoherence .eq. 'Y') then
    DEC_NUMBER=BLOCK_NUMBER*BLOCK_DIMENSION   !Number of decoherent channels
    allocate(dec_pos(DEC_NUMBER))
    allocate(Sigma_D(DEC_NUMBER))
      do i=1,DEC_NUMBER
        dec_pos(i)=i
      enddo
  endif

     !-------- Extra: Nombres de Archivo automáticos
	  if (decoherence .eq. 'Y') then
	    arch_all=trim(adjustR(arch_trans)) // '_dec.dat'
	    arch_tot=trim(adjustR(arch_trans)) // '_tot_dec.dat'
 	  else
  	    arch_all=trim(adjustR(arch_trans)) // '_coh.dat'
   	    arch_tot=trim(adjustR(arch_trans)) // '_tot_coh.dat'
  	  endif

  open(77, file=arch_all)
  open(79, file=arch_tot)

  call asignations_eph
  call DAmatoPastawski

  close(77)
  close(79)
end program MultiDP


subroutine asignations_eph
  use MultiDP_global
  implicit none
  real(kind(0d0)) :: kT,Emin,Emax,eta,hw,vg,ener_0

  kT=0.025d0 !eV
  Emin= 0.d0
  Emax= 1.5d0 
  eta=1.d-8

  do s=1,RESOLUTION
   epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s,kind(0d0))/real(RESOLUTION,kind(0d0)),eta)
  enddo  

!*** Parametros e-ph
  hw=0.2d0
  Vg=0.1d0
  ener_0=0.5d0

!--- Leads asignations -------------------------------------
  !El canal 1-2 corresponde a cero fonones, 3-4 a un fonon, etc
  !Entro siempre por el lead 1.
  do i=1,CHANNEL_NUMBER
    ener_lead(i)=2.05d0- dble(int((i-1)/2))*hw     
    hop_lead(i)=1.d0
    tun_lead(i)=0.1d0
  enddo

!--- Decs asignations --------------------------------------
  if (decoherence .eq. 'Y') then
    do i=1,DEC_NUMBER
      Sigma_D(i)=dcmplx(0.d0,-kT)
    enddo
  endif

!--- Energias tight-binding
  do i=1,block_dimension
    do j=1,block_dimension
      if (i.eq.j) then
        ener(1,i,j)=dcmplx(ener_0+hw*(real(i-1,kind(0d0))+0.5d0)-(dabs(Vg)**2)/hw,0.d0)
      elseif ((i.eq. j+1) .or. (j.eq. i+1)) then
        ener(1,i,j)=dcmplx(dsqrt(real(i+1,kind(0d0)))*Vg,0.d0)
      else
        ener(1,i,j)=dcmplx(0.d0,0.d0)
      endif
    enddo
  enddo
end subroutine asignations_eph

subroutine DAmatoPastawski
  use MultiDP_global
  implicit none
  integer :: i_dim,j_dim
  real(kind(0d0)),allocatable :: ID(:,:),dos(:),T_LL(:,:),T_LD(:,:),T_DD(:,:),	&
	T_ef(:,:)
  complex(kind(0d0)),allocatable :: matv(:,:),math(:,:),matd(:,:),inver(:,:),	&
	Green(:,:,:,:)
    type selfenergies
	complex(kind(0d0)),allocatable :: L(:,:,:),R(:,:,:)
    end type selfenergies
  type (selfenergies) :: delta,sim_delta

  allocate(matv(BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(matd(BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(math(BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(inver(BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(Green(BLOCK_NUMBER,BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(T_LL(CHANNEL_NUMBER,CHANNEL_NUMBER))

  if (decoherence .eq. 'Y') then
    allocate(T_LD(CHANNEL_NUMBER,DEC_NUMBER))
    allocate(T_DD(DEC_NUMBER,DEC_NUMBER))
    allocate(T_ef(CHANNEL_NUMBER,CHANNEL_NUMBER))
  endif

!--- Identity Matrix ----------------------------------------
  allocate(ID(BLOCK_DIMENSION,BLOCK_DIMENSION))
  do i=1,BLOCK_DIMENSION
  do j=1,BLOCK_DIMENSION
	if (i .eq. j) then
	 ID(i,j)=1.d0
	else
	 ID(i,j)=0.d0
  	endif
  enddo !j
  enddo !i
 
!--- Adding self energies of environment --------------------
  if (decoherence .eq. 'Y' ) then
    do i=1,DEC_NUMBER
      m=1+mod(dec_pos(i)-1,BLOCK_DIMENSION)
      k=1+int((dec_pos(i)-1)/BLOCK_DIMENSION)
      ener(k,m,m)=ener(k,m,m)+ sigma_D(i)
    enddo !i
  endif

  allocate(delta%L(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(delta%R(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
!  allocate(sim_delta%L(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(sim_delta%R(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  allocate(dos(RESOLUTION))
  dos=0.d0

  do s=1,RESOLUTION

    do i=1,CHANNEL_NUMBER
      m=1+mod(chan_pos(i)-1,BLOCK_DIMENSION)
      k=1+int((chan_pos(i)-1)/BLOCK_DIMENSION)
      ener(k,m,m)=ener(k,m,m)+ sigma_L(i,s)
    enddo !i

    delta%L=dcmplx(0.d0,0.d0)
    delta%R=dcmplx(0.d0,0.d0)
!    sim_delta%L=dcmplx(0.d0,0.d0)
    sim_delta%R=dcmplx(0.d0,0.d0)

    !*** Delta L (menos) ************************************************
    do l=2,BLOCK_NUMBER

      math(:,:)=hop(l-1,:,:)		
      inver(:,:)=epsi(s)*ID(:,:)-ener(l-1,:,:)-delta%L(l-1,:,:)

      call invert(BLOCK_DIMENSION,inver)
      matd=matmul(matmul(transpose(math),inver),math)

      delta%L(l,:,:)=matd(:,:)
    enddo !l

    !******** Delta R (mas) **************************************
    do l=2,BLOCK_NUMBER 

      math(:,:)=hop(BLOCK_NUMBER+1-l,:,:) 
      inver(:,:)=epsi(s)*ID(:,:)-ener(BLOCK_NUMBER+2-l,:,:)	&
		-delta%R(BLOCK_NUMBER+2-l,:,:)
  
      call invert(BLOCK_DIMENSION,inver)
      matv=matmul(math,inver)
      matd=matmul(matv,transpose(math))

      delta%R(BLOCK_NUMBER+1-l,:,:)=matd(:,:)
      sim_delta%R(BLOCK_NUMBER+1-l,:,:)=matv(:,:)
    enddo !l

    !--- Green Function Matrix: Diagonal elements --------------------
    do i=1,BLOCK_NUMBER
      inver(:,:)=epsi(s)*ID(:,:)-ener(i,:,:)-delta%L(i,:,:)-delta%R(i,:,:)
      call invert(BLOCK_DIMENSION,inver)

      Green(i,i,:,:)=inver(:,:)
    enddo !i

    !******************** DoS **************************************!
    !do l=1,BLOCK_NUMBER						    !
    !  do i_dim=1,BLOCK_DIMENSION				    !
    !    dos(s)=dos(s)-(1.0d0/pi)*aimag(Green(l,l,i_dim,i_dim))	    !
    !  enddo							    !
    !enddo							    !
    !write(300,*)real(epsi(s)),dos(s)				    !
    !****************** End DoS ************************************!

    !------ Green Function Matrix: Off-Diagonal elements ---------------
    do i=2,BLOCK_NUMBER
      do l=1,i-1
        Green(l,i,:,:)=matmul(Green(l,i-1,:,:),sim_delta%R(i-1,:,:))
      enddo !l
    enddo !i

   !--- Transmittances between real channels -------------------------------
    do i=1,CHANNEL_NUMBER
      l=1+int((chan_pos(i)-1)/BLOCK_DIMENSION)
      n=1+mod(chan_pos(i)-1,BLOCK_DIMENSION)
      do j=i,CHANNEL_NUMBER
	k=1+int((chan_pos(j)-1)/BLOCK_DIMENSION)
	m=1+mod(chan_pos(j)-1,BLOCK_DIMENSION)
        if (i .eq. j ) then
	  T_LL(i,j)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(sigma_L(i,s)))*	&
		   Green(l,k,n,m)-1.d0))**2 -1.d0
        elseif (i .lt. j) then
	  T_LL(i,j)=4.d0*dabs(aimag(sigma_L(i,s))*aimag(sigma_L(j,s)))*cdabs(Green(l,k,n,m))**2
          T_LL(j,i)=T_LL(i,j)
	endif
      enddo !j
    enddo !i
    
    if (decoherence .eq. 'Y') then
       !-- Transmittances between real and decoherent channels ----------------
       do i=1,CHANNEL_NUMBER
	 l=1+int((chan_pos(i)-1)/BLOCK_DIMENSION)
         n=1+mod(chan_pos(i)-1,BLOCK_DIMENSION)
         do j=1,DEC_NUMBER
	   k=1+int((dec_pos(j)-1)/BLOCK_DIMENSION)
	   m=1+mod(dec_pos(j)-1,BLOCK_DIMENSION)
	   if (l .le. k) then  !!!!    Only works if Tij  = Tji
             T_LD(i,j)=4.d0*dabs(aimag(sigma_L(i,s))*aimag(sigma_D(j)))*cdabs(Green(l,k,n,m))**2
	   else 
             T_LD(i,j)=4.d0*dabs(aimag(sigma_L(i,s))*aimag(sigma_D(j)))*cdabs(Green(k,l,m,n))**2
           endif
         enddo !j
       enddo !i
      
       !--- Transmittances between decoherent channels -------------------------------
       do i=1,DEC_NUMBER
	 l=1+int((dec_pos(i)-1)/BLOCK_DIMENSION)
         n=1+mod(dec_pos(i)-1,BLOCK_DIMENSION)
         do j=i,DEC_NUMBER
	   k=1+int((dec_pos(j)-1)/BLOCK_DIMENSION)
	   m=1+mod(dec_pos(j)-1,BLOCK_DIMENSION)
           if (i .eq. j ) then
	     T_DD(i,j)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(sigma_d(i)))*	&
		   Green(l,k,n,m)-1.d0))**2 -1.d0
           elseif (i .lt. j) then
	     T_DD(i,j)=4.d0*dabs(aimag(sigma_d(i))*aimag(sigma_d(j)))*cdabs(Green(l,k,n,m))**2
             T_DD(j,i)=T_DD(i,j)
           endif
         enddo !j
       enddo !i

       !---- Effective transmittance ------------------------------------------------------
       call invert_real(DEC_NUMBER,T_DD)
       T_ef= T_LL - matmul(T_LD, matmul(T_DD,transpose(T_LD)) )

        write(77,1177)real(epsi(s)),(T_ef(1,j),j=2,channel_number)
        write(79,1177)real(epsi(s)),sum(T_ef(1,:))-T_ef(1,1)

   else !decoherence .ne. Y
	write(77,1177)real(epsi(s)),(T_LL(1,j),j=2,channel_number)
        write(79,1177)real(epsi(s)),sum(T_LL(1,:))-T_LL(1,1)
   endif 

    !--- Advance info --------------------------------------------------!*
    do n=20,100,20							!*
      if (n .eq. int(real(s)*100./real(RESOLUTION))) then		!*
        if (h .eq. n) then						!*
		goto 118						!*
        else								!*
		print *, n,'% completado'				!*
        	h=n							!*
        endif								!*
      endif								!*
    enddo								!*
   118 continue		


!!! PRUEBA:::::: CORRECCION DE ERROR EN ENER POR SIGMA :::::::........ !!
    do i=1,CHANNEL_NUMBER 						!
      m=1+mod(chan_pos(i)-1,BLOCK_DIMENSION)				!
      k=1+int((chan_pos(i)-1)/BLOCK_DIMENSION)				!
      ener(k,m,m)=ener(k,m,m)- sigma_L(i,s)				!
    enddo !i								!
!!!!!!! FIN PRUEBA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enddo !s
 1177 format(50E16.8)
end subroutine DAmatoPastawski

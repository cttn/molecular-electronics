!################################################################################!
!#  FALTANTES A OPTIMIZAR y/o AGREGAR :                     			#!
!#     * Calculos cuando T(i,j) .ne. T(j,i)                 			#!
!#     * Calculos con matrices matv,math,matd a eliminar         		#!
!#     * Hay un Lead por canal!                         			#!
!#     * Revertir Locacion de memoria, especialmente subrutinas.     	        #! 
!#     * Diferenciar el caso matriz T simétrica o no y elegir subrut            #!
!#        eb funcion de eso.                         			        #! 
!################################################################################!
module MultiDP_global
  implicit none
  integer,save :: BLOCK_DIMENSION, BLOCK_NUMBER, RESOLUTION, CHANNEL_NUMBER, 		&
     DEC_NUMBER
  integer :: i,j,m,s,k,l,H,N,M_comp=0,M_real=0
  integer,allocatable,save :: chan_pos(:),dec_pos(:)
  complex(kind(0d0)),allocatable,save :: ener(:,:,:),hop(:,:,:),             		&
     Sigma_D(:),epsi(:)                        
  real(kind(0d0)),allocatable,save :: ener_lead(:),hop_lead(:),tun_lead(:)
  real(kind(0d0)),parameter :: pi=dacos(-1.d0)
  character,save :: decoherence,arch_trans*20
  contains
!%%% Sigma LEAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  function Sigma_L(var,e)
   complex(kind(0d0)) :: Sigma_L
   integer :: e,var
     if    ((real(epsi(e))-ener_lead(var)) .gt. 2.0d0*dabs(hop_lead(var))) then 
     Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0                           	&
        -dsqrt(((real(epsi(e))-ener_lead(var))/2.0d0)**2-hop_lead(var)**2),0.0d0)   
     elseif(dabs(real(epsi(e))-ener_lead(var)) .le. dabs(2.0d0*hop_lead(var))) then
     Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0,-dsqrt(hop_lead(var)**2    &
        -((real(epsi(e))-ener_lead(var))/2.0d0)**2))
     elseif((real(epsi(e))-ener_lead(var)) .lt. -2.0d0*dabs(hop_lead(var))) then 
     Sigma_L=dcmplx((real(epsi(e))-ener_lead(var))/2.0d0                         	&
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

!--- Interfase ------------------------------------------------------------------!
  decoherence='Y' !Yes/NO
  arch_trans='mtdec.dat'

  RESOLUTION=150        !Number of energy points
  allocate(epsi(RESOLUTION))


  BLOCK_DIMENSION=2      !Dimension of blocks
  BLOCK_NUMBER=704       !Number of blocks
  allocate(ener(BLOCK_NUMBER,BLOCK_DIMENSION,BLOCK_DIMENSION))
  if (block_number .gt. 1) then
   allocate(hop(BLOCK_NUMBER-1,BLOCK_DIMENSION,BLOCK_DIMENSION))
  endif


  CHANNEL_NUMBER=2       !Numbers of real channels
  allocate(chan_pos(CHANNEL_NUMBER))
  allocate(ener_lead(CHANNEL_NUMBER))
  allocate(hop_lead(CHANNEL_NUMBER))
  allocate(tun_lead(CHANNEL_NUMBER))
    chan_pos(1)=1
    chan_pos(2)=BLOCK_NUMBER*BLOCK_DIMENSION !200

  if (decoherence .eq. 'Y') then
    DEC_NUMBER=BLOCK_NUMBER*BLOCK_DIMENSION   !Number of decoherent channels
    allocate(dec_pos(DEC_NUMBER))
    allocate(Sigma_D(DEC_NUMBER))
      do i=1,DEC_NUMBER
        dec_pos(i)=i
      enddo
    !   dec_pos(1)=1
    !   dec_pos(2)=10
  endif

  open(77, file=arch_trans)

!  call asignations
  call asignaciones
  call DAmatoPastawski

  close(77)
end program MultiDP


subroutine asignations
  use MultiDP_global
  implicit none
  real(kind(0d0)) :: kT,Emin,Emax,eta

  kT=0.025d0 !eV
  Emin=-10.0d0 
  Emax=10.0d0 
  eta=1.d-6

  do s=1,RESOLUTION
   epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s,kind(0d0))/real(RESOLUTION,kind(0d0)),eta)
  enddo  

!--- Leads asignations -------------------------------------
  do i=1,CHANNEL_NUMBER
    ener_lead(i)=0.d0
    hop_lead(i)=3.d0
    tun_lead(i)=3.d0
  enddo

!--- Decs asignations --------------------------------------
  if (decoherence .eq. 'Y') then
    do i=1,DEC_NUMBER
      Sigma_D(i)=dcmplx(0.d0,-kT) 
    enddo
  endif

!--- Hamiltonian -------------------------------------------
  hop=dcmplx(0.d0,0.d0)
  ener=dcmplx(0.d0,0.d0)
  do i =1,BLOCK_NUMBER
    do j=1,BLOCK_DIMENSION
      do m=1,BLOCK_DIMENSION
       if (j .eq. m) then
      ener(i,j,m)=dcmplx(0.d0,0.d0)
      if (i .lt. BLOCK_NUMBER) hop(i,j,m) =dcmplx(0.d0,0.d0)
       elseif( j .eq. m-1) then
      ener(i,j,m)=dcmplx(5.d0,0.d0)
      if (i .lt. BLOCK_NUMBER) hop(i,j,m) =dcmplx(5.d0,0.d0)
       elseif (j .eq. m+1) then
      ener(i,j,m)=dcmplx(5.d0,0.d0)
      if (i .lt. BLOCK_NUMBER) hop(i,j,m) =dcmplx(5.d0,0.d0)
       endif
      enddo !m
    enddo !j
  enddo !i

end subroutine asignations

subroutine DAmatoPastawski
  use MultiDP_global
  implicit none
  integer :: i_dim,j_dim
  real(kind(0d0)),allocatable :: ID(:,:),dos(:),T_LL(:,:),T_LD(:,:),T_DD(:,:),     &
     T_ef(:,:)
  complex(kind(0d0)),allocatable :: matv(:,:),math(:,:),matd(:,:),inver(:,:),     &
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
      inver(:,:)=epsi(s)*ID(:,:)-ener(BLOCK_NUMBER+2-l,:,:)     &
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
    do l=1,BLOCK_NUMBER                             !
      do i_dim=1,BLOCK_DIMENSION                     !
        dos(s)=dos(s)-(1.0d0/pi)*aimag(Green(l,l,i_dim,i_dim))         !
      enddo                                 !
    enddo                                 !
    write(300,*)real(epsi(s)),dos(s)                     !
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
       T_LL(i,j)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(sigma_L(i,s)))*    &
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
          T_DD(i,j)=(cdabs(dcmplx(0.0d0,1.0d0)*2.0d0*dabs(aimag(sigma_d(i)))*    &
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

       write(77,1777)real(epsi(s)),T_ef(1,2)

   else !decoherence .ne. Y
     write(77,*)real(epsi(s)),T_LL(1,2)
   endif 

    !--- Advance info --------------------------------------------------!*
    do n=20,100,20                             !*
      if (n .eq. int(real(s)*100./real(RESOLUTION))) then         !*
        if (h .eq. n) then                         !*
         goto 118                        !*
        else                                 !*
         print *, n,'% completado'                !*
             h=n                            !*
        endif                                 !*
      endif                                 !*
    enddo                                 !*
   118 continue         


!!! PRUEBA:::::: CORRECCION DE ERROR EN ENER POR SIGMA :::::::........ !!
    do i=1,CHANNEL_NUMBER                          !
      m=1+mod(chan_pos(i)-1,BLOCK_DIMENSION)                 !
      k=1+int((chan_pos(i)-1)/BLOCK_DIMENSION)                 !
      ener(k,m,m)=ener(k,m,m)- sigma_L(i,s)                 !
    enddo !i                                 !
!!!!!!! FIN PRUEBA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enddo !s
 1777 format(50E16.8)
end subroutine DAmatoPastawski


subroutine asignaciones
  use MultiDP_global
  implicit none
  integer :: numanillos,indice,jo,cont,contbenceno,contquinoide,unidad,         &
     i_toleran,semilla,tolerancia,dist_qui
  integer,allocatable :: registro(:)
  real(kind(0d0)) :: ran3,kT,Emin,Emax,eta,porc_max,porc_min,porcentaje
  complex(kind(0d0)) :: hop_Npara,hop_Npara_q,ener_N,ener_Nmas,enerpara,enerorto,     &
     vop_qui,voo_qui,vop_ben,voo_ben

  kT=0.025d0 !eV
  Emin=-10.0d0 
  Emax=10.0d0 
  eta=1.d-8

  do s=1,RESOLUTION
   epsi(s)=dcmplx(Emin+(Emax-Emin)*real(s,kind(0d0))/real(RESOLUTION,kind(0d0)),eta)
  enddo  

!--- Leads asignations -------------------------------------
  do i=1,CHANNEL_NUMBER
    ener_lead(i)=0.d0
    hop_lead(i)=5.d0
    tun_lead(i)=5.d0
  enddo

!--- Decs asignations --------------------------------------
  if (decoherence .eq. 'Y') then
    do i=1,DEC_NUMBER
      Sigma_D(i)=dcmplx(0.d0,-kT) 
    enddo
  endif


!--- Energias de sitio y hoppings de todos los sitios en BL de PAni-HCl
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

!-- wrap -----------------------------!
  semilla=118171
  numanillos=201
  indice=(numanillos-1)*7/2+4
  porc_max=0.26d0 !Topes para el criterio de aceptacion de configuraciones
  porc_min=0.24d0
  tolerancia= 100000 !Cantidad de pruebas buscando configuraciones
  dist_qui=1 !Distamcia minima entre quinoides en unidades de bencenos
  allocate(registro(numanillos))


  i_toleran=0
 8900 continue !si el sorteo no fue satisfactorio

  cont=0
  contbenceno=0
  contquinoide=0

  h=dist_qui
  do j=1,numanillos
        if (h .lt. dist_qui) then !Si h < dist_qui, pongo un benceno
             registro(j)=0 !benceno
             contbenceno = contbenceno + 1
                 cont=cont+1
         h=h+1
     else !Si h >= dist_qui, sorteo
        unidad = 1 + int(3*ran3(semilla))
           if (unidad .eq. 1) then 
         registro(j)=1 !quinoide
         contquinoide=contquinoide + 1
         cont=cont+1
                h=0
           else
         registro(j)=0 !benceno
             contbenceno = contbenceno + 1
             cont=cont+1
         h=h+1
           endif
       endif
  enddo

! *** Evaluacion del sorteo ******************************
  i_toleran=i_toleran+1
  if (i_toleran .ge. tolerancia) then
    print *
    print *, 'No se encontraron configuraciones que satisfagan ',         &
     'los porcentajes requridos dentro del limite de tolerancia ',    &
     'elegido: ',tolerancia
    stop 'bluff'
  endif

  porcentaje=real(contquinoide)/real(cont)

  if (porcentaje .gt. porc_max .or. porcentaje .lt. porc_min) then
     goto 8900
  endif

!--- Asignacion de lugares en base al sorteo *************
   jo=0
   do m=1,numanillos
      if (registro(m) .eq. 0) then !benceno         !!!!!!!!!!!!!!!!
     if (real(m)/2.0d0 .eq. real(m/2) ) then !jo par

     !-- 1 --
         ener(jo+1,1,1)=enerpara
         ener(jo+1,1,2)=vop_ben
         ener(jo+1,2,1)=vop_ben
         ener(jo+1,2,2)=enerorto

         hop(jo+1,1,1)=vop_ben
         hop(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+1,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+1,2,2)=voo_ben

     !-- 2 --
         ener(jo+2,1,1)=enerorto
         ener(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,2)=enerorto

         hop(jo+2,1,1)=voo_ben
         hop(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,2)=vop_ben

     !-- 3 --
         ener(jo+3,1,1)=enerorto
         ener(jo+3,1,2)=vop_ben
         ener(jo+3,2,1)=vop_ben
         ener(jo+3,2,2)=enerpara

         if ((jo+3) .ne. indice) then
            hop(jo+3,1,1)=dcmplx(0.0d0,0.0d0)
            hop(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
            hop(jo+3,2,1)=hop_Npara
            hop(jo+3,2,2)=dcmplx(0.0d0,0.0d0)
         endif
  
     jo = jo+3
     else !jo impar

     !-- 1 --
               if ( (m .gt. 1) .and. (registro(m-1) .eq. 1)) then
             ener(jo+1,1,1)=ener_Nmas
           else
             ener(jo+1,1,1)=ener_N
           endif
         ener(jo+1,1,2)=hop_Npara
         ener(jo+1,2,1)=hop_Npara
         ener(jo+1,2,2)=enerpara

         hop(jo+1,1,1)=dcmplx(0.0d0,0.0d0)    
         hop(jo+1,1,2)=dcmplx(0.0d0,0.0d0) 
         hop(jo+1,2,1)=vop_ben 
         hop(jo+1,2,2)=vop_ben 

     !-- 2 --
         ener(jo+2,1,1)=enerorto
         ener(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,2)=enerorto

         hop(jo+2,1,1)=voo_ben
         hop(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,2)=voo_ben !vop_ben

     !-- 3 --
         ener(jo+3,1,1)=enerorto
         ener(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+3,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+3,2,2)=enerorto

         hop(jo+3,1,1)=vop_ben
         hop(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+3,2,1)=vop_ben
         hop(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

     !-- 4 --
         ener(jo+4,1,1)=enerpara
         ener(jo+4,1,2)=hop_Npara
         ener(jo+4,2,1)=hop_Npara
             if ( (m .lt. numanillos) .and. (registro(m+1) .eq. 1) ) then
          ener(jo+4,2,2)=ener_Nmas
         else
          ener(jo+4,2,2)=ener_N
         endif

         if ((jo+4) .ne. indice) then
          hop(jo+4,1,1)=dcmplx(0.0d0,0.0d0) !vop_ben
          hop(jo+4,1,2)=dcmplx(0.0d0,0.0d0)
           if ( (m .lt. numanillos) .and. (registro(m+1) .eq. 1) ) then
            hop(jo+4,2,1)=hop_Npara_q
           else
            hop(jo+4,2,1)=hop_Npara
           endif
          hop(jo+4,2,2)=dcmplx(0.0d0,0.0d0)
         endif

     jo = jo+4
     endif

     contbenceno = contbenceno + 1
     
  elseif (registro(m) .eq. 1) then !quinoide  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (real(m)/2.0d0 .eq. real(m/2)) then !jo par

     !-- 1 --
         ener(jo+1,1,1)=enerpara
         ener(jo+1,1,2)=vop_qui
         ener(jo+1,2,1)=vop_qui
         ener(jo+1,2,2)=enerorto 

         hop(jo+1,1,1)=vop_qui
         hop(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+1,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+1,2,2)=voo_qui

     !-- 2 --
         ener(jo+2,1,1)=enerorto
         ener(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,2)=enerorto

         hop(jo+2,1,1)=voo_qui
         hop(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,2)=vop_qui !!!!!!!!!!!!!!!!!!!!!!

     !-- 3 --
         ener(jo+3,1,1)=enerorto
         ener(jo+3,1,2)=vop_qui
         ener(jo+3,2,1)=vop_qui
         ener(jo+3,2,2)=enerpara

         if ((jo+3) .ne. indice) then
         hop(jo+3,1,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+3,2,1)=hop_Npara_q
         hop(jo+3,2,2)=dcmplx(0.0d0,0.0d0)
           endif
         
     jo = jo+3
             
     else !jo impar
         
     !-- 1 --    
             ener(jo+1,1,1)=ener_Nmas
         ener(jo+1,1,2)=hop_Npara_q
         ener(jo+1,2,1)=hop_Npara_q
         ener(jo+1,2,2)=enerpara

         hop(jo+1,1,1)=dcmplx(0.0d0,0.0d0)    !siempre los del lado derecho
         hop(jo+1,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+1,2,1)=vop_qui
         hop(jo+1,2,2)=vop_qui

     !-- 2 --
         ener(jo+2,1,1)=enerorto
         ener(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+2,2,2)=enerorto

         hop(jo+2,1,1)=voo_qui
         hop(jo+2,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,1)=dcmplx(0.0d0,0.0d0)
         hop(jo+2,2,2)=voo_qui !vop_qui

     !-- 3 --
         ener(jo+3,1,1)=enerorto
         ener(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
         ener(jo+3,2,1)=dcmplx(0.0d0,0.0d0)
         ener(jo+3,2,2)=enerorto

         hop(jo+3,1,1)=vop_qui
         hop(jo+3,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+3,2,1)=vop_qui
         hop(jo+3,2,2)=dcmplx(0.0d0,0.0d0)

     !-- 4 --
         ener(jo+4,1,1)=enerpara
         ener(jo+4,1,2)=hop_Npara_q
         ener(jo+4,2,1)=hop_Npara_q
             ener(jo+4,2,2)=ener_Nmas

         if ((jo+4) .ne. indice) then
         hop(jo+4,1,1)=dcmplx(0.0d0,0.0d0)!vop_qui
         hop(jo+4,1,2)=dcmplx(0.0d0,0.0d0)
         hop(jo+4,2,1)=hop_Npara !hop_Npara_q
         hop(jo+4,2,2)=dcmplx(0.0d0,0.0d0)
         endif


     jo = jo+4

     endif
     
     contquinoide=contquinoide + 1

  endif
  
  enddo !m

!*********** EVALUACION DE ASIGNACIONES *************************

!---- Registro de configuracion ------------------
  open(unit=8181, file='configuracion.dat')
   do j=1,numanillos
    write(8181,*)registro(j)
   enddo
  close(8181)

end subroutine asignaciones

!*********************************************************
!*      Generador de numeros aleatorios                  *
!*********************************************************
   FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
!     PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
 11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12        continue
 13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
   END FUNCTION
!  (C) Copr. 1986-92 Numerical Recipes Software.

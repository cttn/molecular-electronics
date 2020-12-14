!##################################################################################!
module MOD_nmrFA
    use mtmod, only: grnd
    implicit none
    contains

    subroutine svd_reconst(M,N,A,rank)
        !*** Propósito: Aproxima una matriz con SVD hasta rango "rank"
        implicit none
        external DGESVD
        integer,        intent(in)    :: M, N, rank
        real(kind(0d0)),intent(inout) :: A(M,N)
        integer     :: INFO, LWORK, min_MN, i,j,k
        real(kind(0d0)),allocatable :: S(:), work(:)
        real(kind(0d0)) :: U(M,M), VT(N,N)

        !-- Alocando el vector de valores singulares
        min_MN = M; if (M .gt. N) min_MN = N; allocate(S(min_MN))
        if(rank .gt. min_MN) stop 'El rango no puede ser meyor que la menor dimension de la matriz a descomponer!'

        !-- Buscando lwork apropiado
        lwork = -1; allocate(work(1)); call DGESVD('S','S',M,N,A,M,S,U,M,VT,N,WORK,LWORK,INFO)
        call check_lapack_info(INFO,'svd_reconst, calculando lwork')
        lwork = int(work(1)); deallocate(work)
        
        !-- Calculando SVD
        allocate(work(lwork)); call DGESVD('S','S',M,N,A,M,S,U,M,VT,N,WORK,LWORK,INFO)   
        call check_lapack_info(INFO,'svd_reconst, calculando SVD')

        !-- Reconstruccion de la matriz hasta rango "rank"        
        A = 0.d0
        do i=1,M
          do j=1,N
            do k=1,rank
              A(i,j)= A(i,j) + s(k)*U(i,k)*VT(k,j)
            enddo
          enddo
        enddo
    end subroutine svd_reconst

    subroutine matrix_gen_1(M,N,A)
        !*** Propósito: Generar una matriz de las dimensiones pedidas
        implicit none
        integer,intent(in) :: M,N
        real(kind(0d0)),intent(out) :: A(M,N)
        integer :: i,j
        real(kind(0d0)) :: ri, rj
        do i=1,M
          do j=1,N
            ri=real(i,kind(0d0))
            rj=real(j,kind(0d0))
            A(i,j)=ri*rj-dsqrt(ri+rj) -ri + 5.d0
          enddo
        enddo
    end subroutine matrix_gen_1

    subroutine check_matrix_dif(M,N,A,B,max_val)
        !*** Propósito: Realiza la diferencia entre dos matrices de iguales 
        ! dimensiones y Entrega el elemento de mayor valor absoluto.
        implicit none
        integer,intent(in)          :: M, N
        real(kind(0d0)),intent(in)  :: A(M,N), B(M,N)
        real(kind(0d0)),intent(out) :: max_val
        real(kind(0d0)) :: dif(M,N)
        integer :: i,j
        do i = 1,M
          do j = 1,N
            dif(i,j) = dabs(A(i,j) - B(i,j))
          enddo
        enddo; max_val = maxval(dif)
    end subroutine check_matrix_dif

    subroutine import_nmr_data(M,N,A,archivo)
        !*** Propósito: Importar los datos de RMN a una matriz tratable por el
        ! programa
        implicit none
        integer,intent(in) :: M,N !Numero de filas y de columnas a leer
        real(kind(0d0)),intent(out) :: A(M,N) !Matriz de datos a rellenar
        character(len=*) :: archivo
        integer :: i,j, iostatus,iostatus_open
         !-- Leyendo datos NMR: El metodo ALS la necesita MxN con:
         !     M dimension del # de experimentos y N del numero de frecuencias
        print *, 'Imporatando datos: ', trim(adjustl(archivo))
        open(24, file = trim(adjustl(archivo)), iostat=iostatus_open, err=10024, status='old')
           do i=1,M 
             read(24,1124,iostat=iostatus) (A(i,j), j=1,N)
             if (iostatus .gt. 0) then
                 print *, 'Error de lectura. Archivo: ', archivo
                 print *, 'Codigo de Error: ', iostatus
                 stop
             elseif( iostatus .lt. 0 ) then
                 print *, 'Error de lectura. Archivo: ', archivo
                 print *, 'Final de archivo imprevisto. Línea: ', i
             endif
           enddo
           !-- Zonceras y zandeces
           write(*,'(A16,I4.4,A16,I4.4,A26)') '    --- Leídos ',i-1, ' experimentos y ', j-1, ' frecuencias en los datos.'
           !read(24,1124,iostat=iostatus)
           !read(24,1124,iostat=iostatus); if (iostatus .ge. 0) print *, '   --- Al parecer hay mas datos en el archivo.'
        close(24)
        10024 if (iostatus_open .ne. 0) then
                print *, 'Error de I/O al tratar de abrir el archivo: ', trim(adjustl(archivo))
                if (iostatus_open .eq. 2) print *, 'No se encontró el archivo...'
                stop 'Corroborar que el archivo de datos indicado sea el correcto'
              endif
     1124 format(4096E16.8)
    end subroutine import_nmr_data

    subroutine ALS_factor_analysis(M,N,nmr_svd,K,mat_S,mat_L)
        !*** Propósito: Implementar el método ALS [Anal. Chem. 76, 1982 (2004)]
        ! como una tecnica de "Factor Analysis" para NMR
        implicit none
        integer,intent(in) :: M,N,K !K es el número de factores propuesto, igual al rango (rank) usado para SVD
        real(kind(0d0)),intent(in)  :: nmr_svd(M,N) !Matriz de datos de RMN aproximada por SVD
        real(kind(0d0)),intent(out) :: mat_S(M,K),mat_L(K,N) !Matriz de concentraciones y composiciones respectivamente
        integer :: i,j,cont, tolerancia, ndig, ncont, nequal
        logical :: convergence_status=.false.
        real(kind(0d0)) :: val,prev_val

        !print *, 'ALS no implementado 100%: Aun no hay un buen criterio de convergencia... '

        !-- Setting del criterio de convergencia
        Ndig = 2            !# de digitos de la precision buscada para converger
        Nequal =  1000       !# de veces que debo obtener el mismo val para converger
        tolerancia = 100000  !# de veces que debo correr antes de desistir si no converge
        val = 0.d0          
        cont = 0            !# de veces que corrió el loop que trata de converger
        ncont = 0           !# de veces que encuentra el mismo val

        !-- Guess inicial para L
        do i=1,K
          do j=1,N
            mat_L(i,j)=grnd()
          enddo
        enddo
        call upgrade_rank(K,N,mat_L,K)
        call solve_eq_S(M,N,K,nmr_svd,mat_S,mat_L)      !Resolver D=SL para S

        do while(.not. convergence_status)
          cont = cont + 1
          call replace_negative_zero(M,K,mat_S)         !Condicion para S
          call upgrade_rank(M,K,mat_S,K)                !Subir el rango si bajo
          call solve_eq_L(M,N,K,nmr_svd,mat_S,mat_L)    !Resolver D=SL para L
          call apply_L_condition(K,M,mat_L)             !Condicion para L
          call upgrade_rank(K,N,mat_L,K)                !Subir el rango si bajo
          call solve_eq_S(M,N,K,nmr_svd,mat_S,mat_L)    !Resolver D=SL para S      

          !Evaluar convergencia por cuadrados minimos
          prev_val=val
          call check_matrix_dif(M,N,nmr_svd,matmul(mat_S,mat_L),val)
          
          if ( int(val*(10.d0**Ndig)) .eq. int(prev_val*(10.d0**Ndig)) ) then
              ncont = ncont + 1
          else
              ncont = 0
          endif

          if (ncont .ge. Nequal) then
             print *, '   --- Convergencia!'
             print *, '   --- Precision de la solucion: ', val
             print *, '   --- Numero de iteraciones: ', cont
             convergence_status = .true.
          endif

          if (cont .ge. tolerancia) then
             print *, '   --- No hubo Convergencia... Umbral de tolerancia excedido'
             print *, '   --- Precision de la solucion: ', val
             print *, '   --- Numero de iteraciones: ', cont
             convergence_status = .true.
          endif
       enddo

!        call solve_eq_L(M,N,K,nmr_svd,mat_S,mat_L)      !Resolver D=SL para L
!        call solve_eq_S(M,N,K,nmr_svd,mat_S,mat_L)      !Resolver D=SL para S
    end subroutine ALS_factor_analysis

    subroutine solve_eq_L(M,N,K,D,S,X)
        !*** Propósito: Utiliza Lapack para resolver X en: D=SX (B=AX en lapack)
        implicit none
        external DGELS
        integer,intent(in) :: M, N, K
        real(kind(0d0)),intent(in)  :: D(M,N), S(M,K)   !!Matrices D y S respectivamente
        real(kind(0d0)),intent(out) :: X(K,N)           !!Matriz L
        integer :: LDB, LWORK, INFO, i
        real(kind(0d0)),allocatable :: work(:)
        real(kind(0d0)) :: B(M,N), A(M,K)               !!Para no sobreescribir los originales
        if(M .gt. N) stop 'La matriz de datos debe ser de dimension MxN con M < N'
        if(K .gt. M) stop 'El rango elegido para SVD, K, nunca puede ser mayor que el numero de filas de la matriz de datos'
        LDB = M; B(:,:) = D(:,:); A(:,:) = S(:,:)

        !-- Calculando lwork
        LWORK = -1; allocate(work(1)); call DGELS('N',M,K,N,A,M,B,LDB,WORK,LWORK,INFO)
        call check_lapack_info(INFO,'solve_eq_L, calculando lwork')       
        lwork = int(work(1)); deallocate(work)

        !-- Calculando X (Primeras K filas de B)
        allocate(work(lwork))
        call DGELS('N',M,K,N,A,M,B,LDB,WORK,LWORK,INFO)
        call check_lapack_info(INFO,'solve_eq_L, calculando X')       

        !-- Rearmando la solucion      
        do i=1,K; X(i,:) = B(i,:); enddo
    end subroutine solve_eq_L

    subroutine solve_eq_S(M,N,K,D,X,L)
        !*** Propósito: Utiliza Lapack para resolver X en: D=XL (B=AX en lapack)
        ! Debe transponerse la ecuacion como DT = LT XT 
        implicit none
        external DGELS
        integer,intent(in) :: M, N, K
        real(kind(0d0)),intent(in)  :: D(M,N), L(K,N)   !! Matrices D y L respectivamente
        real(kind(0d0)),intent(out) :: X(M,K)           !! Resultado: S
        integer :: LDB, LWORK, INFO, i, j
        real(kind(0d0)),allocatable :: work(:)
        real(kind(0d0)) :: B(N,M), A(N,K)               !!Transpuestas de D y L respectivamente

        if(M .gt. N) stop 'La matriz de datos debe ser de dimension MxN con M < N'
        if(K .gt. M) stop 'El rango elegido para SVD, K, nunca puede ser mayor que el numero de filas de la matriz de datos'
        LDB = N
        do i = 1,N
          do j = 1,M
            B(i,j) = D(j,i)
          enddo
          do j=1,K
            A(i,j) = L(j,i)
          enddo
        enddo

        !-- Calculando lwork
        LWORK = -1; allocate(work(1)); call DGELS('N',N,K,M,A,N,B,LDB,WORK,LWORK,INFO)
        lwork = int(work(1)); deallocate(work)
        call check_lapack_info(INFO,'solve_eq_S, calculando lwork')

        !-- Calculando X (Primeras K filas de B, luego transponer)
        allocate(work(lwork)); call DGELS('N',N,K,M,A,N,B,LDB,WORK,LWORK,INFO)
        call check_lapack_info(INFO,'solve_eq_S, calculando X')

        !-- Reordenando (transoniendo) la solucion
        do i=1,K
          do j=1,M
            X(j,i) = B(i,j)
          enddo
        enddo
   end subroutine solve_eq_S

   subroutine check_lapack_info(info,descrip)
       !*** Propósito: chequear variable INFO de las subrutinas de lapack... en desarrollo...
       implicit none
       integer,intent(in) :: INFO            !INFO de lapack
       character(len=*),optional :: descrip !Cadena opcional, descripción del ugar del codigo en la que se hace el chequeo.
       !Se podria agregar un input opcional: Subrutina de lapack a chequear, para mejor debug a partir de INFO...
       if (INFO .ne. 0) then
           print *
           print *, 'Error en lapack... INFO=',INFO
           if (present(descrip)) then
               print *, 'Ubicacion: ', trim(adjustl(descrip))
           else
               print *, 'No se brindaron datos de la ubicacion'
           endif
           stop
        endif
    end subroutine check_lapack_info

    subroutine replace_negative_zero(M,N,A)
        !*** Propósito: Toma una matriz cualquiera y reemplaza sus elementos negativos por ceros
        implicit none
        integer,intent(in) :: M,N
        real(kind(0d0)),intent(inout) :: A(M,N)
        integer :: i,j
        do i=1,M
          do j=1,N
            if (A(i,j) .lt. 0.0d0) A(i,j) = 0.d0
          enddo
        enddo
    end subroutine replace_negative_zero

    subroutine apply_L_condition(M,N,L)
        !*** Porpoósito: Aplicar la condicion particular estipulada para L
        implicit none
        integer,intent(in) :: M, N 
        real(kind(0d0)),intent(inout) :: L(M,N)
        integer :: i,j,jmax
        real(kind(0d0)) :: ref, frac

        frac = 0.25d0 !! Fraccion del elemento con mayor valor absoluto a usar.

        do i=1,M
          !Para cada fila i, busco el elemento de mayor valor absoluto
          ref  = 0.d0
          jmax = 1
          do j=1,N 
            if ( dabs(L(i,j)) .gt. ref ) then
                ref  = dabs(L(i,j))
                jmax = j
            endif
          enddo

          !Si ese elemento es negativo, cambio el signo de la fila
          if ( L(i,jmax) .lt. 0.d0 ) L(i,:) = -L(i,:)

          !Todos los elementos menores que una fraccion del mayor, se hacen cero
          do j = 1,N
             if ( L(i,j) .lt. frac*L(i,jmax) ) L(i,j) = 0.d0
          enddo
        enddo
    end subroutine apply_L_condition

    subroutine upgrade_rank(M,N,A,rnk)
        !*** Propósito: Si l rango de A es menor que rnk, se sube el rango de la
        !  matriz A hasta rnk
        implicit none
        integer,intent(in) :: M, N, rnk
        real(kind(0d0)),intent(inout) :: A(M,N)
        integer :: val, i, j

        call check_rank(M,N,A,val)
        do while (val .lt. rnk)
           do i=1,M
             do j=1,N
               A(i,j) = A(i,j) + grnd()
             enddo
           enddo
           call check_rank(M,N,A,val)
        enddo
    end subroutine upgrade_rank

    subroutine check_rank(M,N,mat,val)
        !*** Propósito: Calcula el rango de una matriz dada, guarda el resultado
        ! en "val"
        implicit none
        external DGESVD
        integer,intent(in) :: M,N
        real(kind(0d0)),intent(in) :: mat(M,N)
        integer,intent(out) :: val
        integer     :: INFO, LWORK, min_MN, i,cont
        real(kind(0d0)),allocatable :: S(:), work(:)
        real(kind(0d0)) :: U(M,M), VT(N,N), A(M,N), umbral

        umbral = 1.d-1 !Umbral para decidir si un valor es o no cero
        A(:,:) = mat(:,:) !! Para no destruir los datos de mat

        !-- Alocando el vector de valores singulares
        min_MN = M; if (M .gt. N) min_MN = N; allocate(S(min_MN))

        !-- Buscando lwork apropiado
        lwork = -1; allocate(work(1)); call DGESVD('N','N',M,N,A,M,S,U,M,VT,N,WORK,LWORK,INFO)
        call check_lapack_info(INFO,'check_rank, calculando lwork')
        lwork = int(work(1)); deallocate(work)
        
        !-- Calculando SVD
        allocate(work(lwork)); call DGESVD('N','N',M,N,A,M,S,U,M,VT,N,WORK,LWORK,INFO)   
        call check_lapack_info(INFO,'check_rank, calculando SVD')

        !Contando valores singulares distintos de cero
        cont=0
        do i=1,min_MN
          if ( dabs(s(i)) .gt. umbral) cont = cont + 1
        enddo

        val = cont
        !print *, val
    end subroutine check_rank
end module MOD_nmrFA




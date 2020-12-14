program nmrFA
    use MOD_nmrFA
    implicit none
    integer, parameter :: M=8, N=1000
    integer :: rank, i,j
    real(kind(0d0)) :: nmr_or(M,N),nmr_svd(M,N), val
    real(kind(0d0)),allocatable :: mat_S(:,:),mat_L(:,:)
    character(len=50) :: sol

    !-- Import NMR DATA ----------------------------------
    call import_nmr_data(M,N,nmr_or,'NMR_spectra_2012.dat')

    !Normalizando un poco el eje y
    nmr_or=nmr_or/1.d4

    !-- Sacando los valores negativos de la matriz de datos
    call replace_negative_zero(M,N,nmr_or)

    !-- Chequeando el rango de la matriz de datos
    call check_rank(M,N,nmr_or,rank)
    print *, '   --- El rango de la matriz de datos es: ', rank

    !-- SVD ----------------------------------------------
    print *, 'Reconstruyendo la matriz de datos por SVD'
    rank = 2
    nmr_svd(:,:) = nmr_or(:,:)
    call svd_reconst(M,N,nmr_svd,rank)

    !-- Calculando Precision de la SVD -------------------
    call check_matrix_dif(M,N,nmr_svd,nmr_or,val)
    print *, '   --- Rango: ', rank
    print *, "   --- Precision SVD: ", val

    !-- Metodo ALS ---------------------------------------
    print *, 'Implementando ALS: CUIDADO algoritmo de convergencia incompleto'
    allocate(mat_S(M,rank));allocate(mat_L(rank,N))
    call ALS_factor_analysis(M,N,nmr_svd,ranK,mat_S,mat_L)

    !-- Escribiendo solucion para L ----------------------
    write(sol,*)rank
    sol= 'componentes_R' // trim(adjustl(sol)) // '.res'
    open(33, file=sol)
    do i=1,rank
      write(33,1133) (mat_L(i,j), j=1,N)
    enddo
   1133 format(5000E16.8)
    close(33)

    !-- Escribiendo solucion para S ----------------------
    write(sol,*)rank
    sol= 'concentraciones_R' // trim(adjustl(sol)) // '.res'
    open(34, file=sol)
    do i=1,M
      write(34,1134) (mat_S(i,j), j=1,rank)
    enddo
   1134 format(5000E16.8)
    close(34)

end program nmrFA

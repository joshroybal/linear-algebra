program driver
    use linearalgebra
    implicit none
    ! data
    real, dimension(4, 2) :: A
    real, dimension(2, 4) :: B
    real, dimension(4, 4) :: C, L, U, D
    ! processing
    call rndmat(4, 2, A)
    call rndmat(2, 4, B)
    C = matmul(a, b)
    write (*,*) 'A'
    call printmatrix(4, 2, A)
    write (*,*) 'B'
    call printmatrix(2, 4, B)
    write (*,*) 'C = A * B'
    call printmatrix(4, 4, C)
    call rndmat(4, 4, D)
    write (*,*) 'D'
    call printmatrix(4, 4, D)
    write (*,*) 'Calling LU factorization subroutine.'
    call lufact(4, D, L, U)
    call printmatrix(4, 4, L)
    call printmatrix(4, 4, U)
    write (*,*) '|D| =', det(4, D)
    !write (*,*) 'Calling LUP factorization subroutine.'
    !call lupfact(4, D, L, U)
    !call printmatrix(4, 4, L)
    !call printmatrix(4, 4, U)
    write (*,*) '|D| =', detp(4, D)    
end program driver

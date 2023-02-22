program driver
    use linearalgebra
    implicit none
    ! data
    integer, parameter :: N = 4
    real :: d
    real, dimension(N) :: x, y
    real, dimension(N,N) :: A, AINV, L, U, P
    ! processing
    call rndmat(N, N, A)
    call rndmat(N, N, A)
    call rndmat(N, 1, x)
    !call random_number(A)
    call lufact(N, A, L, U)
    AINV = matinv(N, A)
    d = det(N, A)
    write (*,*) 'no pivoting'
    write (*,*) 'A'
    call printmatrix(N, N, A)
    write (*,*) 'L'
    call printmatrix(N, N, L)
    write (*,*) 'U'
    call printmatrix(N, N, U)
    write (*,*) 'L * U'
    call printmatrix(N, N, matmul(L, U))
    write (*,*) '|A| =', d
    write (*,*) 'AINV'
    call printmatrix(N, N, AINV)
    write (*,*) 'A * AINV'
    call printmatrix(N, N, matmul(A, AINV))
    write (*,*) 'AINV * A'
    call printmatrix(N, N, matmul(AINV, A))
    ! with partial pivoting
    call lupfact(N, A, L, U, P)
    d = determinant(N, A)
    write (*,*) 'partial pivoting'
    write (*,*) 'L'
    call printmatrix(N, N, L)
    write (*,*) 'U'
    call printmatrix(N, N, U)
    write (*,*) 'P'
    call printmatrix(N, N, P)
    write (*,*) 'P * A'
    call printmatrix(N, N, matmul(P, A))
    write (*,*) 'L * U'
    call printmatrix(N, N, matmul(L, U))
    write (*,*) '|A| =', d
    write (*,*) 'system of linear equations'
    write (*,*) x
    y = linsys(N, A, x)
    write (*,*) y
end program driver

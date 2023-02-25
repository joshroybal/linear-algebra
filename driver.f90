program driver
    use linearalgebra
    implicit none
    ! data
    integer, parameter :: N = 4
    double precision :: d
    double precision, dimension(N) :: x, y, err
    double precision, dimension(N,N) :: A, AINV, L, U, P
    ! processing
    A(:n,1) = (/2.,4.,8.,6./)
    A(:n,2) = (/1.,3.,7.,7./)
    A(:n,3) = (/1.,3.,3.,9./)
    A(:n,4) = (/0.,1.,5.,8./)
    !call random_number(A)
    call random_number(x)
    !A = floor(10. * A)
    x = floor(10. * x)
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
    write (*,*) 'partial pivoting'
    call lupfact(N, A, L, U, P)
    d = determinant(N, A)
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
    write (*,*) 'A'
    call printmatrix(N, N, A)
    write (*,*) 'x'
    call printmatrix(N, 1, x)
    y = linsys(N, A, x)
    write (*,*) 'sol''n vector y'
    call printmatrix(N, 1, y)
    write (*,*) 'A * y'
    call printmatrix(N, 1, matmul(A, y))
    err = abs(matmul(A, y) - x)
    write (*,*) 'error vector'
    call printmatrix(N, 1, err)
end program driver

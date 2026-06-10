program driver
    use linearalgebra
    implicit none
    ! data
    integer, parameter :: N = 4
    double precision :: d
    double precision, dimension(N) :: x, y, err
    double precision, dimension(N,N) :: A, AINV, L, U, P
    ! processing
    !write (*,'(a,/)') 'CONTENT-TYPE: US-ASCII'
    A(:n,1) = (/2.,4.,8.,6./)
    A(:n,2) = (/1.,3.,7.,7./)
    A(:n,3) = (/1.,3.,3.,9./)
    A(:n,4) = (/0.,1.,5.,8./)
    !call random_number(A)
    call random_number(x)
    !A = floor(10. * A)
    x = floor(10. * x)
    call lufact(A, L, U)
    AINV = matinv(A)
    d = det(A)
    write (*,*) 'no pivoting'
    write (*,*) 'A'
    call printmatrix(A)
    write (*,*) 'L'
    call printmatrix(L)
    write (*,*) 'U'
    call printmatrix(U)
    write (*,*) 'L * U'
    call printmatrix(matmul(L, U))
    write (*,*) '|A| =', d
    write (*,*) 'AINV'
    call printmatrix(AINV)
    write (*,*) 'A * AINV'
    call printmatrix(matmul(A, AINV))
    write (*,*) 'AINV * A'
    call printmatrix(matmul(AINV, A))
    ! with partial pivoting
    write (*,*) 'partial pivoting'
    call lupfact(A, L, U, P)
    d = determinant(A)
    write (*,*) 'L'
    call printmatrix(L)
    write (*,*) 'U'
    call printmatrix(U)
    write (*,*) 'P'
    call printmatrix(P)
    write (*,*) 'P * A'
    call printmatrix(matmul(P, A))
    write (*,*) 'L * U'
    call printmatrix(matmul(L, U))
    write (*,*) '|A| =', d
    write (*,*) 'system of linear equations'
    write (*,*) 'A'
    call printmatrix(A)
    write (*,*) 'x'
    call printvector(x)
    y = linsys(A, x)
    write (*,*) 'sol''n vector y'
    call printvector(y)
    write (*,*) 'A * y'
    call printvector(matmul(A, y))
    err = abs(matmul(A, y) - x)
    write (*,*) 'error vector'
    call printvector(err)
end program driver

program driver
    use matrix
    implicit none
    ! data
    integer :: i
    double precision, allocatable, dimension(:,:) :: a, ainv, l, u, p
    ! processing
    do i = 2, 5
        allocate(a(i,i))
        allocate(ainv(i,i))
        allocate(l(i,i))
        allocate(u(i,i))
        allocate(p(i,i))
        call random_number(a)
        a = aint(20 * a - 10)
        write (*,*) 'A'
        call dump(a)

        write (*,*) '|A| =', determinant(a)
        call matinv(a, ainv)

        write (*,*) 'INV(A)'
        call dump(ainv)
        write (*,*) 'A * INV(A)'
        call dump(matmul(a, ainv))

        write (*,*) 'LU factorization'
        call lufact(a, l, u)
        write (*,*) 'L'
        call dump(l)
        write (*,*) 'U'
        call dump(u)
        write (*,*) 'error'
        call dump(matmul(u, l) - a)

        write (*,*) 'LUP factorization'
        write (*,*) 'A'
        call dump(a)
        call lupfact(a, l, u, p)
        write (*,*) 'L'
        call dump(l)
        write (*,*) 'U'
        call dump(u)
        write (*,*) 'P'
        call dump(p)
        write (*,*) 'PA'
        call dump(matmul(a, p))
        write (*,*) 'LU'
        call dump(matmul(u, l))
        write (*,*) 'error'
        call dump(matmul(matmul(u, l), transpose(p)) - a)
        deallocate(a)
        deallocate(ainv)
        deallocate(l)
        deallocate(u)
        deallocate(p)
    end do
end program

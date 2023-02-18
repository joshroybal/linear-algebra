module linearalgebra
    use prng
    implicit none
    contains
        ! print matrix by rows
        subroutine printmatrix(m, n, A)
            ! dummy arguments
            integer, intent(in) :: m, n
            real, dimension(m,n), intent(in) :: A
            ! local data
            integer :: i, j
            ! processing
            do i = 1, m
                j = 0
                !write (*,'(3f16.3)') A(i,1:n)
                write (*,*) (A(i,j),j=1,n)
            end do
            write (*,*)
        end subroutine printmatrix
        
        ! return identity matrix
        function ident(n) result(I)
            ! dummy arguments
            integer, intent(in) :: n
            ! function result location
            real, dimension(n,n) :: I
            ! local data
            integer :: j, k
            ! processing
            do k = 1, n
                do j = 1, n
                    if (j .eq. k) then
                        I(j,k) = 1
                    else
                        I(j,k) = 0
                    end if
                end do
            end do
        end function ident

        ! invariant pseudo-random sequence
        subroutine rndmat(m, n, A)
            ! dummy arguments
            integer, intent(in) :: m, n
            real, dimension(m, n), intent(out) :: A
            ! local data
            integer :: i, j
            ! processing
            call setseed(1000)
            do j = 1, n
                do i = 1, m
                    A(i,j) = int(10.0 * randomreal()) - 5.0
                end do
            end do
        end subroutine rndmat

        ! Exchange real locations
        subroutine exchange(r1, r2)
            ! dummy arguments
            real, intent(inout) :: r1, r2
            ! local data
            real :: tmp
            ! processing
            tmp = r1
            r1 = r2
            r2 = tmp
        end subroutine exchange

        ! Exchange rows of matrix
        subroutine interchange(m, n, i, j, A)
            ! dummy arguments
            integer, intent(in) :: m, n, i, j
            real, dimension(m,n), intent(inout) :: A
            ! local data
            integer :: k
            ! processing
            do k = 1, n
                call exchange(A(i,k), A(j,k))
            end do
        end subroutine interchange

        ! multiply ith row by scalar s
        subroutine scalerow(m, n, i, s, A)
            ! dummy arguments
            integer, intent(in) :: m, n, i
            real, intent(in) :: s
            real, dimension(m,n), intent(inout) :: A
            ! processing
            A(i,1:n) = s * A(i,1:n)
        end subroutine scalerow

        ! Add row i to row j of matrix
        subroutine addrows(m, n, i, j, A)
            ! dummy arguments
            integer, intent(in) :: m, n, i, j
            real, dimension(m,n), intent(inout) :: A
            ! processing
            A(j,1:n) = A(j,1:n) + A(i,1:n)
        end subroutine addrows
        
        ! LU decomposition
        subroutine lufact(n, A, L, U)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n, n), intent(in) :: A
            real, dimension(n, n), intent(out) :: L, U
            ! local data
            integer :: i, j
            real :: multiplier
            ! processing
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n) 
                end do
            end do
        end subroutine lufact
        
        ! LUP factorization - result is no. of pivoting interchanges
        function lupfact(n, A, L, U, P) result(s)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n, n), intent(in) :: A
            real, dimension(n, n), intent(out) :: L, U, P
            ! function result location
            integer :: s
            ! local data
            integer :: i, j
            real :: multiplier
            ! processing
            s = 0
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                s = s + pivotp(n, n, i, U, P)
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n) 
                end do
            end do
            print *, 'L'
            call printmatrix(n, n, L)
            print *, 'U'
            call printmatrix(n, n, U)
            print *, 'P'
            call printmatrix(n, n, P)
        end function lupfact

        ! partial pivoting - row interchange only
        subroutine pivot(m, n, col, A)
            ! dummy arguments
            integer, intent(in) :: m, n, col
            real, dimension(m,n), intent(inout) :: A
            ! local data
            integer :: i, j, maxidx
            ! processing
            do i = col, n - 1
                maxidx = i
                do j = i + 1, n
                    if (abs(A(j,col)) > abs(A(maxidx,col))) maxidx = j
                end do
                if (maxidx .ne. i) call interchange(n, n, i, maxidx, A)
            end do
        end subroutine pivot

        ! partial pivoting with permutation matrix
        ! result is no. of substitutions along the way
        function pivotp(m, n, col, A, P) result(s)
            ! dummy arguments
            integer, intent(in) :: m, n, col
            real, dimension(m,n), intent(inout) :: A, P
            ! function result location
            integer :: s
            ! local data
            integer :: i, j, maxidx
            ! processing
            s = 0
            do i = col, n - 1
                maxidx = i
                do j = i + 1, n
                    if (abs(A(j,col)) > abs(A(maxidx,col))) maxidx = j
                end do
                if (maxidx .ne. i) then
                    call interchange(n, n, i, maxidx, A)
                    call interchange(n, n, i, maxidx, P)
                    s = s + 1
                end if
            end do
        end function pivotp
        
        ! dot product
        function dot(n, A, B) result(r)
            ! dummy argument
            integer, intent(in) :: n
            real, dimension(n), intent(in) :: A, B
            ! function result location
            real :: r
            ! processing
            r = sum(A*B)
        end function dot

        ! determinant
        function det(n, A) result(r)
            ! dummy argument
            integer, intent(in) :: n
            real, dimension(n), intent(in) :: A
            ! function result location
            real :: r
            ! local data
            integer :: i
            real, dimension(n,n) :: L, U
            ! processing
            call lufact(n, A, L, U)
            r = 1.0
            do i = 1, n
                r = r * U(i,i)
            end do
        end function det

        ! determinant using partial pivoting
        function detp(n, A) result(r)
            ! dummy argument
            integer, intent(in) :: n
            real, dimension(n), intent(in) :: A
            ! function result location
            real :: r
            ! local data
            integer :: i, s
            real, dimension(n,n) :: L, U, P
            ! processing
            s = lupfact(n, A, L, U, P)
            r = 1.0
            do i = 1, n
                r = r * U(i,i)
            end do
            r = -1**s * r
        end function detp        
end module linearalgebra

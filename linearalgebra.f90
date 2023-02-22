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
            integer :: i
            !character (len=12) :: nstr, fmtstr
            ! processing
            !write (nstr,*) n
            !nstr = adjustl(nstr)
            !fmtstr = '('//trim(nstr)//'f15.9)'
            do i = 1, m
                !write (*,fmtstr) (A(i,j),j=1,n)
                write (*,*) A(i,1:n)
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
            !call setseed(1000)
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

        ! Row reduce matrix to row echelon form.
        subroutine reducerows(n, A)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n,n), intent(inout) :: A
            ! local data
            integer :: i, j
            real :: multiplier
            ! processing
            do i = 1, n - 1
                call pivot(n, n, i, A)
                do j = i + 1, n
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                end do
            end do
            do i = n, 2, -1
                do j = i - 1, 1, -1
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                end do
            end do
            do i = 1, n
                multiplier = 1./A(i,i)
                A(i,1:n) = multiplier*A(i,1:n)
            end do
        end subroutine reducerows

        ! Row reduce matrix to row echelon form.
        ! Also solves systems of linear equations.
        subroutine rowechelon(m, n, A)
            ! dummy arguments
            integer, intent(in) :: m, n
            real, dimension(m,n), intent(inout) :: A
            ! local data
            integer :: i, j
            real :: multiplier
            ! processing
            do i = 1, m - 1
                !call pivot(m, n, i, A)
                do j = i + 1, m
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                end do
            end do
            do i = m, 2, -1
                do j = i - 1, 1, -1
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                end do
            end do
            do i = 1, m
                multiplier = 1./A(i,i)
                A(i,1:n) = multiplier*A(i,1:n)
            end do
        end subroutine rowechelon

        ! Row reduce matrix A to row echelon form.
        ! Performs like operations on matrix B.
        subroutine reducerows2(n, A, B)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n,n), intent(inout) :: A, B
            ! local data
            double precision, parameter :: epsilon = tiny(1.)
            integer :: i, j
            real :: multiplier
            ! processing
            do i = 1, n - 1
                if (abs(A(i,i)) .le. epsilon) cycle
                call pivot2(n, n, i, A, B)
                do j = i + 1, n
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                    B(j,1:n) = multiplier*B(i,1:n) + B(j,1:n)
                end do
            end do
            do i = n, 2, -1
                if (abs(A(i,i)) .le. epsilon) cycle
                do j = i - 1, 1, -1
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                    B(j,1:n) = multiplier*B(i,1:n) + B(j,1:n)
                end do
            end do
            do i = 1, n
                if (abs(A(i,i)) .le. epsilon) cycle
                multiplier = 1./A(i,i)
                A(i,1:n) = multiplier*A(i,1:n)
                B(i,1:n) = multiplier*B(i,1:n)
            end do
        end subroutine reducerows2

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

        subroutine lupfact(n, A, L, U, P)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n,n), intent(in) :: A
            real, dimension(n,n), intent(out) :: L, U, P
            !local data
            integer :: i, j
            real :: multiplier
            ! processing
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                call pivot2(n, n, i, U, P)
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n)
                end do
            end do
        end subroutine lupfact

        ! LUP factorization and partial pivoting determinant
        function determinant(n, A) result(det)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n, n), intent(in) :: A
            ! function result location
            real :: det
            ! local data
            integer :: i, j
            real :: multiplier, s
            real, dimension(n) :: diagonal
            real, dimension(n,n) :: L, U, P
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
            diagonal = (/ (U(i,i),i=1,n) /)
            det = (-1)**s * product(diagonal)
        end function determinant

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

        ! partial pivoting - row interchange only - move 2nd matrix in parallel
        subroutine pivot2(m, n, col, A, B)
            ! dummy arguments
            integer, intent(in) :: m, n, col
            real, dimension(m,n), intent(inout) :: A, B
            ! local data
            integer :: i, j, maxidx
            ! processing
            do i = col, n - 1
                maxidx = i
                do j = i + 1, n
                    if (abs(A(j,col)) > abs(A(maxidx,col))) maxidx = j
                end do
                if (maxidx .ne. i) then
                    call interchange(n, n, i, maxidx, A)
                    call interchange(n, n, i, maxidx, B)
                end if
            end do
        end subroutine pivot2

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
            r = product( (/ (U(i,i),i=1,n) /) )
        end function det

        ! Matrix inversion, Gaussian elimination with partial pivoting.
        function matinv(n, A) result(AINV)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n,n), intent(in) :: A
            ! function result location
            real, dimension(n,n) :: AINV
            ! local data
            real, dimension(n,n) :: ACOPY
            ! processing
            ACOPY = A
            AINV = ident(n)
            call reducerows2(n, ACOPY, AINV)
        end function matinv

        ! Solve system of linear equations
        function linsys(n, A, x) result(y)
            ! dummy arguments
            integer, intent(in) :: n
            real, dimension(n,n), intent(in) :: A
            real, dimension(n), intent(in) :: x
            ! function result location
            real, dimension(n) :: y
            ! local data
            real, dimension(n,n+1) :: T
            ! processing
            T(1:n,1:n) = A
            T(1:n,n+1) = x
            call rowechelon(n, n+1, T)
            y = T(1:n,n+1)
        end function linsys
end module linearalgebra

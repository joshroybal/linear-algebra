module linearalgebra
    implicit none
    contains
        ! print matrix by rows
        subroutine printmatrix(m, n, A)
            ! dummy arguments
            integer, intent(in) :: m, n
            double precision, dimension(m,n), intent(in) :: A
            ! local data
            integer :: i
            do i = 1, m
                write (*,*) A(i,1:n)
            end do
            write (*,*)
        end subroutine printmatrix

        ! return identity matrix
        function ident(n) result(I)
            ! dummy arguments
            integer, intent(in) :: n
            ! function result location
            double precision, dimension(n,n) :: I
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

        ! Exchange double precision locations
        subroutine exchange(r1, r2)
            ! dummy arguments
            double precision, intent(inout) :: r1, r2
            ! local data
            double precision :: tmp
            ! processing
            tmp = r1
            r1 = r2
            r2 = tmp
        end subroutine exchange

        ! Exchange rows of matrix
        subroutine interchange(m, n, i, j, A)
            ! dummy arguments
            integer, intent(in) :: m, n, i, j
            double precision, dimension(m,n), intent(inout) :: A
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
            double precision, intent(in) :: s
            double precision, dimension(m,n), intent(inout) :: A
            ! processing
            A(i,1:n) = s * A(i,1:n)
        end subroutine scalerow

        ! Add row i to row j of matrix
        subroutine addrows(m, n, i, j, A)
            ! dummy arguments
            integer, intent(in) :: m, n, i, j
            double precision, dimension(m,n), intent(inout) :: A
            ! processing
            A(j,1:n) = A(j,1:n) + A(i,1:n)
        end subroutine addrows

        ! Row reduce matrix to row echelon form.
        ! Also solves systems of linear equations.
        subroutine reducerows(m, n, A)
            ! dummy arguments
            integer, intent(in) :: m, n
            double precision, dimension(m,n), intent(inout) :: A
            ! local data
            integer :: i, j
            double precision :: multiplier
            ! processing
            do i = 1, m - 1
                call pivot(m, n, i, A)
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
        end subroutine reducerows

        ! Row reduce matrix A to row echelon form.
        ! Performs like operations on matrix B.
        subroutine reducerows2(n, A, B)
            ! dummy arguments
            integer, intent(in) :: n
            double precision, dimension(n,n), intent(inout) :: A, B
            ! local data
            integer :: i, j
            double precision :: multiplier
            ! processing
            do i = 1, n - 1
                call pivot2(n, i, A, B)
                do j = i + 1, n
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                    B(j,1:n) = multiplier*B(i,1:n) + B(j,1:n)
                end do
            end do
            do i = n, 2, -1
                do j = i - 1, 1, -1
                    multiplier = -A(j,i)/A(i,i)
                    A(j,1:n) = multiplier*A(i,1:n) + A(j,1:n)
                    B(j,1:n) = multiplier*B(i,1:n) + B(j,1:n)
                end do
            end do
            do i = 1, n
                multiplier = 1./A(i,i)
                A(i,1:n) = multiplier*A(i,1:n)
                B(i,1:n) = multiplier*B(i,1:n)
            end do
        end subroutine reducerows2

        ! LU decomposition
        subroutine lufact(n, A, L, U)
            ! dummy arguments
            integer, intent(in) :: n
            double precision, dimension(n, n), intent(in) :: A
            double precision, dimension(n, n), intent(out) :: L, U
            ! local data
            integer :: i, j
            double precision :: multiplier
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
            double precision, dimension(n,n), intent(in) :: A
            double precision, dimension(n,n), intent(out) :: L, U, P
            !local data
            integer :: i, j
            double precision :: multiplier
            ! processing
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                call lupivot(n, i, L, U, P)
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n)
                end do
            end do
        end subroutine lupfact

        ! LUP decomposition auxiliary pivoting subroutine
        subroutine lupivot(n, col, L, U, P)
            ! dummy arguments
            integer, intent(in) :: n, col
            double precision, dimension(n,n), intent(inout) :: L, U, P
            ! local data
            integer :: i, j, maxidx
            double precision :: tmp
            ! processing
            maxidx = col
            do i = col + 1, n
                if (abs(U(i,col)) > abs(U(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(n, n, col, maxidx, U)
                call interchange(n, n, col, maxidx, P)
                do j = 1, col - 1
                    tmp = L(col,j)
                    L(col,j) = L(maxidx,j)
                    L(maxidx,j) = tmp
                end do
            end if
        end subroutine lupivot

        ! LUP factorization and partial pivoting determinant
        function determinant(n, A) result(det)
            ! dummy arguments
            integer, intent(in) :: n
            double precision, dimension(n, n), intent(in) :: A
            ! function result location
            double precision :: det
            ! local data
            integer :: i, j
            double precision :: multiplier, s
            double precision, dimension(n) :: diagonal
            double precision, dimension(n,n) :: L, U, P
            ! processing
            s = 0
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                if (pivotp(n, i, U, P) .eqv. .true.) s = s + 1
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
            double precision, dimension(m,n), intent(inout) :: A
            ! local data
            integer :: i, maxidx
            ! processing
            maxidx = col
            do i = col + 1, m
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) call interchange(m, n, col, maxidx, A)
        end subroutine pivot

        ! partial pivoting - row interchange only - move 2nd matrix in parallel
        subroutine pivot2(n, col, A, B)
            ! dummy arguments
            integer, intent(in) :: n, col
            double precision, dimension(n,n), intent(inout) :: A, B
            ! local data
            integer :: i, maxidx
            ! processing
            maxidx = col
            do i = col + 1, n
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(n, n, col, maxidx, A)
                call interchange(n, n, col, maxidx, B)
            end if
        end subroutine pivot2

        ! partial pivoting with permutation matrix
        ! result is no. of substitutions along the way
        function pivotp(n, col, A, P) result(q)
            ! dummy arguments
            integer, intent(in) :: n, col
            double precision, dimension(n,n), intent(inout) :: A, P
            ! function result location
            logical :: q
            ! local data
            integer :: i, maxidx
            ! processing
            q = .false.
            maxidx = col
            do i = col + 1, n
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(n, n, col, maxidx, A)
                call interchange(n, n, col, maxidx, P)
                q = .true.
            end if
        end function pivotp

        ! determinant no pivoting
        function det(n, A) result(r)
            ! dummy argument
            integer, intent(in) :: n
            double precision, dimension(n), intent(in) :: A
            ! function result location
            double precision :: r
            ! local data
            integer :: i
            double precision, dimension(n,n) :: L, U
            ! processing
            call lufact(n, A, L, U)
            r = product( (/ (U(i,i),i=1,n) /) )
        end function det

        ! Matrix inversion, Gaussian elimination with partial pivoting.
        function matinv(n, A) result(AINV)
            ! dummy arguments
            integer, intent(in) :: n
            double precision, dimension(n,n), intent(in) :: A
            ! function result location
            double precision, dimension(n,n) :: AINV
            ! local data
            double precision, dimension(n,n) :: ACOPY
            ! processing
            ACOPY = A
            AINV = ident(n)
            call reducerows2(n, ACOPY, AINV)
        end function matinv

        ! Solve system of linear equations
        function linsys(n, A, x) result(y)
            ! dummy arguments
            integer, intent(in) :: n
            double precision, dimension(n,n), intent(in) :: A
            double precision, dimension(n), intent(in) :: x
            ! function result location
            double precision, dimension(n) :: y
            ! local data
            double precision, dimension(n,n+1) :: T
            ! processing
            T(1:n,1:n) = A
            T(1:n,n+1) = x
            call reducerows(n, n+1, T)
            y = T(1:n,n+1)
        end function linsys
end module linearalgebra

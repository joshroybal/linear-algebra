module linearalgebra
    implicit none
    contains
        ! princ column vector
        subroutine printvector(v)
            ! dummy argument
            double precision, dimension(:), intent(in) :: v
            ! local data
            integer :: i
            ! processing
            do i = 1, size(v)
                write (*,*) v(i)
            end do
            write (*,*)
        end subroutine printvector

        ! print matrix by rows
        subroutine printmatrix(A)
            ! dummy argument
            double precision, dimension(:,:), intent(in) :: A
            ! local data
            integer :: i
            ! processing
            do i = 1, size(A, 1)
                write (*,*) A(i,1:size(A,2))
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
        subroutine interchange(i, j, A)
            ! dummy arguments
            integer, intent(in) :: i, j
            double precision, dimension(:,:), intent(inout) :: A
            ! local data
            integer :: k
            ! processing
            do k = 1, size(A, 2)
                call exchange(A(i,k), A(j,k))
            end do
        end subroutine interchange

        ! Row reduce matrix to row echelon form.
        ! Also solves systems of linear equations.
        subroutine reducerows(A)
            ! dummy arguments
            double precision, dimension(:,:), intent(inout) :: A
            ! local data
            integer :: i, j, m, n
            double precision :: multiplier
            ! processing
            m = size(A, 1) ; n = size(A, 2)
            do i = 1, m - 1
                call pivot(i, A)
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
        subroutine reducerows2(A, B)
            ! dummy arguments
            double precision, dimension(:,:), intent(inout) :: A, B
            ! local data
            integer :: i, j, n
            double precision :: multiplier
            ! processing
            n = size(A, 1)
            do i = 1, n - 1
                call pivot2(i, A, B)
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
        subroutine lufact(A, L, U)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: A
            double precision, dimension(:,:), intent(out) :: L, U
            ! local data
            integer :: i, j, n
            double precision :: multiplier
            ! processing
            n = size(A, 1)
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

        subroutine lupfact(A, L, U, P)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: A
            double precision, dimension(:,:), intent(out) :: L, U, P
            !local data
            integer :: i, j, n
            double precision :: multiplier
            ! processing
            n = size(A, 1)
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                call lupivot(i, L, U, P)
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n)
                end do
            end do
        end subroutine lupfact

        ! LUP decomposition auxiliary pivoting subroutine
        subroutine lupivot(col, L, U, P)
            ! dummy arguments
            integer, intent(in) :: col
            double precision, dimension(:,:), intent(inout) :: L, U, P
            ! local data
            integer :: i, j, maxidx
            double precision :: tmp
            ! processing
            maxidx = col
            do i = col + 1, size(U, 1)
                if (abs(U(i,col)) > abs(U(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(col, maxidx, U)
                call interchange(col, maxidx, P)
                do j = 1, col - 1
                    tmp = L(col,j)
                    L(col,j) = L(maxidx,j)
                    L(maxidx,j) = tmp
                end do
            end if
        end subroutine lupivot

        ! LUP factorization and partial pivoting determinant
        function determinant(A) result(d)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: A
            ! function result location
            double precision :: d
            ! local data
            integer :: i, j, n
            double precision :: multiplier, s
            double precision, dimension(size(A, 1)) :: diagonal
            double precision, dimension(size(A, 1),size(A, 2)) :: L, U, P
            ! processing
            n = size(diagonal)
            s = 0
            P = ident(n)    ! permuation matrix
            L = ident(n)    ! lower triangular matrix
            U = A           ! upper triangular matrix
            do i = 1, n - 1
                if (pivotp(i, U, P) .eqv. .true.) s = s + 1
                do j = i + 1, n
                    multiplier = -U(j,i)/U(i,i)
                    L(j,i) = -multiplier
                    U(j,1:n) = multiplier*U(i,1:n) + U(j,1:n)
                end do
            end do
            diagonal = (/ (U(i,i),i=1,n) /)
            d = (-1)**s * product(diagonal)
        end function determinant

        ! partial pivoting - row interchange only
        subroutine pivot(col, A)
            ! dummy arguments
            integer, intent(in) :: col
            double precision, dimension(:,:), intent(inout) :: A
            ! local data
            integer :: i, maxidx
            ! processing
            maxidx = col
            do i = col + 1, size(A, 1)
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) call interchange(col, maxidx, A)
        end subroutine pivot

        ! partial pivoting - row interchange only - move 2nd matrix in parallel
        subroutine pivot2(col, A, B)
            ! dummy arguments
            integer, intent(in) :: col
            double precision, dimension(:,:), intent(inout) :: A, B
            ! local data
            integer :: i, maxidx
            ! processing
            maxidx = col
            do i = col + 1, size(A, 1)
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(col, maxidx, A)
                call interchange(col, maxidx, B)
            end if
        end subroutine pivot2

        ! partial pivoting with permutation matrix
        ! result is no. of substitutions along the way
        function pivotp(col, A, P) result(q)
            ! dummy arguments
            integer, intent(in) :: col
            double precision, dimension(:,:), intent(inout) :: A, P
            ! function result location
            logical :: q
            ! local data
            integer :: i, maxidx
            ! processing
            q = .false.
            maxidx = col
            do i = col + 1, size(A, 1)
                if (abs(A(i,col)) > abs(A(maxidx,col))) maxidx = i
            end do
            if (maxidx .ne. col) then
                call interchange(col, maxidx, A)
                call interchange(col, maxidx, P)
                q = .true.
            end if
        end function pivotp

        ! determinant no pivoting
        function det(A) result(r)
            ! dummy argument
            double precision, dimension(:,:), intent(in) :: A
            ! function result location
            double precision :: r
            ! local data
            integer :: i
            double precision, dimension(size(A,1),size(A,2)) :: L, U
            ! processing
            call lufact(A, L, U)
            r = product( (/ (U(i,i),i=1,size(A,1)) /) )
        end function det

        ! Matrix inversion, Gaussian elimination with partial pivoting.
        function matinv(A) result(AINV)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: A
            ! function result location
            double precision, dimension(size(A,1),size(A,2)) :: AINV
            ! local data
            double precision, dimension(size(A,1),size(A,2)) :: ACOPY
            ! processing
            ACOPY = A
            AINV = ident(size(A,1))
            call reducerows2(ACOPY, AINV)
        end function matinv

        ! Solve system of linear equations
        function linsys(A, x) result(y)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: A
            double precision, dimension(:), intent(in) :: x
            ! function result location
            double precision, dimension(size(x)) :: y
            ! local data
            integer :: n
            double precision, dimension(size(A,1),size(A,2)+1) :: T
            ! processing
            n = size(A,2)
            print *, n
            T(1:n,1:n) = A
            T(1:n,n+1) = x
            call reducerows(T)
            y = T(1:n,n+1)
        end function linsys
end module linearalgebra

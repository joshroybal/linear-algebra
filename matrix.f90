module matrix
    implicit none
    contains
        ! show matrix, perhaps the trickiest subprogam of the module!
        subroutine dump(matrix)
            ! dummy argment
            double precision, dimension(:,:), intent(in) :: matrix
            ! local data
            integer :: j
            ! processing
            do j = 1, size(matrix, 2)
                write (*,*) matrix(:,j)
            end do
            write (*,*)
        end subroutine

        ! return identity matrix as result
        function ident(n) result(i)
            ! dummy argument
            integer, intent(in) :: n
            ! function result location
            integer, dimension(n,n) :: i
            ! local data
            integer :: j
            ! processing
            i = 0
            do j = 1, n
                i(j,j) = 1
            end do
        end function

        ! subprogram to exchange rows (columns actually)
        subroutine exchange(u, v)
            ! dummy arguments
            double precision, dimension(:) :: u, v
            ! local data
            double precision, dimension(size(u)) :: t
            ! processing
            t = u ; u = v ; v = t
        end subroutine

        ! it's been demonstrated to myself that contiguity is leveraged
        ! compute determinant of matrix - partial pivoting
        function determinant(matrix) result(det)
            ! dummy argument
            double precision, dimension(:,:), intent(in) :: matrix
            ! result location
            double precision :: det
            ! local data
            integer :: j, k, n, pividx, s
            double precision :: factor
            double precision, dimension(size(matrix, 1), size(matrix, 2)) :: upper
            ! processsing
            upper = matrix
            n = size(matrix, 1)
            s = 1
            do j = 1, n
                pividx = j - 1 + maxloc(abs(upper(j,j:)), dim=1)
                if (pividx .ne. j) then
                    call exchange(upper(j:,j), upper(j:,pividx))
                    s = -1 * s
                end if
                if (abs(upper(j,j)) .le. 1e-5) then
                    det = 0
                    return
                end if
                do k = j + 1, n
                    factor = upper(j,k) / upper(j,j)
                    upper(j:,k) = upper(j:,k) - factor * upper(j:,j)
                end do
            end do
            det = s * product((/ (upper(j,j), j = 1, n) /))
        end function

        ! find matrix inverse by Gaussuan elimination
        subroutine matinv(a, ainv)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: a
            double precision, dimension(:,:), intent(out) :: ainv
            ! local data
            integer :: i, k, n, pivot
            double precision :: factor
            double precision, dimension(2*size(a,1), size(a,2)) :: work
            ! processing
            n = size(a, 1)
            work(1:n,1:n) = a
            work(n+1:n+n,1:n) = ident(n)

            ! forwards
            do k = 1, n
                pivot = k - 1 + maxloc(abs(work(k,k:)), dim=1)
                if (pivot .ne. k) call exchange(work(k:,k), work(k:,pivot))
                if (abs(work(k,k)) .le. 1e-5) then
                    write (*,*) 'singular matrix'
                    ainv = 0
                    return
                end if
                do i = k + 1, n
                    factor = work(k,i) / work(k,k)
                    work(k:,i) = work(k:,i) - factor * work(k:,k)   
                end do
            end do
            
            ! backwards
            do k = n, 1, -1
                work(k:,k) = work(k:,k) / work(k,k)
                do i = k - 1, 1, -1
                    factor = work(k,i)
                    work(k:,i) = work(k:,i) - factor * work(k:,k)
                end do
            end do
            ainv = work(n+1:n+n,1:n)
        end subroutine

        ! lu factorization - square matrices
        subroutine lufact(matrix, lower, upper)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: matrix
            double precision, dimension(:,:), intent(out) :: lower, upper
            ! local data
            integer :: j, k, n
            ! processing
            n = size(matrix, 1)
            lower = ident(n)
            upper = matrix
            do j = 1, n - 1
                if (abs(upper(j,j)) .le. 1e-5) then
                    write (*,*) 'simple lu factorization not possible'
                    return
                end if
                do k = j + 1, n
                    lower(j,k) = upper(j,k) / upper(j,j)
                    upper(j:,k) = upper(j:,k) - lower(j,k) * upper(j:,j)
                end do
            end do
        end subroutine

        ! LUP factorization of a square mattrix
        subroutine lupfact(a, l, u, p)
            ! dummy arguments
            double precision, dimension(:,:), intent(in) :: a
            double precision, dimension(:,:), intent(out) :: l, u, p
            ! local data
            integer :: j, k, n, pividx
            ! processing
            n = size(a, 1)
            l = ident(n)
            u = a
            p = ident(n)
            do j = 1, n
                pividx = j - 1 + maxloc(abs(u(j,j:)), dim=1)
                if (pividx .ne. j) then
                    call exchange(u(j:,j), u(j:,pividx))
                    call exchange(p(:,j), p(:,pividx))
                    if (j .gt. 1) call exchange(l(1:j-1,j), l(1:j-1,pividx))
                end if
                if (abs(u(j,j)) .le. 1e-5) then
                    u(j,j:) = 0
                    cycle
                end if
                do k = j + 1, n
                    l(j,k) = u(j,k) / u(j,j)
                    u(j:,k) = u(j:,k) - l(j,k) * u(j:,j)
                end do
            end do
        end subroutine
end module

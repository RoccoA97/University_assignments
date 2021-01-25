! 06.f90
! gfortran -Wall -ffree-line-length-0 -o 06.o 06.f90 -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack





! module for shrodinger equation solving utilities
module sch_eq_utils
    implicit none

    contains
        ! initialize hamiltonian operator
        function init_ham(N, L, m, omega) result(ham)
            ! input arguments
            integer(4) :: N     ! hamiltonian dimesnions
            real(8)    :: L     ! width of the box
            real(8)    :: m     ! mass
            real(8)    :: omega ! omega

            real(8), dimension(N+1,N+1) :: ham
            real(8) :: DL
            integer(4) :: ii, jj

            DL = L / N
            do ii=1,N+1
                do jj=1,N+1
                    if (ABS(ii-jj)==1) then
                        ham(ii,jj) = - (1.0d0/(2.0d0*m)) * (1.0d0/(DL**2))
                    else if (ii==jj) then
                        ham(ii,jj) =   (1.0d0/(2.0d0*m)) * (2.0d0/(DL**2)) + (0.5d0*m*omega**2) * (-L/2 + DL*(ii-1))**2
                    else
                        ham(ii,jj) =    0.0d0
                    end if
                end do
            end do
        end function init_ham


        ! diagonalize hamiltonian operator
        subroutine diag_ham(ham, eigs)
            ! input arguments
            real(8), dimension(:,:) :: ham
            real(8), dimension(:)   :: eigs

            real(8), dimension(size(ham,1))  :: dd
            real(8), dimension(size(dd,1)-1) :: sd

            real(8), allocatable :: work(:)
            integer(4) :: N, lwork, info
            integer(4) :: ii

            N = size(ham, 1)

            ! fill diagonal
            do ii=1,N
                dd(ii) = ham(ii,ii)
            end do
            ! fill subdiagonal
            do ii=1,N-1
                sd(ii) = ham(ii,ii+1)
            end do

            N = size(dd, 1)
            lwork = max(1,2*N-2)

            allocate(work(lwork))
            call dstev('V', N, dd, sd, ham, N, work, info)
            deallocate(work)

            eigs = dd
        end subroutine diag_ham


        ! print eigenfunctions on file
        subroutine print_eigfc(xs, ys, filename, unit)
            ! input arguments
            real(8), dimension(:) :: xs
            real(8), dimension(:) :: ys
            character(*) :: filename
            integer(4) :: unit

            integer(4) :: ii

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(xs,1)
                write(unit,'(g0ag0)') xs(ii), char(9), ys(ii)
            end do
            close(unit)
        end subroutine print_eigfc


        ! print eigenvalues on file
        subroutine print_eigvl(xs, es, eigs, filename, unit)
            ! input arguments
            real(8), dimension(:) :: xs
            real(8), dimension(:) :: es
            real(8), dimension(:) :: eigs
            character(*) :: filename
            integer(4) :: unit

            integer(4) :: ii

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(xs,1)
                write(unit,'(g0ag0ag0ag0)') xs(ii), char(9), es(ii), char(9), eigs(ii), char(9), eigs(ii) - es(ii)
            end do
            close(unit)
        end subroutine print_eigvl
end module sch_eq_utils





program main
    use sch_eq_utils
    implicit none

    real(8), dimension(:,:), allocatable :: ham
    real(8), dimension(:), allocatable   :: eigs
    real(8), dimension(:), allocatable   :: xs
    real(8), dimension(:), allocatable   :: ns
    real(8), dimension(:), allocatable   :: es

    integer(4)   :: arg1, arg3
    real(8)      :: arg2, arg4, arg5
    real(8)      :: L, DL, M, O
    integer(4)   :: N, K
    integer(4)   :: ii, jj
    character(4) :: col

    character(100) :: arg1char
    character(100) :: arg2char
    character(100) :: arg3char
    character(100) :: arg4char
    character(100) :: arg5char


    ! get command-line arguments
    CALL GET_COMMAND_ARGUMENT(1,arg1char)   ! 1 : number of discretization intervals
    CALL GET_COMMAND_ARGUMENT(2,arg2char)   ! 2 : size of space interval [-a,a]
    CALL GET_COMMAND_ARGUMENT(3,arg3char)   ! 3 : number of eigenfunctions to store
    CALL GET_COMMAND_ARGUMENT(4,arg4char)   ! 4 : mass
    CALL GET_COMMAND_ARGUMENT(5,arg5char)   ! 5 : omega

    READ(arg1char, *) arg1
    READ(arg2char, *) arg2
    READ(arg3char, *) arg3
    READ(arg4char, *) arg4
    READ(arg5char, *) arg5

    ! rename command-line arguments for a better meaning
    N = arg1
    L = arg2
    K = arg3
    M = arg4
    O = arg5

    DL = L / N

    allocate(ham(N+1,N+1))
    allocate(eigs(N+1))
    allocate(xs(N+1))
    allocate(ns(N+1))
    allocate(es(N+1))


    ! initialize discretized Hamiltonian with
    ! * N discretization intervals
    ! * L length of space interval [-a,a]
    ! * M mass and O omega parameters
    ham  = init_ham(N, L, M, O)

    call diag_ham(ham, eigs)

    do ii=1,size(xs,1)
        xs(ii) = (-L/2) + (ii-1)*DL
    end do

    ! normalize eigenfunctions
    do ii=1,size(ham,1)
        do jj=1,size(ham,2)
            ham(ii,jj) = ham(ii,jj) * SQRT(N/L)
        end do
    end do

    do ii=1,K
        if (ii < 10) then
            write(col, "(i1)") ii
        else if (ii < 100) then
            write(col, "(i2)") ii
        else
            write(col, "(i3)") ii
        end if
        call print_eigfc(xs, ham(:,ii), "results/"//TRIM(col)//"_eigenfunction.dat", 1)
    end do


    do ii=1,size(ns,1)
        ns(ii) =  (ii - 1.0d0)
        es(ii) = ((ii - 1.0d0) + 0.5d0) * O
    end do

    call print_eigvl(ns, es, eigs, "results/eigenvalues.dat", 1)
end program main

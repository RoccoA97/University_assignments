! 09.f90
! gfortran -Wall -ffree-line-length-0 -o 09.o 09.f90 -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack





! module err_handling_mod for error handling
module err_handling
    implicit none

    contains
        ! checkpoint subroutine
        subroutine check(debug, condition, msg_type, msg, trg_stop, var)
            ! mandatory args
            logical, intent(in) :: debug
            logical, intent(in) :: condition
            character(*) :: msg_type
            character(*) :: msg
            ! optional args
            logical, intent(in), optional :: trg_stop
            class(*), optional :: var
            ! dummy arguments
            logical :: trg_stop_


            ! assign defaul value to trg_stop argument gfortran -Wall -ffree-line-length-0 -o 09.o 09.f90 -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapackthrough dummy variable
            if (.NOT.present(trg_stop)) then
                trg_stop_ = .TRUE.
            else
                trg_stop_ = trg_stop
            end if

            ! check if error/warning/... condition is satisfied
            if (condition) then

                ! check if debug is enabled
                if (debug) then

                    ! choose among message type (error/warning) if given
                    ! add color and bold style depending on the type of message
                    ! Error -> bold red ([1;31m)
                    ! Warning -> bold pink ([1;95m)
                    ! Other -> bold green ([1;32m)
                    if (msg_type == "Error") then
                        print *, "["//achar(27)//"[1;31m"//msg_type//achar(27)//"[0m]: "//msg
                    else if (msg_type == "Warning") then
                        print *, "["//achar(27)//"[1;95m"//msg_type//achar(27)//"[0m]: "//msg
                    else
                        print *, "["//achar(27)//"[1;32m"//msg_type//achar(27)//"[0m]: "//msg
                    end if

                    ! print variable if var is given
                    if (present(var)) then
                        ! switch over several cases for the type of the variable
                        select type(var)
                            type is (logical)
                                print *, var, "logical variable"
                            type is (integer(2))
                                print *, var, "integer(2) variable"
                            type is (integer(4))
                                print *, var, "integer(4) variable"
                            type is (real(4))
                                print *, var, "real(4) variable"
                            type is (real(8))
                                print *, var, "real(8) variable"
                            type is (complex(4))
                                print *, var, "complex(8) variable"
                            type is (complex(8))
                                print *, var, "complex(16) variable"
                            type is (character(*))
                                print *, var, "character(*) variable"
                        end select
                    end if

                ! end if block for debug
                end if

                ! stop execution if trg_stop_ flag is enabled
                if (trg_stop_) stop

            ! end if block for condition
            end if
        end subroutine check


        ! check if dimensions of input matrices are the same
        function dim_check(mat1, mat2) result(AreDimEq)
            real(8), dimension(:,:) :: mat1, mat2
            integer(4) :: N1, M1, N2, M2
            logical :: AreDimEq


            ! mat1 and mat2 dimensions
            N1 = size(mat1, 1)
            M1 = size(mat1, 2)
            N2 = size(mat2, 1)
            M2 = size(mat2, 2)

            if ((N1.EQ.N2).AND.(M1.EQ.M2)) then
                AreDimEq = .TRUE.
            else
                AreDimEq = .FALSE.
            end if
        end function dim_check


        ! check if input matrices entries are the same at a precision level thld
        function eq_check(mat1, mat2, thld) result(AreMatEq)
            logical :: AreMatEq
            integer(4) :: N1, M1
            integer(4) :: ii, jj
            ! mandatory arguments
            real(8), dimension(:,:) :: mat1, mat2
            ! optional arguments
            real(8), optional :: thld
            ! dummy arguments
            real(8) :: thld_


            ! mat1 dimensions
            N1 = size(mat1, 1)
            M1 = size(mat1, 2)

            ! 1d-6 default threshold for equality check
            if (.NOT.present(thld)) then
                thld_ = 1d-06
            else
                thld_ = thld
            end if

            AreMatEq = .TRUE.
            ! check if dimensions of mat1 and mat2 are the same
            ! otherwise return false
            if (dim_check(mat1,mat2)) then
                ! loop over all matrix elements
                do jj=1,M1
                    do ii=1,N1
                        ! check if equality holds at a certain precision depending on thld
                        ! otherwise exit the loop and return false
                        if (ABS(mat1(ii,jj) - mat2(ii,jj)) < 0.5*thld_*ABS(mat1(ii,jj) + mat2(ii,jj))) then
                            continue
                        else
                            AreMatEq = .FALSE.
                            exit
                        end if
                    end do
                end do
            else
                AreMatEq = .FALSE.
            end if
        end function eq_check
end module err_handling





! module ising_utils for Ising model utilities
module ising_utils
    implicit none


    contains
        ! function to execute tensor product between matrices
        function ising_tensor_prod(mat1, mat2) result(res)
            ! input arguments
            complex(8), dimension(:,:) :: mat1
            complex(8), dimension(:,:) :: mat2

            ! output
            complex(8), dimension(:,:), allocatable :: res

            integer(4) :: N1, N2, M1, M2
            integer(4) :: ii, jj, lli, llj, uui, uuj

            N1 = size(mat1,1)
            N2 = size(mat2,1)
            M1 = size(mat1,2)
            M2 = size(mat2,2)

            allocate(res(N1*N2,M1*M2))

            ! loop to execute tensor product between mat1 and mat2
            do ii=1,N1
                do jj=1,M1
                    lli = (ii-1)*N2 + 1
                    llj = (jj-1)*M2 + 1
                    uui =  ii   *N2
                    uuj =  jj   *M2
                    res(lli:uui,llj:uuj) = mat1(ii,jj)*mat2
                end do
            end do
        end function ising_tensor_prod


        ! function to create an identity matrix for N particles
        function ising_identity(N) result(id)
            ! input arguments
            integer(4) :: N

            ! output
            complex(8), dimension(2**N,2**N) :: id

            integer(4) :: ii

            id = COMPLEX(0.0d0, 0.0d0)
            do ii=1,size(id,1)
                id(ii,ii) = COMPLEX(1.0d0,0.0d0)
            end do
        end function


        ! function for initialization of Ising Hamiltonian
        function ising_hmat_init(N, L) result(hmat)
            ! input arguments
            integer(4) :: N
            real(8)    :: L

            ! output
            complex(8), dimension(2**N,2**N) :: hmat

            complex(8), dimension(2,2) :: s_x, s_z
            integer(4) :: ii

            s_x      = COMPLEX( 0.0d0, 0.0d0)
            s_x(1,2) = COMPLEX( 1.0d0, 0.0d0)
            s_x(2,1) = COMPLEX( 1.0d0, 0.0d0)
            s_z      = COMPLEX( 0.0d0, 0.0d0)
            s_z(1,1) = COMPLEX( 1.0d0, 0.0d0)
            s_z(2,2) = COMPLEX(-1.0d0, 0.0d0)

            ! initialize Hamiltonian to zero
            hmat = COMPLEX(0.0d0, 0.0d0)
            ! fill with the first piece of Ising Hamiltonian
            do ii=1,N
                hmat = hmat + ising_tensor_prod(                  &
                    ising_tensor_prod(ising_identity(ii-1), s_z), &
                    ising_identity(N-ii)                          &
                )
            end do
            ! multiply by lambda the first piece
            hmat = COMPLEX(L,0.0d0) * hmat
            ! fill with the second piece of Ising Hamiltonian
            do ii=1,N-1
                hmat = hmat + ising_tensor_prod(   &
                    ising_tensor_prod(                                &
                        ising_tensor_prod(ising_identity(ii-1), s_x), &
                        s_x                                           &
                    ),                                                &
                    ising_identity(N-ii-1)                            &
                )
            end do
        end function


        ! diagonalize Ising Hamiltonian
        subroutine ising_hmat_diag(hmat, eigs)
            ! input arguments
            complex(8), dimension(:,:) :: hmat
            real(8),    dimension(:)   :: eigs

            complex(8), dimension(:), allocatable :: work(:)
            complex(8), dimension(:), allocatable :: rwork(:)
            integer(4) :: lwork, info
            integer(4) :: N

            N     = size(hmat, 1)
            lwork = max(1,2*N-1)

            allocate(work(lwork))
            allocate(rwork(max(1, 2*N)))
            call zheev('V', 'U', N, hmat, N, eigs, work, lwork, rwork, info)
            deallocate(work)
            deallocate(rwork)
        end subroutine ising_hmat_diag


        ! print Ising Hamiltonian on standard output
        subroutine ising_print_hmat_std(hmat, formatted)
            ! input arguments
            complex(8), dimension(:,:) :: hmat
            logical                    :: formatted

            integer(4) :: ii, jj

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(hmat,1)
                if (formatted) then
                    write(*, "(*('('sf7.4xspf7.4'i)':x))") (hmat(ii,jj), jj=1,size(hmat,2))
                else
                    write(*, *) (hmat(ii,jj), jj=1,size(hmat,2))
                end if
            end do
        end subroutine ising_print_hmat_std


        ! print Ising Hamiltonian eigenvalues on standard output
        subroutine ising_print_hmat_eigs_std(eigs, formatted)
            ! input arguments
            real(8), dimension(:) :: eigs
            logical               :: formatted

            integer(4) :: ii

            ! write vector elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(eigs,1)
                if (formatted) then
                    write(*, "(*(sf7.4:x))") eigs(ii)
                else
                    write(*, *) eigs(ii)
                end if
            end do
        end subroutine ising_print_hmat_eigs_std
end module ising_utils






program cmdline
    use err_handling
    use ising_utils
    implicit none

    character(20) :: arg

    integer(4)    :: N, K                   ! meaning explained later
    real(8)       :: L
    integer(4)    :: io_stat                ! flag signaling correct reading
    logical       :: def_N, def_L, def_K    ! flags for default args
    character(20) :: str_N, str_L, str_K
    integer(4)    :: ii

    ! set default argument for all arguments
    ! then parse on the arguments given
    ! and unset default argument for arguments
    ! whose value is given
    def_N = .TRUE.
    def_L = .TRUE.
    def_K = .TRUE.
    ! just a machinery to prettify output
    str_N = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_L = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_K = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"

    do ii=1,command_argument_count()
        call get_command_argument(ii, arg)

        select case (arg)
        case ('-N', '--nspins')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) N
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for nspins", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_N = .FALSE.

            case ('-L',  '--lambda')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) L
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for lambda", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_L = .FALSE.

            case ('-K',  '--levels')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) K
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for levels", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_K = .FALSE.

            case ('-h',  '--help')
                call print_help()
                stop
        end select
    end do

    ! set default argument if no option is given for an argument
    if (def_N .EQV. .TRUE.) N = 2
    if (def_L .EQV. .TRUE.) L = 1.0d0
    if (def_K .EQV. .TRUE.) K = 1
    ! set red color for output if default is not given
    ! set green color for output if default is given
    if (def_N .EQV. .TRUE.) str_N = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_L .EQV. .TRUE.) str_L = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_K .EQV. .TRUE.) str_K = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"


    call main(N, L, K)


    contains
        ! print help on stdout (explanation of the possible arguments)
        subroutine print_help()
            print *, 'usage: cmdline [OPTIONS]'
            print *, ''
            print *, 'cmdline options:'
            print *, ''
            print *, char(9), '-N, --nspins    [default=2]    number N of spins'
            print *, char(9), '-L, --lambda    [default=1]    value of interaction strength lambda'
            print *, char(9), '-K, --levels    [default=1]    first K levels (eigenvalues) to print'
        end subroutine print_help


        ! print arguments on stdout
        subroutine print_args()!(NX, L, NT, T, TF, M, O)
            print *, 'Running with the following cmdline options:'
            print *, ''
            write(*, "(a,a,g8.2,a,a,g0)") 'N =', char(9), N, char(9), 'default [T/F]: ', str_N
            write(*, "(a,a,g8.2,a,a,g0)") 'L =', char(9), L, char(9), 'default [T/F]: ', str_L
            write(*, "(a,a,g8.2,a,a,g0)") 'K =', char(9), K, char(9), 'default [T/F]: ', str_K
            write(*, "(a)") ''
        end subroutine print_args


        ! mode 1 program
        subroutine main(N, L, K)
            ! input arguments
            integer(4) :: N
            real(8)    :: L
            integer(4) :: K

            complex(8), dimension(2**N,2**N) :: hmat
            real(8),    dimension(2**N)      :: eigs
            integer(4) :: ii

            call print_args()

            eigs = 0.0d0
            hmat = ising_hmat_init(N,L)
            ! call ising_print_hmat_std(hmat, .TRUE.)

            call ising_hmat_diag(hmat, eigs)
            ! write(*, "(a)") ''
            ! call ising_print_hmat_std(hmat, .TRUE.)
            ! write(*, "(a)") ''
            ! call ising_print_hmat_eigs_std(eigs, .TRUE.)

            do ii=1,K
                if (ii<=2**N) then
                    write(*, "(g0)") eigs(ii)
                else
                    write(*, "(a)")  "****************"
                end if
            end do
        end subroutine main

end program cmdline

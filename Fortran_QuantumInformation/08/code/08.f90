! 08.f90
! gfortran -Wall -ffree-line-length-0 -o 08.o 08.f90





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


            ! assign defaul value to trg_stop argument through dummy variable
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





module dmat
    implicit none

    contains
        ! initialize random state
        function dmat_init_rand_state(N, D, isSep) result(state)
            ! input arguments
            integer(4) :: N
            integer(4) :: D
            logical    :: isSep

            integer(4) :: dim
            complex(8), dimension(:), allocatable :: state
            real(8),    dimension(:), allocatable :: re, im
            real(8)    :: norm

            if (isSep) then
                dim = D*N
            else
                dim = D**N
            end if

            allocate(state(dim))
            allocate(re(dim))
            allocate(im(dim))

            ! fill real and imaginary part vectors
            call RANDOM_NUMBER(re)
            call RANDOM_NUMBER(im)
            re = re*2.0d0 - 1.0d0
            im = im*2.0d0 - 1.0d0

            ! fill state vector with coefficients and normalize
            norm  = SUM(re**2 + im**2)
            state = (re*COMPLEX(1.0d0,0.0d0) + im*COMPLEX(0.0d0,1.0d0)) / COMPLEX(SQRT(norm),0.0d0)

            ! do ii=1,dim
            !     call RANDOM_NUMBER(re)
            !     call RANDOM_NUMBER(im)
            !     state(ii) = COMPLEX(re*2.0d0-1.0d0, im*2.0d0-1.0d0)
            ! end do
        end function dmat_init_rand_state


        ! initialize density matrix for a random state
        function dmat_init_rho(state) result(rho)
            ! input arguments
            complex(8), dimension(:) :: state

            complex(8), dimension(:,:), allocatable :: rho
            complex(8), dimension(:,:), allocatable :: bra, ket
            integer(4) :: N


            N = size(state,1)

            allocate(rho(N,N))
            allocate(bra(1,N))
            allocate(ket(N,1))

            bra(1,:) = CONJG(state)
            ket(:,1) =       state

            rho = MATMUL(ket,bra)
        end function dmat_init_rho


        ! initialize density matrix for a random state
        function dmat_rho_reduce(rho, nsub, N, D) result(rho_red)
            ! input arguments
            complex(8), dimension(:,:) :: rho
            integer(4) :: nsub
            integer(4) :: N
            integer(4) :: D

            complex(8), dimension(D**(N-1),D**(N-1)) :: rho_red
            complex(8) :: sum
            integer(4) :: lli, llj, mmi, mmj, kk, ii, jj, iir, jjr

            ! loop over rho elements of all subsystems except nsub
            do lli=0,D**(nsub-1)-1
                do mmi=0,(D**(N-nsub))-1
                    do llj=0,D**(nsub-1)-1
                        do mmj=0,(d**(N-nsub))-1
                            ! trace the nsub subsystem
                            sum = COMPLEX(0.0d0,0.0d0)
                            do kk=0,D-1
                                ii = lli + 1 + (kk + mmi*D)*D**(nsub-1)
                                jj = llj + 1 + (kk + mmj*D)*D**(nsub-1)
                                sum = sum + rho(ii,jj)
                            end do
                            ! fill reduced density matrix
                            iir = lli + 1 + mmi*D**(nsub-1)
                            jjr = llj + 1 + mmj*D**(nsub-1)
                            rho_red(iir, jjr) = sum
                        end do
                    end do
                end do
            end do
        end function dmat_rho_reduce


        ! print random state on standard output
        subroutine dmat_print_state_std(state, formatted)
            ! input arguments
            complex(8), dimension(:) :: state
            logical                  :: formatted

            integer(4) :: ii

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(state,1)
                if (formatted) then
                    write(*, "('('sf7.4xspf7.4'i)')") state(ii)
                else
                    write(*, *) state(ii)
                end if
            end do
        end subroutine dmat_print_state_std


        ! print random state density matrix on standard output
        subroutine dmat_print_rho_std(rho, formatted)
            ! input arguments
            complex(8), dimension(:,:) :: rho
            logical                    :: formatted

            integer(4) :: ii, jj

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(rho,1)
                if (formatted) then
                    write(*, "(*('('sf7.4xspf7.4'i)':x))") (rho(ii,jj), jj=1,size(rho,2))
                else
                    write(*, *) (rho(ii,jj), jj=1,size(rho,2))
                end if
            end do
        end subroutine dmat_print_rho_std
end module dmat




program cmdline
    use err_handling
    use dmat
    implicit none

    character(20) :: arg

    integer(4)    :: N, D                                 ! meaning explained later
    integer(4)    :: S, M
    integer(4)    :: K
    integer(4)    :: io_stat                              ! flag signaling correct reading
    logical       :: def_N, def_D, def_S, def_M, def_K    ! flags for default args
    character(20) :: str_N, str_D, str_S, str_M, str_K
    integer(4) :: ii

    ! set default argument for all arguments
    ! then parse on the arguments given
    ! and unset default argument for arguments
    ! whose value is given
    def_N = .TRUE.
    def_D = .TRUE.
    def_S = .TRUE.
    def_M = .TRUE.
    def_K = .TRUE.
    ! just a machinery to prettify output
    str_N = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_D = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_S = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_M = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_K = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"

    do ii=1,command_argument_count()
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-N', '--subsysN')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) N
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for subsysN", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_N = .FALSE.

            case ('-D',  '--dimD')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) D
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for dimD", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_D = .FALSE.

            case ('-S',  '--separable')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) S
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for separable", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_S = .FALSE.

            case ('-M',  '--mode')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) M
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for mode", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_M = .FALSE.

            case ('-K',  '--reduceK')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) K
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for reduceK", &
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
    if (def_D .EQV. .TRUE.) D = 2
    if (def_S .EQV. .TRUE.) S = 0
    if (def_M .EQV. .TRUE.) M = 1
    if (def_K .EQV. .TRUE.) K = 1
    ! set red color for output if default is not given
    ! set green color for output if default is given
    if (def_N .EQV. .TRUE.) str_N = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_D .EQV. .TRUE.) str_D = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_S .EQV. .TRUE.) str_S = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_M .EQV. .TRUE.) str_M = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_K .EQV. .TRUE.) str_K = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"


    select case (M)
        case (1)
            call mode_1(N, D, S)
        case (2)
            call mode_2(N, D, S, K)
        case (3)
            call mode_3()
        case default
            print *, "Wrong execution mode given"
            call print_help()
            stop
    end select


    contains
        ! print help on stdout (explanation of the possible arguments)
        subroutine print_help()
            print *, 'usage: cmdline [OPTIONS]'
            print *, ''
            print *, 'cmdline options:'
            print *, ''
            print *, char(9), '-N, --subsysN      [default=10]    number N of subsystems'
            print *, char(9), '-D, --dimD         [default=10]    dimension D of every subsystem wavefunction'
            print *, char(9), '-S, --separable    [default= 0]    istantiate non-separable (0) or separable (else) state'
            print *, char(9), '-M, --mode         [default= 1]    select execution mode (1/2/3)'
            print *, char(9), '-K, --reduceK      [default= 1]    index of subsystem to trace (from 1 to N)'
        end subroutine print_help


        ! print arguments on stdout
        subroutine print_args()!(NX, L, NT, T, TF, M, O)
            print *, 'Running with the following cmdline options:'
            print *, ''
            write(*, "(a,a,g8.2,a,a,g0)") 'N =', char(9), N, char(9), 'default [T/F]: ', str_N
            write(*, "(a,a,g8.2,a,a,g0)") 'D =', char(9), D, char(9), 'default [T/F]: ', str_D
            write(*, "(a,a,g8.2,a,a,g0)") 'S =', char(9), S, char(9), 'default [T/F]: ', str_S
            write(*, "(a,a,g8.2,a,a,g0)") 'M =', char(9), M, char(9), 'default [T/F]: ', str_M
            write(*, "(a,a,g8.2,a,a,g0)") 'K =', char(9), K, char(9), 'default [T/F]: ', str_K
            write(*, "(a)") ''
        end subroutine print_args


        ! mode 1 program
        subroutine mode_1(N, D, S)
            ! input arguments
            integer(4) :: N
            integer(4) :: D
            integer(4) :: S

            complex(8), dimension(:),   allocatable :: state
            ! complex(8), dimension(:,:), allocatable :: rho
            ! complex(8), dimension(:,:), allocatable :: rho_red
            real(8)    :: t_i, t_f
            integer(4) :: dim, dim_red


            call print_args()

            if (S==1) then
                dim     = D*N
                dim_red = D*(N-1)
            else
                dim     = D**N
                dim_red = D**(N-1)
            end if

            allocate(state(dim))

            call CPU_TIME(t_i)
            state   = dmat_init_rand_state(N, D, S.NE.0)
            call CPU_TIME(t_f)

            print *, "CPU time:"
            print *, t_f-t_i
        end subroutine mode_1


        ! mode 2 program
        subroutine mode_2(N, D, S, K)
            ! input arguments
            integer(4) :: N
            integer(4) :: D
            integer(4) :: S
            integer(4) :: K

            complex(8), dimension(:),   allocatable :: state
            complex(8), dimension(:,:), allocatable :: rho
            complex(8), dimension(:,:), allocatable :: rho_red
            integer(4) :: dim, dim_red


            call print_args()

            if (S==1) then
                dim     = D*N
                dim_red = D*(N-1)
            else
                dim     = D**N
                dim_red = D**(N-1)
            end if

            allocate(state(dim))

            state = dmat_init_rand_state(N, D, S.NE.0)

            write(*, "(a)") "State coefficients:"
            call dmat_print_state_std(state, .TRUE.)
            write(*, "(a)") ''
            if (S.EQ.0) then
                rho     = dmat_init_rho(state)
                rho_red = dmat_rho_reduce(rho, K, N, D)

                write(*, "(a)") "Density matrix:"
                call dmat_print_rho_std(rho, .TRUE.)
                write(*, "(a)") ''
                write(*, "(ai2a)") "Reduced density matrix (k=", K, "):"
                call dmat_print_rho_std(rho_red, .TRUE.)
            end if
        end subroutine mode_2


        ! mode 2 program
        subroutine mode_3()
            integer(4) :: N
            integer(4) :: D
            integer(4) :: S

            complex(8), dimension(:),   allocatable :: state
            complex(8), dimension(:,:), allocatable :: rho
            complex(8), dimension(:,:), allocatable :: rho_red_l
            complex(8), dimension(:,:), allocatable :: rho_red_r

            N = 2
            D = 2
            S = 0

            allocate(state(N**D))

            state = dmat_init_rand_state(N, D, S.NE.0)

            write(*, "(a)") ''
            write(*, "(a)") "Running for two spin one-half qubits:"
            write(*, "(a)") ''
            write(*, "(a)") "State coefficients:"
            call dmat_print_state_std(state, .TRUE.)
            write(*, "(a)") ''
            if (S.EQ.0) then
                rho       = dmat_init_rho(state)
                rho_red_l = dmat_rho_reduce(rho, 1, N, D)
                rho_red_r = dmat_rho_reduce(rho, 2, N, D)

                write(*, "(a)") "Density matrix:"
                call dmat_print_rho_std(rho, .TRUE.)
                write(*, "(a)") ''
                write(*, "(ai2a)") "Reduced density matrix (on left system):"
                call dmat_print_rho_std(rho_red_l, .TRUE.)
                write(*, "(a)") ''
                write(*, "(ai2a)") "Reduced density matrix (on right system):"
                call dmat_print_rho_std(rho_red_r, .TRUE.)
            end if
        end subroutine mode_3

end program cmdline

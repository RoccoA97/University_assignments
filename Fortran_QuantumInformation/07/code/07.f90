! 07.f90
! gfortran -Wall -ffree-line-length-0 -o 07.o 07.f90 -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack -L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3f
! ./07.0 [options] (tipe -h for help)





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





! module for time-dependent Schrodinger equation utilities
module tdep_se_utils
    use, intrinsic :: iso_c_binding
    implicit none

    contains
        ! initialize hamiltonian operator
        function ham_init(N, L, M, O, t0, T) result(ham)
            ! input arguments
            integer(4) :: N     ! hamiltonian dimesnions
            real(8)    :: L     ! width of the box
            real(8)    :: M     ! mass
            real(8)    :: O     ! omega
            real(8)    :: t0    ! initial istant of time
            real(8)    :: T     ! characteristic time

            real(8), dimension(N+1,N+1) :: ham
            real(8) :: DX
            integer(4) :: ii, jj

            DX = L / N
            do ii=1,N+1
                do jj=1,N+1
                    if (ABS(ii-jj)==1) then
                        ham(ii,jj) = - (1.0d0/(2.0d0*M)) * (1.0d0/(DX**2))
                    else if (ii==jj) then
                        ham(ii,jj) =   (1.0d0/(2.0d0*M)) * (2.0d0/(DX**2)) + (0.5d0*M*O**2) * (-L/2 + DX*(ii-1) - t0/T)**2
                    else
                        ham(ii,jj) =    0.0d0
                    end if
                end do
            end do
        end function ham_init


        ! diagonalize hamiltonian operator
        subroutine ham_diag(ham, eigs)
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
        end subroutine ham_diag


        ! kinetic term of Hamiltonian
        function ham_T(N, L, M) result(T)
            ! input arguments
            integer(4) :: N
            real(8)    :: L
            real(8)    :: M

            real(8), dimension(N+1) :: T
            real(8)    :: PI
            real(8)    :: K, DP, DX
            integer(4) :: ii

            PI = ACOS(-1.0d0)
            K  = 2*PI / L
            DP = K / N
            DX = L / N

            do ii=1,INT(N/2)
                T(ii) = (1.0d0/(2.0d0*M)) * (2.0d0*PI*(ii)/L)**2
            end do

            do ii=INT(N/2)+1,N+1
                T(ii) = (1.0d0/(2.0d0*M)) * (2.0d0*PI*(ii-N-1)/L)**2
            end do
        end function ham_T


        ! potential term of Hamiltonian
        function ham_V(N, L, M, O, time, TT) result(V)
            ! input arguments
            integer(4) :: N
            real(8)    :: L, M, O, time, TT


            real(8), dimension(N+1) :: V
            real(8)    :: DX
            integer(4) :: ii

            DX = L / N

            do ii=1,size(V,1)
                V(ii) = (1.0d0/2.0d0) * M * O**2 * (-L/2.0 + (ii-1)*DX - time/TT)**2
            end do

            ! well potential
            ! do ii=1,size(V,1)
            !     if ((ii==1).OR.(ii==size(V,1))) then
            !         V(ii) = 0.0d0
            !     else
            !         V(ii) = 1.0d30
            !     end if
            ! end do
        end function ham_V


        ! fast fourier transform
        function fft(psi) result(FFT_psi)
            ! input arguments
            complex(8), dimension(:) :: psi

            complex(8), dimension(size(psi,1)) :: FFT_psi
            integer(4) :: N
            integer(8) :: plan

            N = size(psi,1)

            call dfftw_plan_dft_1d(plan, N, psi, FFT_psi, -1, 64)
            call dfftw_execute_dft(plan, psi, FFT_psi)
            call dfftw_destroy_plan(plan)
        end function fft


        ! anti fast fourier transform
        function afft(psi) result(AFFT_psi)
            ! input arguments
            complex(8), dimension(:) :: psi

            complex(8), dimension(size(psi,1)) :: AFFT_psi
            integer(4) :: N
            integer(8) :: plan

            N = size(psi,1)

            call dfftw_plan_dft_1d(plan, N, psi, AFFT_psi, +1, 64)
            call dfftw_execute_dft(plan, psi, AFFT_psi)
            call dfftw_destroy_plan(plan)
        end function afft


        ! time evolution of a step DT
        subroutine time_ev(psi, DT, L, M, O, time, TT)
            ! input arguments
            complex(8), dimension(:) :: psi
            real(8) :: DT, time, TT
            real(8) :: L, M, O

            real(8), dimension(size(psi,1)) :: V
            real(8), dimension(size(psi,1)) :: T
            integer(4) :: N
            integer(4) :: ii

            N = size(psi,1) - 1

            ! compute potential and kinetic terms
            V = ham_V(N, L, M, O, time, TT)     ! potential term
            T = ham_T(N, L, M)                  ! kinetic   term

            ! apply exponential containing the potential term
            do ii=1,N+1
                psi(ii) = ZEXP(COMPLEX(0.0d0,-0.5d0*DT*V(ii))) * psi(ii)
            end do
            ! fft of wave function and normalization
            psi = fft(psi) / SQRT(1.0d0*N+1)
            ! apply exponential containing the kinetic term
            do ii=1,N+1
                psi(ii) = ZEXP(COMPLEX(0.0d0,-1.0d0*DT*T(ii))) * psi(ii)
            end do
            ! afft of wave function and normalization
            psi = afft(psi) / SQRT(1.0d0*N+1)
            ! apply exponential containing the potential term
            do ii=1,N+1
                psi(ii) = ZEXP(COMPLEX(0.0d0,-0.5d0*DT*V(ii))) * psi(ii)
            end do
        end subroutine time_ev


        ! print wavefunction on file
        subroutine print_wvfc(xs, ys, filename, unit)
            ! input arguments
            real(8),    dimension(:) :: xs
            complex(8), dimension(:) :: ys
            character(*) :: filename
            integer(4) :: unit

            integer(4) :: ii

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(xs,1)
                write(unit,'(g0ag0ag0)') xs(ii), char(9), REAL(ys(ii)), char(9), AIMAG(ys(ii))
            end do
            close(unit)
        end subroutine print_wvfc


        ! print wavefunction evolution on file in cmap format
        subroutine print_wvfc_cmap(xs, ys, zs, filename, unit)
            ! input arguments
            real(8),    dimension(:)   :: xs
            real(8),    dimension(:)   :: ys
            complex(8), dimension(:,:) :: zs
            character(*) :: filename
            integer(4) :: unit

            integer(4) :: ii, jj

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(ys,1)
                do jj=1, size(xs,1)
                    write(unit,'(g0ag0ag0)') xs(jj), char(9), ys(ii), char(9), REAL(zs(jj,ii)*CONJG(zs(jj,ii)))
                end do
            end do
            close(unit)
        end subroutine print_wvfc_cmap


        ! print wavefunction exp. value evolution on file
        subroutine print_wvfc_expv(xs, ys, zs, filename, unit)
            ! input arguments
            real(8),    dimension(:)   :: xs
            real(8),    dimension(:)   :: ys
            complex(8), dimension(:,:) :: zs
            character(*) :: filename
            integer(4) :: unit

            real(8)    :: expv, flct
            integer(4) :: ii, jj

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(ys,1)
                ! compute expectation value
                expv = 0.0d0
                do jj=1, size(xs,1)
                    expv = expv + xs(jj) * REAL(zs(jj,ii)*CONJG(zs(jj,ii))) * (xs(2)-xs(1))
                end do
                ! compute fluctuation
                flct = 0.0d0
                do jj=1,size(xs,1)
                    flct = flct + (REAL(zs(jj,ii)*CONJG(zs(jj,ii))) * (xs(2)-xs(1))) * (xs(jj) - expv)**2
                end do
                flct = SQRT(flct)
                ! write expectation value and expectation value +/- sigma on file
                write(unit,'(g0ag0ag0ag0ag0ag0ag0ag0)') ys(ii), char(9), expv, char(9), expv-      flct, char(9), expv+      flct, char(9), &
                                                                                        expv-2.0d0*flct, char(9), expv+2.0d0*flct, char(9), &
                                                                                        expv-3.0d0*flct, char(9), expv+3.0d0*flct
            end do
            close(unit)
        end subroutine print_wvfc_expv
end module tdep_se_utils





program cmdline
    use tdep_se_utils
    use err_handling
    implicit none

    character(20) :: arg

    integer(4)    :: NX, NT                                                ! meaning explained later
    real(8)       :: L, T, TF                                              ! meaning explained later
    real(8)       :: M, O                                                  ! meaning explained later
    integer(4)    :: io_stat                                               ! flag signaling correct reading
    logical       :: def_NX, def_L, def_NT, def_T, def_TF, def_M, def_O    ! flags for default args
    character(20) :: str_NX, str_L, str_NT, str_T, str_TF, str_M, str_O
    integer(4) :: ii

    ! set default argument for all arguments
    ! then parse on the arguments given
    ! and unset default argument for arguments
    ! whose value is given
    def_NX = .TRUE.
    def_L  = .TRUE.
    def_NT = .TRUE.
    def_T  = .TRUE.
    def_TF = .TRUE.
    def_M  = .TRUE.
    def_O  = .TRUE.
    ! just a machinery to prettify output
    str_NX = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_L  = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_NT = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_T  = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_TF = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_M  = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"
    str_O  = ""//achar(27)//"[1;31m"//"F"//achar(27)//"[0m"

    do ii=1,command_argument_count()
        call get_command_argument(ii, arg)

        select case (arg)
            case ('-NX', '--stepsX')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) NX
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for stepsX", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_NX = .FALSE.

            case ('-L',  '--lengthX')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) L
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for lengthX", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_L  = .FALSE.

            case ('-NT', '--stepsT')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) NT
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for stepsT", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_NT = .FALSE.

            case ('-T',  '--lengthT')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) T
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for lengthT", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_T  = .FALSE.

            case ('-TF', '--finalT')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) TF
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for finalT", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_TF = .FALSE.

            case ('-M',  '--mass')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) M
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for mass", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_M  = .FALSE.

            case ('-O',  '--omega')
                call get_command_argument(ii+1, arg)
                READ(arg, *, iostat=io_stat) O
                ! check validity of argument
                call check(debug     = .TRUE., &
                           condition = (io_stat.NE.0), &
                           msg_type  = "Error", &
                           msg       = "Invalid parameter value for omega", &
                           trg_stop  = .TRUE.)
                ! disable default value for argument
                def_O  = .FALSE.

            case ('-h',  '--help')
                call print_help()
                stop
        end select
    end do

    ! set default argument if no option is given for an argument
    if (def_NX .EQV. .TRUE.) NX = 1000
    if (def_L  .EQV. .TRUE.) L  = 10.0d0
    if (def_NT .EQV. .TRUE.) NT = 1000
    if (def_T  .EQV. .TRUE.) T  = 5.0d0
    if (def_TF .EQV. .TRUE.) TF = 5.0d0
    if (def_M  .EQV. .TRUE.) M  = 1.0d0
    if (def_O  .EQV. .TRUE.) O  = 1.0d0
    ! set red color for output if default is not given
    ! set green color for output if default is given
    if (def_NX .EQV. .TRUE.) str_NX = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_L  .EQV. .TRUE.) str_L  = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_NT .EQV. .TRUE.) str_NT = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_T  .EQV. .TRUE.) str_T  = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_TF .EQV. .TRUE.) str_TF = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_M  .EQV. .TRUE.) str_M  = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"
    if (def_O  .EQV. .TRUE.) str_O  = ""//achar(27)//"[1;32m"//"T"//achar(27)//"[0m"

    call main(NX, L, NT, T, TF, M, O)


    contains
        ! print help on stdout (explanation of the possible arguments)
        subroutine print_help()
            print *, 'usage: cmdline [OPTIONS]'
            print *, ''
            print *, 'cmdline options:'
            print *, ''
            print *, char(9), '-NX, --stepsX     [default=1000]    number of spatial discretization intervals'
            print *, char(9), '-L,  --lengthX    [default=  10]    length of symmetric interval of simulation'
            print *, char(9), '-NT, --stepsT     [default=1000]    number of time discretization intervals'
            print *, char(9), '-T,  --lengthT    [default=   5]    characteristic time of potential (denominator)'
            print *, char(9), '-TF, --finalT     [default=   5]    final time of simulation'
            print *, char(9), '-M,  --mass       [default=   1]    mass parameter of hamiltonian'
            print *, char(9), '-O,  --omega      [default=   1]    omega parameter of hamiltonian'
        end subroutine print_help


        ! print arguments on stdout
        subroutine print_args()!(NX, L, NT, T, TF, M, O)
            print *, 'Running with the following cmdline options:'
            print *, ''
            write(*, "(a,a,g8.2,a,a,g0)") 'NX =', char(9), NX, char(9), 'default [T/F]: ', str_NX
            write(*, "(a,a,g8.2,a,a,g0)") 'L  =', char(9), L,  char(9), 'default [T/F]: ', str_L
            write(*, "(a,a,g8.2,a,a,g0)") 'NT =', char(9), NT, char(9), 'default [T/F]: ', str_NT
            write(*, "(a,a,g8.2,a,a,g0)") 'T  =', char(9), T,  char(9), 'default [T/F]: ', str_T
            write(*, "(a,a,g8.2,a,a,g0)") 'TF =', char(9), TF, char(9), 'default [T/F]: ', str_TF
            write(*, "(a,a,g8.2,a,a,g0)") 'M  =', char(9), M,  char(9), 'default [T/F]: ', str_M
            write(*, "(a,a,g8.2,a,a,g0)") 'O  =', char(9), O,  char(9), 'default [T/F]: ', str_O
        end subroutine print_args


        ! main program
        subroutine main(NX, L, NT, T, TF, M, O)
            ! use tdep_se_utils
            ! implicit none

            real(8),    dimension(:,:), allocatable :: ham
            real(8),    dimension(:),   allocatable :: eigs
            real(8),    dimension(:),   allocatable :: xs
            real(8),    dimension(:),   allocatable :: ts
            complex(8), dimension(:),   allocatable :: psi
            complex(8), dimension(:,:), allocatable :: psi_all

            ! integer(4) :: arg1, arg3
            ! real(8)    :: arg2, arg4, arg5, arg6
            real(8)    :: L, DX
            real(8)    :: M, O
            real(8)    :: T0, tt, T, TF, DT
            integer(4) :: NX, NT
            integer(4) :: ii, jj
            character(4) :: jj_str

            call print_args()

            T0 = 0.0d0

            DX = L  / NX
            DT = (TF-T0) / NT


            ! allocate memory to store vectors and matrices
            allocate(ham(NX+1,NX+1))
            allocate(eigs(NX+1))
            allocate(psi(NX+1))
            allocate(xs(NX+1))
            allocate(ts(NT+1))
            allocate(psi_all(NX+1,NT+1))


            ! initialize discretized Hamiltonian with
            ! * N discretization intervals
            ! * L length of space interval [-a,a]
            ! * M mass and O omega parameters
            ham = ham_init(NX, L, M, O, t0, T)

            ! diagonalize hamiltonian
            call ham_diag(ham, eigs)

            ! define space and time grids
            do ii=1,size(xs,1)
                xs(ii) = (-L/2) + (ii-1)*DX
            end do
            do ii=1,size(ts,1)
                ts(ii) = (ii-1)*DT
            end do

            ! normalize eigenfunctions
            do ii=1,size(ham,1)
                do jj=1,size(ham,2)
                    ham(ii,jj) = ham(ii,jj) * SQRT(NX/L)
                end do
            end do

            ! extract ground state at t0 and save points
            tt = 0.0d0
            jj = 0
            do ii=1,size(ham,1)
                psi(ii) = COMPLEX(ham(ii,1), 0.0d0)
            end do
            psi_all(:,1) = psi
            write(jj_str, "(i1)") jj
            call print_wvfc(xs, psi, "results/"//TRIM(jj_str)//"_step_wvfc.dat", 1)

            ! iterate and save points at time step ii
            do ii=1,NT
                tt = tt  + DT
                jj = jj + 1
                if (ii < 10) then
                    write(jj_str, "(i1)") jj
                else if (ii < 100) then
                    write(jj_str, "(i2)") jj
                else if (ii < 1000) then
                    write(jj_str, "(i3)") jj
                else
                    write(jj_str, "(i4)") jj
                end if
                call time_ev(psi, DT, L, M, O, tt, T)
                psi_all(:,ii+1) = psi
                call print_wvfc(xs, psi, "results/"//TRIM(jj_str)//"_step_wvfc.dat", 1)
            end do

            ! save data for colormap plot
            ! through which time evolution of wavefunction
            ! can be visualized in the entire space grid
            call print_wvfc_cmap(xs, ts, psi_all, "results/wvfc_cmap.dat", 1)
            call print_wvfc_expv(xs, ts, psi_all, "results/wvfc_expv.dat", 1)
        end subroutine main

end program cmdline

! 05.f90
! gfortran -Wall -ffree-line-length-0 -o 05.o 05.f90 -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack





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





! module for hermiatian matrix handling and eigenvalues analysis
module hmat
    use err_handling
    implicit none

    contains
        ! generate random normal number
        function random_normal() result(e)

            real(8) :: e

            ! local variables
            real(8) :: s  =  0.449871
            real(8) :: t  = -0.386595
            real(8) :: a  =  0.19600
            real(8) :: b  =  0.25472
            real(8) :: r1 =  0.27597
            real(8) :: r2 =  0.27846
            real(8) :: u, v, x, y, q

            ! generate P = (u,v) uniform in rectangle enclosing acceptance region
            do
                call RANDOM_NUMBER(u)
                call RANDOM_NUMBER(v)
                v = 1.7156 * (v - 0.5)

                ! evaluate the quadratic form
                x = u - s
                y = abs(v) - t
                q = x**2 + y*(a*y - b*x)

                ! accept P if inside inner ellipse
                if (q < r1) exit
                ! reject P if outside outer ellipse
                if (q > r2) cycle
                ! reject P if outside acceptance region
                if (v**2 < -4.0*log(u)*u**2) exit
            end do

            ! return ratio of P's coordinates as the normal deviate
            e = v/u
            return

        end function random_normal


        ! init hermitian random matrix
        function hmat_init(dim) result(a)
            integer(4) :: dim
            integer(4) :: ii, jj
            complex(8), dimension(dim,dim) :: a


            ! pre-condition: check invalid dimensions for hmat
            call check(debug     = .TRUE., &
                       condition = (dim<=0), &
                       msg_type  = "Error", &
                       msg       = "Invalid dimensions at initialization of hermitian matrix", &
                       trg_stop  = .TRUE.)
            ! hmat generic random initialization
            do jj=1,dim
                do ii=1,dim
                    if (jj>ii) then
                        a(ii,jj) = COMPLEX(random_normal(), random_normal())
                        a(jj,ii) = CONJG(a(ii,jj))
                    else if (ii==jj) then
                        a(ii,jj) = COMPLEX(rand(), 0d0)
                    else
                        continue
                    end if
                end do
            end do
        end function hmat_init


        ! init hermitian random diagonal matrix
        function hmat_diag_init(dim) result(a)
            integer(4) :: dim
            integer(4) :: ii
            complex(8), dimension(dim,dim) :: a

            ! pre-condition: check invalid dimensions for hmat
            call check(debug     = .TRUE., &
                       condition = (dim<=0), &
                       msg_type  = "Error", &
                       msg       = "Invalid dimensions at initialization of hermitian matrix", &
                       trg_stop  = .TRUE.)
            ! hmat diagonal random initialization
            a = COMPLEX(0d0,0d0)
            do ii=1,dim
                a(ii,ii) = COMPLEX(random_normal(), 0d0)
            end do
        end function hmat_diag_init


        ! compute and return eigenvalues of hermitian matrix a
        function hmat_eigs(a) result(eigs)
            complex(8), dimension(:,:) :: a

            real(8), dimension(lbound(a,dim=1):ubound(a,dim=1)) :: eigs
            complex(8), dimension(:), allocatable :: work(:)
            real(8),    dimension(:), allocatable :: rwork(:)
            integer(4) :: N
            integer(4) :: lwork, info

            N     = size(a, 1)
            lwork = max(1,2*N-1)

            allocate(work(lwork))
            allocate(rwork(max(1, 3*N-2)))
            call zheev('N', 'U', N, a, N, eigs, work, lwork, rwork, info)
            deallocate(work)
            deallocate(rwork)
        end function hmat_eigs


        ! compute and return the normalized spacings between eigenvalues of a
        ! global mean version
        function hmat_s(eigs) result(s)
            real(8), dimension(:)  :: eigs

            real(8), dimension(lbound(eigs,dim=1):(ubound(eigs,dim=1)-1)) :: s
            real(8) :: s_m
            integer(4) :: ii

            do ii=1,size(eigs,1)-1
                s(ii) = eigs(ii+1) - eigs(ii)
            end do

            s_m = (eigs(size(eigs,1)) - eigs(1)) / size(s,1)

            do ii=1,size(s,1)
                s(ii) = s(ii) / s_m
            end do
        end function hmat_s


        ! compute and return the normalized spacings between eigenvalues of a
        ! local mean version
        function hmat_s_loc(eigs, level) result(s)
            real(8), dimension(:)  :: eigs

            integer(4) :: level
            real(8), dimension(lbound(eigs,dim=1):(ubound(eigs,dim=1)-1)) :: s
            integer(4) :: ii, ll, uu


            do ii=1,size(s,1)
                ll = MAX(1, ii-level)
                uu = MIN(size(s,1), ii+level)
                s(ii) = (eigs(ii+1) - eigs(ii)) / ((eigs(uu) - eigs(ll)) / (uu-ll+1))
            end do
        end function hmat_s_loc


        ! make histogram data
        subroutine hist2pdf(data, min, max, nbins, filename, unit)
            ! input arguments
            real(8), dimension(:) :: data
            real(8) :: min, max
            integer(4) :: nbins
            character(*) :: filename
            integer(4) :: unit

            real(8), dimension(nbins) :: binc, binh
            real(8) :: binw
            integer(4) :: ii, jj
            logical :: bin_inf, bin_sup, x_inf, x_sup

            ! set to zero bin centers and heights
            binc = 0
            binh = 0

            ! define bin centers
            binw = (max - min) / nbins
            do ii=1,size(binc)
                binc(ii) = min + (ii-0.5)*binw
            end do

            ! fill bins (~define bin heights)
            do ii=1,size(data,1)
                do jj=1,nbins
                    bin_inf = data(ii) > (min + (jj-1)*binw)
                    bin_sup = data(ii) < (min +  jj   *binw)
                    x_inf   = data(ii) > min
                    x_sup   = data(ii) < max
                    if (bin_inf .AND. bin_sup .AND. x_inf .AND. x_sup) then
                        binh(jj) = binh(jj) + 1
                    end if
                end do
            end do

            ! normalize bin heights to get a pdf
            do ii=1,size(binh,1)
                binh(ii) = binh(ii)/(size(data,1)*binw)
            end do

            ! write results on file to use for plotting
            open(unit, file=filename, status='replace')
            do ii=1,size(binc,1)
                write(unit,'(g0ag0)') binc(ii), char(9), binh(ii)
            end do
            close(unit)
        end subroutine hist2pdf


        ! find <r>
        subroutine hmat_r_avg(eigs, filename, unit)
            ! input arguments
            real(8), dimension(:) :: eigs
            character(*) :: filename
            integer(4) :: unit

            real(8), dimension(size(eigs,1)-1) :: s
            real(8), dimension(size(eigs,1)-1) :: r
            real(8) :: num, den, r_avg, r_avg_std
            integer(4) :: ii


            ! compute spacings
            s = 0
            do ii=1,size(eigs,1)-1
                s(ii) = eigs(ii+1) - eigs(ii)
            end do

            ! compute r
            do ii=1,size(s,1)-1
                num   = MIN(s(ii),s(ii+1))
                den   = MAX(s(ii),s(ii+1))
                r(ii) = num / den
            end do

            ! compute average of r
            r_avg = 0
            do ii=1,size(r,1)
                r_avg = r_avg + r(ii)
            end do
            r_avg = r_avg / size(r,1)

            ! compute std of average of r
            r_avg_std = 0
            do ii=1,size(r,1)
                r_avg_std = (r_avg - r(ii))**2
            end do
            r_avg_std = SQRT(r_avg_std / (size(r,1)-1)) / SQRT(1.0d0*size(r,1))

            ! write results on file
            open(unit, file=filename, status='replace')
            write(unit,'(g0ag0)') r_avg, char(9), r_avg_std
            close(unit)
        end subroutine hmat_r_avg


        ! print matrix on standard output (terminal)
        subroutine hmat_print_std(a, formatted)
            complex(8), dimension(:,:) :: a
            integer(4) :: ii, jj
            logical    :: formatted

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(a,1)
                if (formatted) then
                    write(*, "(*('('sf9.6xspf9.6'i)':x))") (a(ii,jj), jj=1,size(a,2))
                else
                    write(*, *) (a(ii,jj), jj=1,size(a,2))
                end if
            end do
        end subroutine hmat_print_std


        ! print matrix on standard output (terminal)
        subroutine hmat_print_file(a, fname, unit, formatted)
            complex(8), dimension(:,:) :: a
            integer(4) :: ii, jj
            integer(4) :: unit
            character(*) :: fname
            logical :: formatted

            open(unit=unit, file=fname, action='write', status='replace')

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            do ii=1,size(a,1)
                if (formatted) then
                    write(unit, "(*('('sf9.6xspf9.6'i)':x))") (a(ii,jj), jj=1,size(a,2))
                else
                    write(unit, *) (a(ii,jj), jj=1,size(a,2))
                end if
            end do

            close(unit)
        end subroutine hmat_print_file
end module hmat





program test
    use hmat
    implicit none


    complex(8), dimension(:,:), allocatable :: a
    real(8), dimension(:), allocatable :: eigvl
    real(8), dimension(:), allocatable :: ss_temp
    real(8), dimension(:), allocatable :: ss

    character(100) :: arg1char
    character(100) :: arg2char
    character(100) :: arg3char
    character(100) :: arg4char
    character(100) :: arg5char
    character(100) :: arg6char
    character(100) :: arg7char
    character(100) :: arg8char

    integer(4) :: arg2, arg3, arg4, arg7, arg8, samples, N, nbins, local, level
    real(8) :: arg5, arg6, xmin, xmax
    character(100) :: arg1, h_type

    integer(4) :: ii, ll, uu

    ! get command-line arguments
    CALL GET_COMMAND_ARGUMENT(1,arg1char)   ! 1 : matrix type (herm/diag)
    CALL GET_COMMAND_ARGUMENT(2,arg2char)   ! 2 : number of samples
    CALL GET_COMMAND_ARGUMENT(3,arg3char)   ! 3 : dimension of matrices N
    CALL GET_COMMAND_ARGUMENT(4,arg4char)   ! 4 : number of bins
    CALL GET_COMMAND_ARGUMENT(5,arg5char)   ! 5 : xmin for plot range
    CALL GET_COMMAND_ARGUMENT(6,arg6char)   ! 6 : xmax for plot range
    CALL GET_COMMAND_ARGUMENT(7,arg7char)   ! 7 : enable local average mode
    CALL GET_COMMAND_ARGUMENT(8,arg8char)   ! 8 : level for local average

    READ(arg1char, *) arg1
    READ(arg2char, *) arg2
    READ(arg3char, *) arg3
    READ(arg4char, *) arg4
    READ(arg5char, *) arg5
    READ(arg6char, *) arg6
    READ(arg7char, *) arg7
    READ(arg8char, *) arg8

    ! rename command-line arguments
    h_type  = arg1
    samples = arg2
    N       = arg3
    nbins   = arg4
    xmin    = arg5
    xmax    = arg6
    local   = arg7
    level   = arg8

    ! allocate memory
    allocate(a(N,N))
    allocate(ss((N-1)*samples))

    ! random hermitian matrix
    if (TRIM(h_type)=="herm") then
        ! local version
        if (local==1) then
            do ii=1,samples
                a       = hmat_init(N)
                eigvl   = hmat_eigs(a)
                ss_temp = hmat_s_loc(eigvl, level)

                ll = (ii-1)*(N-1) + 1
                uu = ii*(N-1)
                ss(ll:uu) = ss_temp
            end do

            call hist2pdf(ss, xmin, xmax, nbins, "results/herm_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//TRIM(arg8char)//"_level_"//"loc.dat", 1)
            call hmat_r_avg(eigvl,         "results/r_avg_herm_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//TRIM(arg8char)//"_level_"//"loc.dat", 1)
        ! global version
        else
            do ii=1,samples
                a       = hmat_init(N)
                eigvl   = hmat_eigs(a)
                ss_temp = hmat_s(eigvl)

                ll = (ii-1)*(N-1) + 1
                uu = ii*(N-1)
                ss(ll:uu) = ss_temp
            end do

            call hist2pdf(ss, xmin, xmax, nbins, "results/herm_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//"glob.dat", 1)
            call hmat_r_avg(eigvl,         "results/r_avg_herm_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//"glob.dat", 1)
        end if
    ! random diagonal real matrix
    else if (TRIM(h_type)=="diag") then
        ! local version
        if (local==1) then
            do ii=1,samples
                a       = hmat_diag_init(N)
                eigvl   = hmat_eigs(a)
                ss_temp = hmat_s_loc(eigvl, level)

                ll = (ii-1)*(N-1) + 1
                uu = ii*(N-1)
                ss(ll:uu) = ss_temp
            end do

            call hist2pdf(ss, xmin, xmax, nbins, "results/diag_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//TRIM(arg8char)//"_level_"//"loc.dat", 1)
            call hmat_r_avg(eigvl,         "results/r_avg_diag_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//TRIM(arg8char)//"_level_"//"loc.dat", 1)
        ! global version
        else
            do ii=1,samples
                a       = hmat_diag_init(N)
                eigvl   = hmat_eigs(a)
                ss_temp = hmat_s(eigvl)

                ll = (ii-1)*(N-1) + 1
                uu = ii*(N-1)
                ss(ll:uu) = ss_temp
            end do
            call hist2pdf(ss, xmin, xmax, nbins, "results/diag_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//"glob.dat", 1)
            call hmat_r_avg(eigvl,         "results/r_avg_diag_"//TRIM(arg2char)//"_samples_"//TRIM(arg3char)//"_N_"//"glob.dat", 1)
        end if
    else
        print *, "Invalid matrix initialization given"
    end if

end program test

! Ex02-Ardino-CODE.f90

module zmatmod
    implicit none

    ! define new type
    type zmat
        integer*2, dimension(2) :: dim
        complex*16, dimension(:,:), allocatable :: elem
        complex*16 :: tr
        complex*16 :: det
    end type zmat



    ! operator interfaces
    ! trace operator
    interface operator(.tr.)
        module procedure zmat_tr
    end interface
    ! adjoint operator
    interface operator(.adj.)
        module procedure zmat_adj
    end interface
    ! determinant operator (not yet implemented)
    interface operator(.det.)
        module procedure zmat_det
    end interface



    contains
        ! init zmat element with all zero entries
        function zmat_init_zero(dim) result(a)
            integer*2, dimension(2) :: dim
            integer*2 :: ii, jj
            type(zmat) :: a

            if ((dim(1)>0).AND.(dim(2)>0)) then
                a%dim(1) = dim(1)
                a%dim(2) = dim(2)
                allocate(a%elem(a%dim(1), a%dim(2)))

                do jj=1,dim(2)
                    do ii=1,dim(1)
                        a%elem(ii,jj) = COMPLEX(0d0, 0d0)
                    end do
                end do

                a%tr = zmat_tr(a)
                a%det = zmat_det(a)
            else
                print *, "Invalid dimensions at initialization of zmat type"
            end if
        end function zmat_init_zero


        ! init zmat element with random entries
        function zmat_init_rand(dim) result(a)
            integer*2, dimension(2) :: dim
            integer*2 :: ii, jj
            type(zmat) :: a

            if ((dim(1)>0).AND.(dim(2)>0)) then
                a%dim(1) = dim(1)
                a%dim(2) = dim(2)
                allocate(a%elem(a%dim(1), a%dim(2)))

                do jj=1,dim(2)
                    do ii=1,dim(1)
                        a%elem(ii,jj) = COMPLEX(rand(), rand())
                    end do
                end do

                a%tr = zmat_tr(a)
                a%det = zmat_det(a)
            else
                print *, "Invalid dimensions at initialization of zmat type"
            end if
        end function zmat_init_rand


        ! compute trace
        function zmat_tr(a) result(z_tr)
            integer*2 :: dd
            complex*16 :: z_tr
            type(zmat), intent(in) :: a

            z_tr = COMPLEX(0d0, 0d0)
            if (a%dim(1).EQ.a%dim(2)) then
                do dd=1,a%dim(1)
                    z_tr = z_tr + a%elem(dd,dd)
                end do
            else
                print *, "Cannot compute trace: zmat not a square matrix"
            end if
        end function zmat_tr


        ! TODO: compute determinant
        function zmat_det(a) result(z_det)
            complex*16 :: z_det
            type(zmat), intent(in) :: a

            z_det = COMPLEX(1d0, 0d0)
            if (a%dim(1).EQ.a%dim(2)) then
                ! TODO
            else
                print *, "Cannot compute determinant: zmat not a square matrix"
            end if
        end function zmat_det


        ! compute adjoint
        function zmat_adj(a) result(z_adj)
            type(zmat), intent(in) :: a
            type(zmat) :: z_adj

            z_adj%dim(1) = a%dim(2)
            z_adj%dim(2) = a%dim(1)

            allocate(z_adj%elem(z_adj%dim(1), z_adj%dim(2)))

            z_adj%elem = CONJG(TRANSPOSE(a%elem))

            if (a%dim(1).EQ.a%dim(2)) then
                z_adj%tr = CONJG(a%tr)
                z_adj%det = CONJG(a%det)
            else
                print *, "Cannot compute trace: zmat not a square matrix"
            end if
        end function zmat_adj


        ! print matrix on standard output (terminal)
        subroutine zmat_print_std(a, formatted)
            type(zmat) :: a
            integer*2 :: ii, jj
            logical :: formatted

            ! write matrix dimensions
            write(*, *) "Matrix dimensions:"
            write(*, *) a%dim(1), "x", a%dim(2)

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            write(*, *) "Matrix elements:"
            do ii=1,a%dim(1)
                if (formatted) then
                    write(*, "(*('('sf9.6xspf9.6'i)':x))") (a%elem(ii,jj), jj=1,a%dim(2))
                else
                    write(*, *) (a%elem(ii,jj), jj=1,a%dim(2))
                end if
            end do

            ! write matrix trace
            write(*, *) "Matrix trace:"
            write(*, *) a%tr

            ! write matrix determinant (=1 since it is not implemented)
            write(*, *) "Matrix determinant (not yet implemented):"
            write(*, *) a%det
        end subroutine zmat_print_std


        ! print matrix on standard output (terminal)
        subroutine zmat_print_file(a, fname, unit, formatted)
            type(zmat) :: a
            integer*2 :: ii, jj
            integer*4 :: unit
            character(*) :: fname
            logical :: formatted

            open(unit=unit, file=fname, action='write', status='replace')

            ! write matrix dimensions
            write(unit, *) "Matrix dimensions:"
            write(unit, *) a%dim(1), "x", a%dim(2)

            ! write matrix elements in a reduced format is formatted=.TRUE., otherwise in a raw way
            write(unit, *) "Matrix elements:"
            do ii=1,a%dim(1)
                if (formatted) then
                    write(unit, "(*('('sf9.6xspf9.6'i)':x))") (a%elem(ii,jj), jj=1,a%dim(2))
                else
                    write(unit, *) (a%elem(ii,jj), jj=1,a%dim(2))
                end if
            end do

            ! write matrix trace
            write(unit, *) "Matrix trace:"
            write(unit, *) a%tr

            ! write matrix determinant (=1 since it is not implemented)
            write(unit, *) "Matrix determinant (not yet implemented):"
            write(unit, *) a%det

            close(unit)
        end subroutine zmat_print_file
end module zmatmod





program test
    use zmatmod
    implicit none

    integer*2, dimension(2) :: dim
    type(zmat) :: a                                 ! declare zmat var of name 'a'
    type(zmat) :: a_adj                             ! declare zmat var for adj(a)

    dim(1) = 5                                      ! fix a%elem dimensions to (5,5)
    dim(2) = 5

    a = zmat_init_rand(dim)                         ! initialize randomly 'a'
    a_adj = .adj.a                                  ! adjoint of 'a'

    print *, "zmat variable 'a':"
    call zmat_print_std(a, .TRUE.)                  ! print zmat type on terminal
    print *, NEW_LINE('a'), "Adjoint of 'a':"
    call zmat_print_std(a_adj, .TRUE.)              ! print zmat adjoint on terminal

    print *, NEW_LINE('a'), "Writing 'a' variable on file..."
    call zmat_print_file(a, "zmat.txt", 20, .TRUE.) ! write zmat type on file

    print *, NEW_LINE('a'), "Check if new operators work fine..."
    print *, "Trace of a:           ", .tr.a
    print *, "Trace of adj(a):      ", .tr.a_adj
    print *, "Determinant of a:     ", .det.a       ! Not yet implemented
    print *, "Determinant of adj(a):", .det.a_adj   ! Not yet implemented
end program test

! 04.f90


! module err_handling_mod for error handling
module err_handling_mod
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
end module err_handling_mod





! module matrixmod for matmul functions
module matrixmod
    use err_handling_mod
    implicit none

    logical :: DEBUG
    logical :: TRIGGER_STOP

    contains
        ! matrix matrix multiplication by column
        function matmul_col(mat1, mat2) result(res)
            integer(4) :: ii, jj, kk
            integer(4) :: N, M, P, Q
            real(8), dimension(:,:) :: mat1, mat2
            real(8), dimension(lbound(mat1,dim=1):ubound(mat1,dim=1), &
                               lbound(mat2,dim=2):ubound(mat2,dim=2)) :: res


            ! mat1 and mat2 dimensions
            N = size(mat1, 1)
            M = size(mat1, 2)
            P = size(mat2, 1)
            Q = size(mat2, 2)

            ! check if the two input matrices have valid dimensions
            ! otherwise call the check routine
            ! check then if the dimensions are compatible for the multiplication
            ! pre-condition: check invalid dimensions for mat1
            call check(debug     = DEBUG, &
                       condition = (N<=0).OR.(M<=0), &
                       msg_type  = "Error", &
                       msg       = "Cannot multiply mat1*mat2: invalid dimensions for mat1", &
                       trg_stop  = TRIGGER_STOP)
            ! pre-condition: check invalid dimensions for mat2
            call check(debug     = DEBUG, &
                       condition = (P<=0).OR.(Q<=0), &
                       msg_type  = "Error", &
                       msg       = "Cannot multiply mat1*mat2: invalid dimensions for mat2", &
                       trg_stop  = TRIGGER_STOP)
            ! pre-condition: check if mat2 can multiply mat1 (mat1.cols=?=mat2.rows)
            call check(debug     = DEBUG, &
                       condition = M.NE.P, &
                       msg_type  = "Error", &
                       msg       = "Cannot do mat1*mat2: size(mat1,2) =/= size(mat2,1)", &
                       trg_stop  = TRIGGER_STOP)
            ! loop for matrix matrix multiplication
            do jj=1,Q
                do kk=1,M
                    do ii=1,N
                        res(ii,jj) = res(ii,jj) + mat1(ii,kk)*mat2(kk,jj)
                    end do
                end do
            end do
        end function


        ! matrix matrix multiplication by row
        function matmul_row(mat1, mat2) result(res)
            integer(4) :: ii, jj, kk
            integer(4) :: N, M, P, Q
            real(8), dimension(:,:) :: mat1, mat2
            real(8), dimension(lbound(mat1,dim=1):ubound(mat1,dim=1), &
                               lbound(mat2,dim=2):ubound(mat2,dim=2)) :: res


            ! mat1 and mat2 dimensions
            N = size(mat1, 1)
            M = size(mat1, 2)
            P = size(mat2, 1)
            Q = size(mat2, 2)

            ! check if the two input matrices have valid dimensions
            ! otherwise call the check routine
            ! check then if the dimensions are compatible for the multiplication
            ! pre-condition: check invalid dimensions for mat1
            call check(debug     = DEBUG, &
                       condition = (N<=0).OR.(M<=0), &
                       msg_type  = "Error", &
                       msg       = "Cannot multiply mat1*mat2: invalid dimensions for mat1", &
                       trg_stop  = TRIGGER_STOP)
            ! pre-condition: check invalid dimensions for mat2
            call check(debug     = DEBUG, &
                       condition = (P<=0).OR.(Q<=0), &
                       msg_type  = "Error", &
                       msg       = "Cannot multiply mat1*mat2: invalid dimensions for mat2", &
                       trg_stop  = TRIGGER_STOP)
            ! pre-condition: check if mat2 can multiply mat1 (mat1.cols=?=mat2.rows)
            call check(debug     = DEBUG, &
                       condition = M.NE.P, &
                       msg_type  = "Error", &
                       msg       = "Cannot do mat1*mat2: size(mat1,2) =/= size(mat2,1)", &
                       trg_stop  = TRIGGER_STOP)
            ! loop for matrix matrix multiplication
            do ii=1,N
                do kk=1,M
                    do jj=1,Q
                        res(ii,jj) = res(ii,jj) + mat1(ii,kk)*mat2(kk,jj)
                    end do
                end do
            end do
        end function

end module matrixmod





program performance
    use matrixmod
    implicit none


    ! square matrices dimension
    integer(4) :: N
    ! matrices dinamic declaration
    real(8), dimension(:,:), allocatable :: mat1, mat2, res1, res2, res3
    ! reals for tracking start and finish time of the function used
    real(8) :: start, finish
    real(8) :: start_1, finish_1
    real(8) :: start_2, finish_2
    real(8) :: start_3, finish_3


    ! get command-line arguments (square matrix size)
    character(len=:), allocatable :: argument
    integer*4 :: arglen
    call GET_COMMAND_ARGUMENT(1, length=arglen)
    allocate(character(arglen) :: argument)
    call GET_COMMAND_ARGUMENT(1, value=argument)

    ! convert command line argument to integer (matrix size)
    read(argument(:),'(i5)') N

    DEBUG = .TRUE.
    TRIGGER_STOP = .TRUE.


    ! Pre-condition: check if the first matrix is not too big in dimensions
    call check(debug     = DEBUG, &
               condition = (N>10000), &
               msg_type  = "Warning", &
               msg       = "matrices size too big", &
               trg_stop  = TRIGGER_STOP)


    ! allocate memory for the two matrices to multiply and for the results
    allocate(mat1(N,N))
    allocate(mat2(N,N))
    allocate(res1(N,N))
    allocate(res2(N,N))
    allocate(res3(N,N))


    ! fill input matrices with random numbers
    call RANDOM_NUMBER(mat1)
    call RANDOM_NUMBER(mat2)
    ! initialize output matrices to zero
    res1 = 0.0
    res2 = 0.0
    res3 = 0.0


    ! call matrix-matrix multiplication function and print the results
    ! print cpu time and write on file the result
    ! do it for matmul_col function...
    call CPU_TIME(start)
    res1 = matmul_col(mat1, mat2)
    call CPU_TIME(finish)
    print *, finish-start

    ! ...for matmul_row function...
    call CPU_TIME(start)
    res2 = matmul_row(mat1, mat2)
    call CPU_TIME(finish)
    print *, finish-start

    ! ...and for intrinsic matmul fortran function
    call CPU_TIME(start)
    res3 = MATMUL(mat1, mat2)
    call CPU_TIME(finish)
    print *, finish-start


    ! post-condition: check if res1 and res2 are equal
    call check(debug     = DEBUG, &
               condition = .NOT.eq_check(res1, res2, 1d-14), &
               msg_type  = "Warning", &
               msg       = "res1 and res2 are not equal at given precision", &
               trg_stop  = .FALSE.)
    ! post-condition: check if res1 and res3 are equal
    call check(debug     = DEBUG, &
               condition = .NOT.eq_check(res1, res3, 1d-14), &
               msg_type  = "Warning", &
               msg       = "res1 and res3 are not equal at given precision", &
               trg_stop  = .FALSE.)
    ! post-condition: check if res2 and res3 are equal
    call check(debug     = DEBUG, &
               condition = .NOT.eq_check(res2, res3, 1d-14), &
               msg_type  = "Warning", &
               msg       = "res2 and res3 are not equal at given precision", &
               trg_stop  = .FALSE.)

    ! deallocate memory before ending the program
    deallocate(mat1)
    deallocate(mat2)
    deallocate(res1)
    deallocate(res2)
    deallocate(res3)

end program

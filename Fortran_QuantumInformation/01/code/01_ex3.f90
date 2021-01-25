! 01_ex3.f90

! module matrixmod for matmul functions
module matrixmod
    implicit none

contains
    ! matrix matrix multiplication by column
    function matmul_col(mat1, mat2) result(mat3)
        integer*2 :: ii, jj, kk
        integer*2 :: N, M, L
        real*4, dimension(:,:) :: mat1, mat2
        real*4, dimension(size(mat1, 1),size(mat2, 2)) :: mat3

        N = size(mat1, 1)
        M = size(mat1, 2)
        L = size(mat2, 2)

        do jj=1,M
            do kk=1,L
                do ii=1,N
                    mat3(ii,jj) = mat3(ii,jj) + mat1(ii,kk)*mat2(kk,jj)
                end do
            end do
        end do
    end function

    ! matrix matrix multiplication by row
    function matmul_row(mat1, mat2) result(mat3)
        integer*2 :: ii, jj, kk
        integer*2 :: N, M, L
        real*4, dimension(:,:) :: mat1, mat2
        real*4, dimension(size(mat1, 1),size(mat2, 2)) :: mat3

        N = size(mat1, 1)
        M = size(mat1, 2)
        L = size(mat2, 2)

        do ii=1,N
            do kk=1,L
                do jj=1,M
                    mat3(ii,jj) = mat3(ii,jj) + mat1(ii,kk)*mat2(kk,jj)
                end do
            end do
        end do
    end function

end module matrixmod





program performance
    use matrixmod
    implicit none


    integer*2 :: nn
    integer*2 :: ii, jj
    real*4, dimension (:,:), allocatable :: mat1, mat2, res1, res2, res3
    real*4 :: start, finish

    ! get command-line arguments (matrix size)
    character(len=:), allocatable :: argument
    integer*4 :: arglen
    call GET_COMMAND_ARGUMENT(1, length=arglen)
    allocate(character(arglen) :: argument)
    call GET_COMMAND_ARGUMENT(1, value=argument)


    ! convert command line argument to integer (matrix size)
    read(argument(:),'(i5)') nn
    if (nn > 10000) then
        print *, "Matrix size too big"
        stop
    else if (nn==0) then
        nn = 500 ! default value
    end if

    allocate(mat1(nn,nn))
    allocate(mat2(nn,nn))
    allocate(res1(nn,nn))
    allocate(res2(nn,nn))
    allocate(res3(nn,nn))

    call RANDOM_NUMBER(mat1)
    call RANDOM_NUMBER(mat2)

    print *, "***************************"
    call CPU_TIME(start)
    res1 = matmul_col(mat1, mat2)
    call CPU_TIME(finish)
    open (unit=20, file="matmul_col_results.txt", position="append", action="write")
    write (20,*) nn, finish-start
    close (20)
    print *, "matmul_col time:", finish-start

    print *, "***************************"
    call CPU_TIME(start)
    res2 = matmul_row(mat1, mat2)
    call CPU_TIME(finish)
    open (unit=20, file="matmul_row_results.txt", position="append", action="write")
    write (20,*) nn, finish-start
    close (20)
    print *, "matmul_row time:", finish-start

    print *, "***************************"
    call CPU_TIME(start)
    res3 = MATMUL(mat1, mat2)
    call CPU_TIME(finish)
    open (unit=20, file="matmul_results.txt", position="append", action="write")
    write (20,*) nn, finish-start
    close (20)
    print *, "matmul time:", finish-start
    print *, "***************************"


    deallocate(mat1)
    deallocate(mat2)
    deallocate(res1)
    deallocate(res2)
    deallocate(res3)

end program

! 01_ex2.f90

program precision

    implicit none

    ! integer vars declaration
    integer*2 :: int_16_1, int_16_2, sum_int_16
    integer*4 :: int_32_1, int_32_2, sum_int_32
    ! real vars declaration
    real*4 :: real_32_1, real_32_2, sum_real_32
    real*8 :: real_64_1, real_64_2, sum_real_64


    ! integer sums
    int_16_1 = 2000000
    int_16_2 = 1
    int_32_1 = 2000000
    int_32_2 = 1

    sum_int_16 = int_16_1 + int_16_2
    sum_int_32 = int_32_1 + int_32_2

    print*, "INTEGER*2 sum:", sum_int_16
    print*, "INTEGER*4 sum:", sum_int_32


    ! real sums
    real_32_1 = 4.0e0*atan(1.0e0) * 1.0e32
    real_32_2 = sqrt(1.0e0) * 1.0e21
    real_64_1 = 4.0d0*atan(1.0d0) * 1.0d32
    real_64_2 = sqrt(1.0d0) * 1.0d21

    sum_real_32 = real_32_1 + real_32_2
    sum_real_64 = real_64_1 + real_64_2

    print*, "REAL*4 sum:", sum_real_32
    print*, "REAL*8 sum:", sum_real_64

end program

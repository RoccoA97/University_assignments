# Ex03: Documentation

## `err_handling_mod` module
Module used to handle errors and debugging inside the code of a program. Its core is composed of:
* **`check` subroutine**: it takes in input:
    * **`debug`** (logical): enable debugging output messages if enabled;
    * **`condition`** (logical): condition to verify in order to print debugging messages and/or stop the execution at checkpoint;
    * **`msg_type`** (string): flag of debugging message (i.e. Error/Warning/etc...);
    * **`msg`** (string): debugging message;
    * **`trg_stop`** (logical, optional): enable stop at checkpoint;
    * **`var`** (type given, optional): varible to print if debug is enabled.

* **`dim_check` function**: given two input matrices, it returns true if their dimensions are the same, otherwise it returns false.

* **`eq_check` function**: given two input matrices, it returns true if their entries are the same at a precision level encoded in the input threshold `thld`, otherwise it returns false.



## `matrixmod` module
Module containing several implementations of matrix matrix multiplication algorithm. It contains:
* **`matmul_col` function**: it multiplies the two input matrices (if possible) with the inner loop cycling over contiguous elements of the memory (optimal operation), and returns the result.

* **`matmul_row` function**: it multiplies the two input matrices (if possible) with the outer loop cycling over contiguous elements of the memory (non-optimal operation), and returns the result.



## `performance` program
It tests the cpu time needed for matrix matrix multiplication for every algorithm of `matrixmod` and for intrinsic Fortran `MATMUL` function. The size of the input matrices are given by the user at the beginning of the execution, so run-time. If `-debug` flag is given, several checks are done before the execution of the algorithms from `matrixmod` and after the execution in order to check the correctness of the results obtained.

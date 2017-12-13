// DiagonalizationRoutines.cc
#define lapack_complex_double std::complex<double> /* Workaround to get LAPACKE to take 
std::complex types...Redefines latter type by former type (?) Has to appear ***before*** 
lapacke.h. */
#include <lapacke.h> // Linear algebra routines (here: diagonalization)
#include <complex> // Complex numbers
#include <iostream> // Input from/output to command line
#include "diag_routines.h" // Include for consistency check with header

void ErrorMessage(const int info);

int simple_zheev(const int MatOrder, std::complex<double>*const a, double*const w, 
                 const bool OutputEvecs)
{
    /*
    https://software.intel.com/en-us/mkl-developer-reference-c-heev
    
    Finds the eigenvalues and optionally eigenvectors of a Hermitian matrix in full 
    storage layout. The first argument is the order (size) of the Hermitian matrix. The 
    second argument is the matrix in row-major full storage layout (i.e. a 
    MatOrder*MatOrder array) (double-check this, and the lower triangle thing). The evals 
    are written to the array w (of length MatOrder) in ascending order. If OutputEvecs is 
    true, the orthonormal evecs are written to a (column vectors? row vectors?).
    */
    /* We presume we are working with full arrays and not subarrays, so leading dimension 
    is the same as MatOrder. */
    const int lda=MatOrder;
    // Initialize info to unlikely value.
    int info=-99;
    // Choose the job type
    if (!OutputEvecs)
        char jobz = 'N';
    else
        char jobz = 'V';
    // The routine uses the lower triangle of the Hermitian matrix.
    const char uplo = 'L';
    
    info = LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, MatOrder, a, lda, w);
    
    // Give error message if 'info' is nonzero.
    ErrorMessage(info);
    return info;
}

void ErrorMessage(const int info)
{
    if (info!=0)
        std::cout << "**** Warning -- LAPACKE error\
                      \nSee documentation for error code meaning.\
                      \ninfo = " << info << std::endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// PALET.cpp: Implementation of the p-adic Lifting with Early Termination    //
// (PALET) algorithm.                                                        //
//                                                                           //
// (c) 2024 by Kelsey Kitzmiller, Texas A&M University. All rights reserved. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Headers
#include <iostream>
#include <givaro/modular.h>
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include <fflas-ffpack/fflas/fflas.h>
#include <givaro/modular-balanced.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <sys/resource.h>
#include <ctime>
#include <omp.h>

// Macros, etc.
#define Matrix_2D(A,i,j,n) (A[(i*n)+j])
typedef Givaro::Modular<int64_t> Field;
using namespace std;

// Function declarations
double get_wall_time();
mpq_t* ReadFile(char* fileName, int* n);
void ScaleSystem(mpq_t* rA, mpq_t* rb, int n, mpz_t* A, mpz_t* b);
int ComputeT(mpz_t* A, mpz_t* b, int n, int64_t pinit);
void ComputeDixonBound(int n, mpz_t sigma, mpz_t* B);
int PLUQfactorization(int n, Field F, size_t* P, int64_t* LU, size_t* Q);
void PALET(mpz_t* A, Field F, size_t* P, int64_t* LU, size_t* Q, mpz_t* b, mpz_t B, int n, int pinit, mpz_t* pfinal, int* iters, int* flag, int t, mpz_t* x, double* reconTime);
void DLCMRecon(mpz_t* IntSoln, int n, mpz_t prime, mpz_t DixBound, int* flag, mpq_t* soln);
void RRNumCalc(mpz_t num, mpz_t x, mpz_t p);
int CheckNumMatrix(mpz_t* numMatrix, int i, int t);
int CheckFlags(int* vectorFlags, int t);
void FreeMPZArray(mpz_t* x, int n);
void FreeMPQArray(mpq_t* x, int n);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ALGORITHM DESCRIPTION: Implementation of the PALET algorithm. GMP (6.2.1) and fflas-ffpack (2.5.0) are used to handle integers of arbitrary
// size, and the finite field factorization, respectively.

int main(int argc, char* argv[])
{
    //*********************************************************** Declare all vars ***********************************************************//
    // Input vars
    int n;                     // The dimension of the system to be solved
    mpq_t* rA;                 // The rational matrix A
    mpq_t* rb;                 // The rational RHS vector b

    // Vars
    mpz_t* A;                  // The integral (scaled) matrix A
    mpz_t* b;                  // The integral (scaled) RHS vector b
    mpz_t sigma;               // The largest entry (absolute value) of the scaled system
    int64_t pinit;             // The prime modulus p
    int t;                     // Parameter t (for New Dixon Alg only) to be computed
    mpz_t B;                   // Dixon's bound

    // Timer vars
    double wallTimer=0;        // Used to compute total time spent in PALET alg (exclusive of pre/post-processing steps)
    double wallFactorTime=0;   // Used to compute total time spent in  factorization step
    double wallLiftTime=0;     // Used to compute total time spent in lifting loop (exclusive of numerator reconstructions)
    double wallNumReconTime=0; // Used to compute total time spent in numerator reconstruction steps
    double wallDLCMtime=0;     // Used to compute total time spent in final DLCM reconstruction

    // Factorization vars
    size_t* P;                 // Permutation P used in PLUQ factorization
    int64_t* LU;               // Matrix L*U used in PLUQ factorization
    size_t* Q;                 // Permutation Q used in PLUQ factorization

    // Solution vars
    int iters;                 // The current p-adic lifting iteration
    mpz_t* x;	               // The current p-adic lifting modular solution
    mpz_t pfinal;              // The final value of p^k when lifting loop is terminated
    mpq_t* soln;               // The final rational solution to the system Ax=b

    // Flags
    int flag;                  // Indicate if early termination was achieved (1) or not (0)
    int factorCheck;           // Indicate if factorization was successful (1) or not (0) 
    int reconCheck;            // Indicate if final reconstruction is successful (1) or not (0)

    // GMP Initializations
    mpz_init(sigma);           // Initialize sigma=0
    mpz_init(B);               // Initialize B=0
    mpz_init(pfinal);          // Initialize pfinal=0


    //************************************************************* Check input **************************************************************//
    if(argc!=2)             // Too many or too few arguments. Print error and terminate program
    {
        (argc>2) ? printf("ERROR: Too many arguments. Usage: ./PALET [Matrixfile.txt]\n") : printf("ERROR: Too few arguments. Usage: ./PALET [Matrixfile.txt]\n");
        return 0;
    }


    //************************************************************** Read file ***************************************************************//
    rA = ReadFile(argv[1], &n);

    // Check if read was successful
    if(rA==NULL)
    {
        // Read failed, print error and terminate program
        printf("Read failed\n");
        return 0;
    }


    //******************************************************** Set RHS to unit vector ********************************************************//
    rb = (mpq_t*) calloc(n, sizeof(mpq_t));     // Initialize vector rb
    mpq_set_ui(rb[0], 1, 1);                    // Set first entry to 1
    for(int i=1; i<n; i++)                      // Set all other entries to 0
    {
        mpq_set_ui(rb[i], 0, 1);
    }


    //***************************************************** Scale A and b, compute sigma *****************************************************//
    A = (mpz_t*) calloc(n*n, sizeof(mpz_t));    // Initialize matrix A
    b = (mpz_t*) calloc(n, sizeof(mpz_t));      // Initialize vector b
    
    ScaleSystem(rA, rb, n, A, b);               // Scale rA and rb by appropriate constant, return integral A, b

    // Compute sigma
    for(int i=0; i<n; i++)
    {
        // Compare current sigma value to entry of b. Replace sigma if b[i]>sigma.
        if(mpz_cmp(sigma, b[i])<0)
    {
            mpz_set(sigma, b[i]);
        }

        // Compare current sigma value to each entry in the ith row of A. Replace sigma if A[i,j]>sigma.
        for(int j=0; j<n; j++)
        {
            if(mpz_cmp(sigma, Matrix_2D(A,i,j,n))<0)
            {
                mpz_set(sigma, Matrix_2D(A,i,j,n));
            }
        }
    }

    x = (mpz_t*) calloc(n, sizeof(mpz_t));      // Initialize vector x
    soln = (mpq_t*) calloc(n, sizeof(mpq_t));   // Initialize solution vector 


    //**************************************************** Choose initial prime, compute t ***************************************************//
    pinit = 1000000007;                         // Set pinit to chosen prime

    // Create field
    Field F(pinit);

    // Compute t
    t = ComputeT(A, b, n, pinit);


    //******************************************************** Compute Dixon's bound *********************************************************//
    ComputeDixonBound(n, sigma, &B);            // Set B=2(sigma^2n)(n^n)


    //**************************************************** Initialize factorization vars *****************************************************//
    // Copy A mod pinit  into LU
    LU = FFLAS::fflas_new<int64_t>(n*n);        // Initialize LU matrix
    mpz_t temp;                                 // Declare temp variable
    mpz_init(temp);                             // Initialize temp variable

    for(int i=0; i<n*n; i++)
    {
        // Set each entry LU[i] = A[i] mod pinit
        mpz_mod_ui(temp, A[i], (unsigned int)pinit);
        F.init(LU[i], (int64_t)mpz_get_ui(temp));
    }
    mpz_clear(temp);

    // Initialize permutation vectors
    P = FFLAS::fflas_new<size_t>(n);
    Q = FFLAS::fflas_new<size_t>(n);


    //************************************************************* Start timer **************************************************************//
    wallTimer = get_wall_time();


    //*************************************************************** Factor A ***************************************************************//
    // Start factorization timer
    wallFactorTime = get_wall_time();

    // Factor A
    factorCheck = PLUQfactorization(n, F, P, LU, Q);

    // Check if factorization was successful
    if(factorCheck==0)
    {
        // Factorization failed, print error and terminate program
        printf("Factorization failed\n");
        return 0;
    }

    // End factorization timer
    wallFactorTime = get_wall_time() - wallFactorTime;


    //******************************************************** Perform p-adic lifting ********************************************************//
    // Start lifting timer
    wallLiftTime = get_wall_time();

    // Initialize iteration count and flag to 0
    iters = 0;
    flag = 0; // Assume early termination is not satisfied (0), otherwise it will be set to 1 by PALET function

    // Perform p-adic lifting with PALET criteria
    PALET(A, F, P, LU, Q, b, B, n, pinit, &pfinal, &iters, &flag, t, x, &wallNumReconTime);

    // End lifting timer
    wallLiftTime = (get_wall_time() - wallLiftTime) - wallNumReconTime; // Lift time (exclusive of time spent on numerator reconstructions)


    //***************************************************** Perform DLCM reconstruction ******************************************************//
    // Assume reconstruction is successful (1), otherwise it will be set to 0 by DLCM function
    reconCheck = 1;

    // Start reconstruction timer
    wallDLCMtime = get_wall_time();

    // Perform reconstruction
    DLCMRecon(x, n, pfinal, B, &reconCheck, soln);

    // End reconstruction timer
    wallDLCMtime = get_wall_time() - wallDLCMtime;

    // Check if reconstruction was successful
    if(reconCheck==0)
    {
        // Reconstruction failed, print error and terminate program
        printf("ERROR: Final reconstruction failed\n");
        return 0;
    }


    //************************************************************** End timer ***************************************************************//
    wallTimer = get_wall_time() - wallTimer;


    //************************************************************ Print results *************************************************************//
    printf("\n****************   Results   ****************\n");
    printf("Matrix File: \t%s\n", argv[1]);
    printf("Dimension: \t%d\n", n);
    printf("Method:\t\tPALET Algorithm\n");
    printf("Modulus: \t%ld\n", pinit);
    printf("t value:\t%d\n", t);
    printf("---------------------------------------------\n");
    (flag==1) ? printf("Early termination criteria satisfied\n") : printf("Dixon's bound achieved\n");
    printf("Solution found at iteration:\t\t %d\n", iters);
    printf("Total time (seconds): \t\t%.2f      100.00%%\n", wallTimer);
    printf("  Factorization: \t\t%.2f       %.2f%%\n", wallFactorTime, 100*(wallFactorTime / wallTimer));
    printf("  Lifting: \t\t\t%.2f       %.2f%%\n", wallLiftTime, 100*(wallLiftTime / wallTimer));
    printf("  Numerator Reconstructions:\t%.2f     %.2f%%\n", wallNumReconTime, 100*(wallNumReconTime / wallTimer));
    printf("  DLCM Reconstruction:\t\t%.2f       %.2f%%\n", wallDLCMtime, 100*(wallDLCMtime / wallTimer));
    printf("*********************************************\n\n");

    // Free vars
    FreeMPQArray(rA, n*n);
    FreeMPQArray(rb, n);
    FreeMPZArray(A, n*n);
    FreeMPZArray(b, n);
    mpz_clear(sigma);
    mpz_clear(B);
    FreeMPZArray(x, n);
    FreeMPQArray(soln, n);
    mpz_clear(pfinal);
    free(P);
    free(LU);
    free(Q);

    // Return success
    return 1;
}


// Description: This function returns the current system time as a double. If error occurs, return 0.
double get_wall_time()
{
    // Declare vars
    struct timeval time; // Used to hold current system time

    // Get system time
    if (gettimeofday(&time,NULL))
    {
        // Report error and return 0.
        printf("ERROR GETTING SYSTEM TIME\n");
        return 0;
    }

    // Return time as double
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


// Description: Read given matrix file, assign n and return the matrix
mpq_t* ReadFile(
    char* fileName, // The file to read. Not modified by function.
    int* n          // The dimension of the given matrix. This will be modified by function.
    )
{
    // Declare vars
    int flag;       // Indicate if the read was successful (1) or not (0)
    int dim1;       // The number of rows of the matrix
    int dim2;       // The number of columns of the matrix
    FILE* file;     // The file name
    mpq_t* A;       // The matrix to read into

    // Open file
    file = fopen(fileName, "r");
    if (file == NULL)
    {
        // Error opening file, print error and return failure
        printf("ERROR: File not found\n");
        return NULL;
    }

    // Read matrix dimensions
    flag = fscanf(file, "%d %d", &dim1, &dim2); // First two entries of file should be row and col dimension of A
    if(flag == 0)
    {
        // The file is empty, print error and return failure
        printf("ERROR: Empty file\n");
        return NULL;
    }
    else if(dim1!=dim2 || dim1<1)
    {
        // The matrix is not square, print error and return failure
        printf("ERROR: Invalid matrix dimensions\n");
        return NULL;
    }

    // Initialize A matrix
    A = (mpq_t*) calloc(dim1*dim2, sizeof(mpq_t));

    // Read matrix entries
    for(int i=0; i<dim1*dim2; i++)
    {
        mpq_inp_str(A[i], file, 10);  // Read rational number from character string
        mpq_canonicalize(A[i]);       // Simplify rational number
    }

    // Close file
    fclose(file);

    // Report dimension and return the matrix
    *n = dim1;
    return A;

    // Free vars
    FreeMPQArray(A, dim1*dim2);

}


// Description: This function takes a rational matrix rA and vector rb and computes the lcm of the denomiators of each entry.
// The system is scaled by the lcm value, and the integral matrix A and vector b are computed.
void ScaleSystem(
    mpq_t* rA,   // The rational matrix to be scaled. Not modified by function.
    mpq_t* rb,   // The rational vector to be scaled. Not modified by function.
    int n,       // The number of rows/cols of A and the length of b. Not modified by function.
    mpz_t* A,    // The scaled matrix. This will be modified by function.
    mpz_t* b     // The scaled vector. This will be modified by function.
    )
{
    // Declare vars
    mpz_t temp;  // The temporary variable used to hold the denominator of each entry of rA and rb
    mpz_t LCM;   // The lcm of the denominator of rA and rb

    // Initialize vars to 0
    mpz_init(temp);
    mpz_init(LCM);

    // Set lcm to 1
    mpz_set_ui(LCM, 1);

    // Compute lcm of denominators of rA and rb row-wise
    for(int i=0; i<n; i++)
    {
        mpz_set(temp, mpq_denref(rb[i])); // Set temp to the denominator of rb[i]
        for(int j=0; j<n; j++)
        {
            // For each entry in the ith row of rA, compute lcm of temp and the entry
            mpz_lcm(temp, temp, mpq_denref(Matrix_2D(rA,i,j,n)));
        }

        // Update lcm with the lcm of the ith row of rA, rb
        mpz_lcm(LCM, LCM, temp);
    }

    // Scale system
    for(int i=0; i<n; i++)
    {
        // Scale the ith entry of rb
        mpz_divexact(temp, LCM, mpq_denref(rb[i]));          // Divide lcm by denominator of rb[i], save as temp
        mpz_mul(b[i], mpq_numref(rb[i]), temp);              // Multiply numerator of rb[i] by temp
        
        // Scale the ith row of rA
        for(int j=0; j<n; j++)
        {
            mpz_divexact(temp, LCM, mpq_denref(rA[i*n+j]));  // Divide lcm by denominator of rA{i,j], save as temp
            mpz_mul(A[i*n+j], mpq_numref(rA[i*n+j]), temp);  // Multiply numerator of rA[i,j] by temp
        }
    }

    // Free vars
    mpz_clear(temp);
    mpz_clear(LCM);
}


// Description: This function computes and returns the parameter t required for PALET criteria
int ComputeT(
    mpz_t* A,       // The (integral) matrix of the system to solve. Not modified by function.
    mpz_t* b,       // The (integral) RHS of the system to solve. Not modified by function.
    int n,          // The dimension of the system. Not modified by function.
    int64_t pinit   // The chosen modulus. Not modified by function.
    )
{
    // Declare vars
    mpz_t max;      // This will hold the max(sum|aij|, |bi|) computed up to the current row i
    mpz_t a;        // This will hold the absolute value of aij
    mpz_t sum;      // This will hold the sum of |aij| for the current row i
    mpz_t rhs;      // This will hold |bi| for the current row i
    mpz_t C;        // This will hold the threshold for p^t
    int t=1;        // The value of t
    mpz_t power;    // Temporary variable to hold value of p^t

    // Initialize vars to 0
    mpz_init(max);
    mpz_init(a);
    mpz_init(sum);
    mpz_init(rhs);
    mpz_init(C);
    mpz_init(power);

    // Compute max(sum(|aij|), |bi|)
    for(int i=0; i<n; i++)
    {
        // Compute sum of row i
        mpz_set_ui(sum, 0);
        for(int j=0; j<n; j++)
        {
            mpz_set(a, Matrix_2D(A, i, j, n));    // Get the entry A[i,j]
            mpz_abs(a, a);                        // Take the absolute value of A[i,j]
            mpz_add(sum, sum, a);                 // Add the absolute value to the current sum
        }

        // Compute |bi|
        mpz_set(rhs, b[i]);
        mpz_abs(rhs, rhs);

        // Find max of computed Row-sum, RHS value, and the current maximum
        if(mpz_cmp(sum, rhs)>=0)   // The Row-sum of A is >= RHS value. Compare Row-sum with current maximum
        {
            if(mpz_cmp(sum, max)>=0)
            {
                mpz_set(max, sum); // The row-sum is >= the current max. Replace max with row-sum
            }
        }
        else                       // RHS value is > Row-sum of A. Compare RHS with current maximum.
        {
            if(mpz_cmp(rhs,max)>=0)
            {
                mpz_set(max, rhs); // The RHS is >= the current max. Replace max with RHS value
            }
        }
    }

    // Multiply max by sqrt(2), save as C (NOTE: Cannot do C = max * sqrt(2) becuase mpz_sqrt will truncate sqrt(2) to 1)
    mpz_pow_ui(C, max, 2);// Set C = max^2
    mpz_mul_ui(C, C, 2);  // Set C = 2*max^2
    mpz_sqrt(C, C);       // Set C = floor(sqrt(2*max^2))
    mpz_add_ui(C, C, 1);  // Add 1 to get C = ceil(sqrt(2*max^2))

    // If C is greater than pinit, compute t>1
    if(mpz_cmp_ui(C, (unsigned int)pinit)>0)
    {
        // Set power = p^1
        mpz_set_ui(power, (unsigned int)pinit);

        // Compute t
        while(mpz_cmp(C, power)>0)   // Increment t until C<=p^t
        {
            t++;
            mpz_ui_pow_ui(power, (unsigned int)pinit, (unsigned int)t); // Compute new value of p^t
        }
    }

    // Return t
    return t;

    // Free vars
    mpz_clear(max);
    mpz_clear(a);
    mpz_clear(sum);
    mpz_clear(rhs);
    mpz_clear(C);
    mpz_clear(power);
}



// Description: This function computes Dixon's bound for the given system (2 * sigma^(2n) * n^n)
void ComputeDixonBound(
    int n,        // The dimension (# of rows/columns) of the square matrix A. Not modified by function.
    mpz_t sigma,  // The bound on the size of entries in A, b. Not modified by function.
    mpz_t* B      // This will hold the value of Dixon's bound. This will be modified by function.
    )
{
    // Declare and initialize vars
    mpz_t gmpn;         // This will hold n as a gmp integer
    mpz_t temp1;        // Temporary variable to hold n^n
    unsigned int temp2; // Temporary variable to hold 2*n
    mpz_t temp3;        // Temporary variable to hold sigma^(2n)
    mpz_t result;       // This will hold the computed value of Dixon's bound

    // Initialize vars to 0
    mpz_init(gmpn);
    mpz_init(temp1);
    mpz_init(temp3);
    mpz_init(result);

    // Set gmpn
    mpz_set_ui(gmpn, (unsigned int)n);

    // Compute n^n
    mpz_pow_ui(temp1, gmpn, (unsigned int)n);

    // Compute sigma^(2n)
    temp2 = (unsigned int) 2*n;         // Compute 2n
    mpz_pow_ui(temp3, sigma, temp2);    // Raise sigma to 2n power

    // Compute Dixon's bound
    mpz_mul(result, temp1, temp3);      // Result = sigma^(2n) * n^n
    mpz_mul_ui(result, result, 2);      // Result = 2 * sigma^(2n) * n^n

    // Record result
    mpz_set(*B, result);

    // Free all vars
    mpz_clear(gmpn);
    mpz_clear(temp1);
    mpz_clear(temp3);
    mpz_clear(result);
}


// Description: This function computes the finite field PLUQ factorization of a given matrix using FFPACK. Success (1) or failure (0) is returned.
int PLUQfactorization(
    int n,        // The dimension of the matrix to be factored. Not modified by function.
    Field F,      // The finitie field. Not modified by function.
    size_t * P,   // The permutation P of the PLUQ factorization. This will be modified by function.
    int64_t* LU,  // The matrix L*U of the PLUQ factorization. Initially, this is the matrix to factor. This will be modified by function.
    size_t * Q    // The permutation Q of the PLUQ factorization. This will be modified by function.
    )
{
    // Declare vars
    size_t rank;  // This will hold the rank of the factored matrix
    
    // Factor A
    rank = FFPACK::PLUQ (F, FFLAS::FflasNonUnit, n, n, LU, n, P, Q);
    
    // Check if factorization was successful
    if((int)rank<n)
    {
        // Factorization failed. Print error and return failure
        printf("ERROR: Matrix is rank-deficient\n");
        return 0;
    }

    // Return success
    return 1;
}



// Description: This function performs p-adic lifting using the PALET criteria.
void PALET(
    mpz_t* A,              // The integral matrix of the system to be solved. Not modified by function.
    Field F,               // The finite field to use. Not modified by function.
    size_t* P,             // The permutation P of PLUQ factorization. Not modified by function.
    int64_t* LU,           // The L*U matrix of PLUQ factorization. Not modified by function.
    size_t* Q,             // The permutation Q of PLUQ factorization. Not modified by function.
    mpz_t* b,              // The integral RHS of the system to be solved. Not modified by function.
    mpz_t B,               // Dixon's bound for the system to be solved. Not modified by function.
    int n,                 // Dimension of system to be solved. Not modified by function.
    int pinit,             // The chosen modulus. Not modified by function.
    mpz_t* pfinal,         // The final value of p^k when lifting is terminated. This will be modified by the function.
    int* iters,            // The total number of p-adic lifting iterations performed. This will be modified by the function.
    int* flag,             // Indicates if early termination criteria is satisfied (1) or not (0). This will be modified by function.
    int t,                 // The parameter t reuqired for PALET criteria. Not modified by function.
    mpz_t* x,              // The integral solution vector. This will be modified by function.
    double* reconTime      // The time spent performing numerator reconstructions. This will be modified by function.
    )
{
    // Declare vars
    mpz_t* d;              // The vector d, which will equal (b-Ax_i)/p^(i+1) at each iteration i
    int64_t* y;            // The vector y, which will equal (A^(-1)*d) MOD p at each iteration
    int64_t* newy;         // Temporary vector to hold the values of y at the next iteration
    mpz_t* mpztemp;        // Temporary variable used to compute x_i
    mpz_t p;               // This will hold the value of p^i at each iteration i
    int iter;              // This will hold the iteration number i
    mpz_t RRp;             // This will hold the value of p^(i+1) for numerator reconstructions
    int RRindex;           // The entry of xhat to reconstruct initially (first phase of practical approach). This is chosen below.
    mpz_t* numVector;      // Vector to hold the new reconstructed numerator values for each entry of xhat computed at the current iteration.
    mpz_t* numMatrix;      // Matrix to hold the reconstructed numerator values for each entry of xhat, over the past t+1 iterations.
    int* vectorFlags;      // Vector of indicator variables for each entry of xhat to signal if entry has converged (1) or not (0).
    double tempWallTime=0; // Temporary variable to compute time spent in numerator reconstructions.

    // Initialize vars
    d = (mpz_t*) calloc(n, sizeof(mpz_t));
    y = (int64_t*) calloc(n, sizeof(int64_t));
    newy = (int64_t*) calloc(n, sizeof(int64_t));
    numVector = (mpz_t*) calloc(n, sizeof(mpz_t));
    numMatrix = (mpz_t*) calloc((t+1)*n, sizeof(mpz_t));
    vectorFlags = (int*) calloc(n, sizeof(int));
    mpztemp = (mpz_t*) calloc(n, sizeof(mpz_t));
    mpz_init(p);
    mpz_init(RRp);

    // Set initial values
    iter = 0;
    *flag = 0;
    RRindex = 0;           // First entry is chosen to reconstruct initially
    mpz_set_ui(p, 1);      // First iteration sets p = p^i = p^0 = 1

    // Set initial values for d and y vectors
    for(int i=0; i<n; i++)
    {
        // Set d = b
        mpz_set(d[i], b[i]);

        // Initialize mpztemp
        mpz_set_ui(mpztemp[i], 0);

        // Set y = b mod pinit for initial backsolve (see note below)
        mpz_mod_ui(mpztemp[i], b[i], (unsigned int)pinit);
        y[i] = (int64_t) mpz_get_ui(mpztemp[i]); // Mod maps to positive number [0,...,pinit) where pinit < 2^64, so this should not overflow
    }

    // Declare/initialize trivial finite field for backsolve
    Field::Element one;
    F.init (one, 1);

    // Lifting loop
    while(mpz_cmp(p, B)<0 && *flag==0)   // Iterate until p>=B or early termination criteria satisfied
    {
        // Backsolve: Set y = A^(-1)d MOD pinit (NOTE: fflas-ffpack functions overwrite RHS vector with solution so y must be equal to d before calling functions)
        FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, y, 1, P);                                             // Apply P
        FFLAS::ftrsm(F, FFLAS::FflasLeft , FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, n, 1, one, LU, n, y, 1);   // Apply L
        FFLAS::ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, n, 1, one, LU, n, y, 1); // Apply U
        FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, y, 1, Q);                                            // Apply Q
        // y is now A^(-1)d MOD pinit

        // Lifting step: Update x, d, and y
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for(int i=0; i<n; i++)
            {
                // Set x = x + yp
                mpz_mul_ui(mpztemp[i], p, (unsigned int)y[i]);
                mpz_add(x[i], x[i], mpztemp[i]);

                // Set d = (d - Ay)/pinit
                for(int j=0; j<n; j++)
                {
                    mpz_submul_ui(d[i], Matrix_2D(A,i,j,n), (unsigned int)y[j]);
                }
                mpz_divexact_ui(d[i], d[i], (unsigned int)pinit);

                // Set newy = d mod pinit. Must be saved as a new vector (newy) because each entry of y is used in step above. Cannot overwrite y[i] yet.
                mpz_mod_ui(mpztemp[i], d[i], (unsigned int)pinit);
                newy[i] = (int64_t) mpz_get_ui(mpztemp[i]); // Mod maps to positive number [0,...,pinit) where pinit < 2^64, so this should not overflow
            }
        }

        // Replace y with newy. This sets y to the RHS for next iteration.
        for(int i=0; i<n; i++)
        {
            y[i]=newy[i];
        }
        // y is now d MOD pinit

        // Perform numerator reconstruction on each entry of x
        {
            // Start reconstruction timer
            tempWallTime = get_wall_time();

            // Compute reconstruction modulus p^(i+1)
            mpz_mul_ui(RRp, p, (unsigned int)pinit);

            if(vectorFlags[RRindex]!=1) // First entry HAS NOT converged
            {
                RRNumCalc(numVector[RRindex], x[RRindex], RRp);                                    // Reconstruct numerator of first entry
                mpz_set(Matrix_2D(numMatrix, RRindex, ((iter)%(t+1)), (t+1)), numVector[RRindex]); // Add to list of reconstructed numerators
                vectorFlags[RRindex] = CheckNumMatrix(numMatrix, RRindex, t);                      // Check the list of numerators to see if entry has converged
            }
            else // First entry HAS converged, reconstruct the rest of the entries
            {
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic)
                    for(int i=0; i<n; i++)
                    {
                        // Check if reconstruction is necessary for entry i
                        if(vectorFlags[i]!=1) // The entry has not converged, perform reconstruction
                        {
                            RRNumCalc(numVector[i], x[i], RRp);                                    // Reconstruct numerator of entry i
                            #pragma omp critical
                            mpz_set(Matrix_2D(numMatrix, i, ((iter)%(t+1)), (t+1)), numVector[i]); // Add to list of reconstructed numerators
                            #pragma omp critical
                            vectorFlags[i] = CheckNumMatrix(numMatrix, i, t);                      // Check list of numerators to see if entry i has converged
                        }
                    }
                }
                // Check if all entries have converged (if so, flag will be set to 1)
                *flag = CheckFlags(vectorFlags, n);
            }
            // End reconstruction timer
            *reconTime += get_wall_time() - tempWallTime;
        }

        // Increment p^i and iter
        iter++;
        mpz_mul_ui(p, p, (unsigned int)pinit);
    }

    // Save pfinal and total iterations performed
    mpz_set(*pfinal, p);
    *iters = iter;

    // Free vars
    free(y);
    free(newy);
    FreeMPZArray(d, n);
    FreeMPZArray(numVector, n);
    FreeMPZArray(numMatrix, n*(t+1));
    free(vectorFlags);
    FreeMPZArray(mpztemp, n);
    mpz_clear(p);
    mpz_clear(RRp);
}



// Description: This function computes performs DLCM Reconstruction on an integral vector (modulo a given prime) and saves the result as soln.
void DLCMRecon(
    mpz_t* IntSoln, // The integral vector to be reconstructed. Not modified by function.
    int n,          // The length of the vector to be reconstructed. Not modified by function.
    mpz_t prime,    // The given prime modulus. Not modified by function.
    mpz_t DixBound, // Dixon's bound on the system the reconstructed vector should solve. Not modified by function.
    int* flag,      // Flag indicating if reconstruction is successful (1) or not (0). This will be modified by function.
    mpq_t* soln     // The reconstructed vector. This will be modified by function.
    )
{
    // Declare/Initialize vars
    mpz_t LCM;      // This will hold the LCM of the computed denominators
    mpz_t a;        // This will hold the value of (IntSoln[i]*lcm) MOD prime
    mpz_t b;        // The RR bound sqrt(prime/2) on the numerator value w0.
    mpz_t u;        // The bound prime-b (used to avoid RR iterations)
    mpz_t v0;       // First entry of vector v (initially = prime)
    mpz_t v1;       // Second entry of vector v (initially = 0)
    mpz_t w0;       // First entry of vector w (initially = a)
    mpz_t w1;       // Second entry of vector w (initially = 1)
    mpz_t q;        // The quotient q
    mpz_t z0;       // First entry of vector z
    mpz_t z1;       // Second entry of vector z
    mpz_t HadBound; // This will hold the value of Hadamard's bound for the system solved.
    mpz_t d;        // The RR bound on denominator value w1 (b divided by lcm). When w1 is multiplied by lcm, the reconstructed numerator will satisfy b.
    mpz_t temp;     // Temporary variable used to update lcm

    // Initialize vars to 0
    mpz_init(LCM);
    mpz_init(a);
    mpz_init(b);
    mpz_init(u);
    mpz_init(v0);
    mpz_init(v1);
    mpz_init(w0);
    mpz_init(w1);
    mpz_init(q);
    mpz_init(z0);
    mpz_init(z1);
    mpz_init(HadBound);
    mpz_init(d);
    mpz_init(temp);

    // Compute the Rational Reconstruction bound b
    mpz_cdiv_q_2exp(b, prime, 1);           // Sets b = ceil(p/2)
    mpz_sqrt(b, b);                         // Sets b = floor(sqrt(...)) but we need ceil
    mpz_add_ui(b, b, 1);                    // This gives b = ceil(sqrt(...))

    // Compute u = prime - b
    mpz_sub(u, prime, b);

    // Compute Hadamard's bound
    mpz_cdiv_q_2exp(HadBound, DixBound, 1); // Sets HadBound = ceil(DixBound/2)
    mpz_sqrt(HadBound, HadBound);           // Sets HadBound = floor(sqrt(...)) but we need ceil
    mpz_add_ui(HadBound, HadBound, 1);      // This gives HadBound = ceil(sqrt(...))

    // Initialize lcm and flag to 1
    mpz_set_si(LCM, 1);
    *flag = 1;                              // Assume reconstruction is successful. If not, it will be set to zero below.

    // Perform DLCM Reconstruction entry-wise (cannot be parallel because of lcm)
    for(int i=0; i<n; i++)
    {
        // Compute a: Multiply current entry by LCM and reduce modulo prime
        mpz_mul(a, IntSoln[i], LCM);
        mpz_mod(a, a, prime);

        // If a<b, set x=a/LCM
        if(mpz_cmp(a, b)<0)
        {
            mpq_set_num(soln[i], a);
            mpq_set_den(soln[i], LCM);
        }

        // If a>u, set x=(a-prime)/LCM
        else if(mpz_cmp(a,u)>0)
        {
            mpz_sub(mpq_numref(soln[i]), a, prime);
            mpq_set_den(soln[i], LCM);
        }

        // Otherwise, perform RR iterations
        else
        {
            // Initialize v and w
            mpz_set(v0, prime);
            mpz_set_si(v1, 0);    // Initial v = [prime, 0]
            mpz_set(w0, a);
            mpz_set_si(w1, 1);    // Initial w = [a_i, 1]

            // Compute bound on w1 as d=ceil(b/lcm)
            mpz_cdiv_q(d, b, LCM);

            // Perform RR iterations
            while(mpz_cmp(w0,b)>=0 && mpz_cmpabs(w1, d)<=0)
            {
                // Compute q (quotient) and z0 (remainder) of division v0/w0
                mpz_fdiv_qr(q, z0, v0, w0);

                // Replace v0 and w0
                mpz_set(v0, w0);
                mpz_set(w0, z0);

                // Compute z1 = v1 - q*w1
                mpz_mul(z1, q, w1);
                mpz_sub(z1, v1, z1);

                // Replace v1 and w1
                mpz_set(v1, w1);
                mpz_set(w1, z1);
            }

            // Record the result
            mpq_set_num(soln[i], w0);               // Reconstructed numerator = w0
            mpz_mul(mpq_denref(soln[i]), w1, LCM);  // Reconstructed denominator = LCM*w1

            // Update LCM
            mpz_abs(temp, w1);
            mpz_mul(LCM, LCM, temp);

            // If LCM exceeds HadBound, set flag and terminate reconstruction
            if(mpz_cmp(LCM, HadBound)>0)
            {
                *flag = 0; // Report failure
                break;
            }
        }

        // Simplify entry i
        mpq_canonicalize(soln[i]);
    }

    // Free all vars
    mpz_clear(LCM);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(u);
    mpz_clear(v0);
    mpz_clear(v1);
    mpz_clear(w0);
    mpz_clear(w1);
    mpz_clear(q);
    mpz_clear(z0);
    mpz_clear(z1);
    mpz_clear(HadBound);
    mpz_clear(d);
    mpz_clear(temp);
}



// Description: This function computes the Rational Reconstruction (RR) numerator of a given integer (modulo a given modulus).
void RRNumCalc(
    mpz_t num,      // This will hold the computed numerator value. This will be modified by function.
    mpz_t x,        // The integer to reconstruct. Not modified by function.
    mpz_t p         // The given prime modulus. Not modified by function.
    )
{
    // Declare vars
    mpz_t bound;        // This will hold the RR bound b
    mpz_t u;            // Temporary variable to compute the numerator value
    mpz_t v;            // Temporary variable to compute the numerator value
    mpz_t result;       // This will hold the computed numerator value
    mpz_t modbound;     // This will hold the symmetric representation bound

    // Initialize vars to 0
    mpz_init(bound);
    mpz_init(result);
    mpz_init(u);
    mpz_init(v);
    mpz_init(modbound);

    // Compute the RR bound b
    mpz_set(bound, p);                 // Set b = p^i
    mpz_cdiv_q_2exp(bound, bound, 1);  // Set b = (p^i)/2
    mpz_sqrt(bound, bound);            // This will set b = floor (sqrt{(p^i)/2})
    mpz_add_ui(bound, bound, 1);       // Add 1 to get b = ceil (sqrt{(p^i)/2})

    // Initialize v and result
    mpz_set(v, p);                     // Set v = prime^i
    mpz_abs(result, x);                // Set result = x

    // Compute numerator via mod operations (using unsymmetric respresentation)
    while(mpz_cmp(result, bound)>=0)
    {
        mpz_set(u, v);
        mpz_set(v, result);
        mpz_mod(result, u, v);         // Compute z_k = z_{k-2} MOD z_{k-1}. First iteration computes z_1 = prime^i MOD x
    }

    // Convert final answer to symmetric rep
    mpz_cdiv_q_2exp(modbound, v, 1);
    if(mpz_cmp(result, modbound)>0)
    {
        mpz_sub(result, result, v);
    }

    // Record the result
    mpz_abs(num, result);

    // Free vars
    mpz_clear(bound);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(result);
    mpz_clear(modbound);
}



// Description: This function checks if at least t entries of row i of the given matrix match. If entries match, return 1. Otherwise, return 0.
int CheckNumMatrix(
    mpz_t* numMatrix,   // The matrix to check. Not modified by function.
    int i,              // The row of the matrix to check. Not modified by function.
    int t               // The length of the row to check. Not modified by function.
    )
{
    // Check if next entry is the same as the current
    for(int j = 0; j<t; j++)
    {
        if(mpz_cmp(Matrix_2D(numMatrix, i, j, (t+1)),Matrix_2D(numMatrix, i, j+1, (t+1)))!=0) // Next entry does not match, return failure.
        {
            return 0;
        }
    }

    // All entries match, report success
    return 1;
}


// Description: This function checks if all entries of the given vector of flags are equal to 1. If all entries are equal to 1, return 1. Otherwise, return 0.
int CheckFlags(
    int* vectorFlags,   // The vector of flags to check. Not modified by function.
    int t               // The length of the vector to check. Not modified by function.
    )
{
    // Check if each entry is equal to 1
    for(int i = 0; i<t; i++)
    {
        if(vectorFlags[i]!=1) // Current entry does not equal 1, return failure.
        {
            return 0;
        }
    }

    // All entries equal to 1, report success
    return 1;
}


// Description: This function frees a mpz_t array.
void FreeMPZArray(
    mpz_t* x,  // The array to be freed. This will be modified by function.
    int n      // The length of the array to be freed. Not modified by function.
    )
{
    // Clear the array entry-wise
    for (int i = 0; i < n; i++)
    {
        mpz_clear(x[i]); // Clear entry i
    }

    // Free var
    free(x);
}


// Description: This function frees an mpq_t array
void FreeMPQArray(
    mpq_t* x,  // The array to be freed. This will be modified by the function.
    int n      // The length of the array to be freed. Not modified by function.
    )
{
    // Clear the array entry-wise
    for (int i = 0; i < n; i++)
    {
        mpq_clear(x[i]); // Clear entry i
    }

    // Free var
    free(x);
}

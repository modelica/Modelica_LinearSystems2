within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dgesvx
  "Solve real system of linear equations A**T*X=B with a B matrix with LAPACK routine DGESVX"
  extends Modelica.Icons.Function;
  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),:];
  output Real X[size(A, 1),size(B, 2)];
  output Integer info;
  output Real rcond;

protected
  Real Awork[size(A, 1),size(A, 2)]=A;
  Real Bwork[size(B, 1),size(B, 2)]=B;
  Real AF[size(A, 1),size(A, 2)];
  Real R[size(A, 1)];
  Real C[size(A, 1)];
  Real ferr[size(B, 2)];
  Real berr[size(B, 2)];
  Real work[4*size(A, 1)];
  Integer ipiv[size(A, 1)];
  Integer iwork[size(A, 1)];


external "FORTRAN 77" dgesvx(
    "N",
    "T",
    size(A, 1),
    size(B, 2),
    Awork,
    size(A, 1),
    AF,
    size(A, 1),
    ipiv,
    "N",
    R,
    C,
    B,
    size(A, 1),
    X,
    size(A, 1),
    rcond,
    ferr,
    berr,
    work,
    iwork,
    info)                 annotation (Library="lapack");
  annotation (Documentation(info="Lapack documentation:
    Purpose   
    =======   
    DGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   
    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is 
  
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.   
    Arguments   
    =========   
    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   
    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   
    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   
    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   
    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   
    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   
    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   
    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization 
  
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   
"));
end dgesvx;

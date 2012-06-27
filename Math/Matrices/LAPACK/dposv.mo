within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dposv
  "Compute the solution to A * X = B, where A is a symmetric positive definite matrix"

  extends Modelica.Icons.Function;
  input Real A[:, size(A,1)] "Real symmetric positive definite matrix A";
  input Real B[size(A,2),:] "Right hand side of A*X = B";
  input Boolean upper=true "True if the upper triangle of A is provided";

  output Real X[size(B,1),size(B,2)]=B "Solution of A*X = B";
  output Integer info;
protected
  String uplo=if upper then "U" else "L";
  Integer n=size(A,1);
  Integer nrhs=size(B,2);
  Integer lda=max(1,n);
  Integer ldb=max(1,n);

  external "FORTRAN 77" dposv(uplo, n, nrhs, A, lda, B, ldb, info) annotation (Library="Lapack");

  annotation (
    Documentation(info="Lapack documentation
 

    Purpose   
    =======   

    DPOSV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric positive definite matrix and X and B   
    are N-by-NRHS matrices.   

    The Cholesky decomposition is used to factor A as   
       A = U**T* U,  if UPLO = 'U', or   
       A = L * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is a lower triangular   
    matrix.  The factored form of A is then used to solve the system of   
    equations A * X = B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i of A is not   
                  positive definite, so the factorization could not be   
                  completed, and the solution has not been computed.   

    =====================================================================   
"));

end dposv;

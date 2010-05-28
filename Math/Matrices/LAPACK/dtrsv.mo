within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dtrsv
  "Solve one of the matrix equations  op( A )*x = B where A is upper or lower triangular matrix. BLAS routine"

  input Real A[:,:] "Input matrix A";
  input Real b[size(A,2)] "Input vector b";
  input Boolean upper=true "True if A is upper triangular";
  input Boolean trans=false "True if op(A) means transposed(A)";
  input Boolean unitTriangular=false
    "True if A is unit triangular, i.e. all diagonal elements of A are equal to 1";

  output Real x[size(b,1)]=b "Solution x";
protected
  String uplo=if upper then "U" else "L";
  String transA=if trans then "T" else "N";
  String diag=if unitTriangular then "U" else "N";
  Integer n=size(A, 2) "Number of columns of B";
  Integer lda=max(1,n) "Leading dimension of A";
  Integer incx=1;
external "Fortran 77" dtrsv(uplo, transA, diag, n, A, lda, x, incx) annotation(Library = {"lapack"});

  annotation (Documentation(info="
    
 Purpose   
    =======   
    DTRSV  solves one of the systems of equations   
       A*x = b,   or   A'*x = b,   
    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   
    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   
    Arguments   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   
                TRANS = 'N' or 'n'   A*x = b.   
                TRANS = 'T' or 't'   A'*x = b.   
                TRANS = 'C' or 'c'   A'*x = b.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular matrix and the strictly lower triangular part of   
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular matrix and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of   
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   
    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten   
             with the solution vector x.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.     "));
end dtrsv;

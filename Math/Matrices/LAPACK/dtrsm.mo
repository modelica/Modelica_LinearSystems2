within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dtrsm
  "Solve one of the matrix equations  op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, where A is triangular matrix. BLAS routine"

  input Real A[:,:] "Input matrix A";
  input Real B[:,:] "Input matrix B";
  input Real alpha=1 "Factor alpha";
  input Boolean right=true "True if A is right multiplication";
  input Boolean upper=true "True if A is upper triangular";
  input Boolean trans=false "True if op(A) means transposed(A)";
  input Boolean unitTriangular=false
    "True if A is unit triangular, i.e. all diagonal elements of A are equal to 1";

  output Real X[size(B,1),size(B,2)]=B "Solution X";
protected
  String side=if right then "R" else "L";
  String uplo=if upper then "U" else "L";
  String transA=if trans then "T" else "N";
  String diag=if unitTriangular then "U" else "N";
  Integer m=size(B, 1) "Number of rows of B";
  Integer n=size(B, 2) "Number of columns of B";
  Integer lda=if right then max(1,n) else max(1,m) "First dimension of A";
  Integer ldb=max(1,m) "First dimension of B";

external "Fortran 77" dtrsm(side, uplo, transA, diag, m, n, alpha, A, lda, X, ldb) annotation(Library = {"lapack"});

  annotation (Documentation(info="
    
   Purpose   
    =======   
    DTRSM  solves one of the matrix equations   
       op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,   
    where alpha is a scalar, X and B are m by n matrices, A is a unit, or   
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
       op( A ) = A   or   op( A ) = A'.   
    The matrix X is overwritten on B.   
    Arguments   
    ==========   
    SIDE   - CHARACTER*1.   
             On entry, SIDE specifies whether op( A ) appears on the left   
             or right of X as follows:   
                SIDE = 'L' or 'l'   op( A )*X = alpha*B.   
                SIDE = 'R' or 'r'   X*op( A ) = alpha*B.   
             Unchanged on exit.   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n'   op( A ) = A.   
                TRANSA = 'T' or 't'   op( A ) = A'.   
                TRANSA = 'C' or 'c'   op( A ) = A'.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular   
             as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at   
             least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be   
             at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
             zero then  A is not referenced and  B need not be set before   
             entry.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m   
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
             upper triangular part of the array  A must contain the upper   
             triangular matrix  and the strictly lower triangular part of   
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
             lower triangular part of the array  A must contain the lower   
             triangular matrix  and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
             A  are not referenced either,  but are assumed to be  unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
             LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
             then LDA must be at least max( 1, n ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must   
             contain  the  right-hand  side  matrix  B,  and  on exit  is   
             overwritten by the solution matrix  X.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   LDB  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   "));
end dtrsm;

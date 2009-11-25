within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dgecon
  "Estimates the reciprocal of the condition number of a general real matrix A"

  input Real LU_of_A[:,:] "LU factroization of a real matrix A";
  input Boolean inf=false
    "Is true if infinity norm is used and false for 1-norm";
  input Real anorm "norm of A";
  output Real rcond "Reciprocal condition number of A";
  output Integer info;
protected
  Integer n=size(LU_of_A, 2);
  Integer lda=max(1,size(LU_of_A,1));
  Real work[4*size(LU_of_A,2)];
  Integer iwork[size(LU_of_A,2)];
  String norm = if inf then "I" else "1";

  annotation (Documentation(info=" 
 
  Purpose   
    =======   

    DGECON estimates the reciprocal of the condition number of a general   
    real matrix A, in either the 1-norm or the infinity-norm, using   
    the LU factorization computed by DGETRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as   
       RCOND = 1 / ( norm(A) * norm(inv(A)) ).   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    ANORM   (input) DOUBLE PRECISION   
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.   

    RCOND   (output) DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================  
"));
external "Fortran 77" dgecon(
    norm,
    n,
    LU_of_A,
    lda,
    anorm,
    rcond,
    work,
    iwork,
    info) annotation(Library = {"lapack"});

end dgecon;

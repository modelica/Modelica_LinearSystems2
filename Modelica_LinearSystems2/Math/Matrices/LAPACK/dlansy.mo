within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dlansy "Norm of a symmetric matrix"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Real symmetric matrix A";
  input String norm="1" "Specifies the norm, i.e. 1, I, F, M";
  input Boolean upper = true
    "Specifies whether the upper or lower triangular part of A is referenced";
  output Real anorm "Norm of A";
protected
  String uplo=if upper then "U" else "L";
  Integer n=size(A,1);
  Integer lda=max(1,n);
  Real work[2*n];

external "Fortran 77" dlansy2(norm, uplo, n, A, lda, work, anorm) annotation (Include="
  #include<f2c.h>
   #include <stdio.h>

extern  doublereal dlansy_(char *, char *, integer *, doublereal *, integer *, doublereal *);

int dlansy2_(char *norm, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *work, doublereal *anorm)
{
  *anorm=dlansy_(norm, uplo, n, a, lda, work);
  return 0;
}", Library={"lapack"});
annotation ( Documentation(info="Lapack documentation:

   Purpose
   =======

   DLANSY  returns the value of the one norm,  or the Frobenius norm, or
   the  infinity norm,  or the  element of  largest absolute value  of a
   real symmetric matrix A.

   Description
   ===========

   DLANSY returns the value

      DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
               (
               ( norm1(A),         NORM = '1', 'O' or 'o'
               (
               ( normI(A),         NORM = 'I' or 'i'
               (
               ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

   where  norm1  denotes the  one norm of a matrix (maximum column sum),
   normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
   normF  denotes the  Frobenius norm of a matrix (square root of sum of
   squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.

   Arguments
   =========

   NORM    (input) CHARACTER*1
           Specifies the value to be returned in DLANSY as described
           above.

   UPLO    (input) CHARACTER*1
           Specifies whether the upper or lower triangular part of the
           symmetric matrix A is to be referenced.
           = 'U':  Upper triangular part of A is referenced
           = 'L':  Lower triangular part of A is referenced

   N       (input) INTEGER
           The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
           set to zero.

   A       (input) DOUBLE PRECISION array, dimension (LDA,N)
           The symmetric matrix A.  If UPLO = 'U', the leading n by n
           upper triangular part of A contains the upper triangular part
           of the matrix A, and the strictly lower triangular part of A
           is not referenced.  If UPLO = 'L', the leading n by n lower
           triangular part of A contains the lower triangular part of
           the matrix A, and the strictly upper triangular part of A is
           not referenced.

   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(N,1).

   WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
           where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
           WORK is not referenced.

   =====================================================================
"));

end dlansy;

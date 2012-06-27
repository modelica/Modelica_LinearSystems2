within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dtrtri "Computes the inverse of a triangular real matrix A"

  extends Modelica.Icons.Function;
  input Real A[:, size(A,1)] "Real symmetric positive definite matrix A";
  input Boolean upper=true "True if the upper triangle of A is provided";

  output Real invA[size(A,1),size(A,1)]=A "Inverse of A";
  output Integer info;
protected
  String uplo=if upper then "U" else "L";
  String diag="N";
  Integer n=size(A,1);
  Integer lda=max(1,n);
  external "FORTRAN 77" dtrtri(uplo, diag, n, invA, lda, info) annotation (Library="Lapack");

  annotation (
    Documentation(info="Lapack documentation
   Purpose  
   =======  

   DTRTRI computes the inverse of a real upper or lower triangular  
   matrix A.  

   This is the Level 3 BLAS version of the algorithm.  

   Arguments  
   =========  

   UPLO    (input) CHARACTER*1  
           = 'U':  A is upper triangular;  
           = 'L':  A is lower triangular.  

   DIAG    (input) CHARACTER*1  
           = 'N':  A is non-unit triangular;  
           = 'U':  A is unit triangular.  

   N       (input) INTEGER  
           The order of the matrix A.  N >= 0.  

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)  
           On entry, the triangular matrix A.  If UPLO = 'U', the  
           leading N-by-N upper triangular part of the array A contains  
           the upper triangular matrix, and the strictly lower  
           triangular part of A is not referenced.  If UPLO = 'L', the  
           leading N-by-N lower triangular part of the array A contains  
           the lower triangular matrix, and the strictly upper  
           triangular part of A is not referenced.  If DIAG = 'U', the  
           diagonal elements of A are also not referenced and are  
           assumed to be 1.  
           On exit, the (triangular) inverse of the original matrix, in  
           the same storage format.  

   LDA     (input) INTEGER  
           The leading dimension of the array A.  LDA >= max(1,N).  

   INFO    (output) INTEGER  
           = 0: successful exit  
           < 0: if INFO = -i, the i-th argument had an illegal value  
           > 0: if INFO = i, A(i,i) is exactly zero.  The triangular  
                matrix is singular and its inverse can not be computed.  

   =====================================================================  "));

end dtrtri;

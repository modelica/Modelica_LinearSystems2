within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dlange "Norm of a matrix"

  input Real A[:,:] "Real matrix A";
  input String norm="1" "specifies the norm, i.e. 1, I, F, M";
  output Real anorm "norm of A";
protected
  Integer m=size(A, 1);
  Integer n=size(A,2);
  Integer lda=max(1,m);
  Real work[2*m];

external "Fortran 77" dlange2(norm, m, n, A, lda, work, anorm)
  annotation (Include="
  #include<f2c.h>
   #include <stdio.h> 

extern  doublereal dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *);

int dlange2_(char *norm, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work, doublereal *anorm) 
{
  *anorm=dlange_(norm, m, n, a, lda, work);
  return 0;
}", Library={"lapack"});
annotation ( Documentation(info="<html>


    Purpose   
    =======   

    DLANGE  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real matrix A.   

    Description   
    ===========   

    DLANGE returns the value   

       DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANGE as described   
            above.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.  When M = 0,   
            DLANGE is set to zero.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.  When N = 0,   
            DLANGE is set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(M,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),   
            where LWORK >= M when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   =====================================================================   


</html>"));

end dlange;

within Modelica_LinearSystems2.Math.Matrices.Internal;
function solveSymRight_C
  "Solve real system of linear equations X*A=B where A is symmetrix positive definite"

  extends Modelica.Icons.Function;
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real A[:,size(A, 1)] "Matrix A of X*A = B";
  input Real B[:,:] "Matrix B of X*op(A) = B";
  input Boolean isTriangular=false "True if the A is already lower triangular";
  input Boolean upper=false "True if A is upper triAngular";
  output Real X[size(B, 1),size(B, 2)]=B "Matrix X such that X*A = B";
  output Integer info;
//  Real H[size(A, 1),size(A, 2)]=A;
//  Real AA[size(A,1),size(A,1)]=if upper then A else transpose(A);

//  Real X2[size(B, 1),size(B, 2)];
protected
  String trian=if isTriangular then "T" else "N";
  String uplo=if upper then "U" else "L";
  Integer m=size(B, 1) "Number of rows of B";
  Integer n=size(B, 2) "Number of columns of B";
  Integer lda=max(1,n) "First dimension of A";
  Integer ldb=max(1,m) "First dimension of B";

 external "FORTRAN 77" c_solve2rSym(A, X, trian, uplo, m, n, lda, ldb, info)
  annotation (Include="
  #include<f2c.h>
extern  int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
extern logical lsame_(char *, char *);
//extern int dpotrf_(char *, integer *, doublereal *, integer *, integer *);
extern int dpotrf_(const char* , int  *, double  *, int  *, int  *);

 #include <stdio.h>

int c_solve2rSym_(doublereal *a, doublereal *b, char *trian, char *uplo, integer *m, integer *n, integer *lda, integer *ldb, integer *info)
{
   static logical upp;
   static logical tri;
   static doublereal alpha = 1.0;


   doublereal *aa;


   integer nn=*n;
   integer mm=*m;
   integer llda=*lda;
   integer lldb=*ldb;

   integer i,j;


//    FILE *fileptr;
//    fileptr = fopen(\"test.txt\",\"w\");
//    fprintf(fileptr,\"a = %f, %f, %f, %f, %f, %f\\n\",a[0],a[1],a[2],a[3],a[4],a[5]);

   aa = (doublereal *) malloc((nn*nn+1)*sizeof(doublereal));

   upp = lsame_(uplo, \"U\");
   tri = lsame_(trian, \"T\");

   if(upp)
   {
     for(i=0;i<nn*nn;i++)
       aa[i]=a[i];
   }
   else
   {
     for(i=0;i<nn;i++)
       for(j=0;j<nn;j++)
         aa[j*nn+i]=a[i*nn+j];
   }

// fprintf(fileptr,\"aa = %f, %f, %f\\n %f, %f, %f\\n %f, %f, %f\\n\",aa[0],aa[3],aa[6],aa[1],aa[4],aa[7],aa[2],aa[5],aa[8]);

   if(! tri)
     dpotrf_(\"U\",n,aa,lda,info);


   for(i=1;i<nn;i++)
     for(j=i;j<nn;j++)
       aa[(i-1)*nn+j]=aa[j*nn+i-1];

   dtrsm_(\"R\", \"U\", \"N\", \"N\", m, n, &alpha, aa, lda, b, ldb);
   dtrsm_(\"R\", \"L\", \"N\", \"N\", m, n, &alpha, aa, lda, b, ldb);
// fprintf(fileptr,\"b2 = %f, %f, %f\\n %f, %f, %f\\n %f, %f, %f\\n\",b[0],b[3],b[6],b[1],b[4],b[7],b[2],b[5],b[8]);



   free(aa);
//   fclose(fileptr);
  return 0;
}", Library={"lapack"});
  annotation(Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>", info="<html>
This function solves the equation
<blockquote><pre>

  <strong>X</strong>*<strong>A</strong> = <strong>B</strong>

</pre></blockquote>
where matrix <strong>A</strong> with symmetric positiv definite matrix. The calculation is rather efficient since symmetrie and decomposition of positive definite matrices is exploited.
<p>
Due to symmetrie, Matrix <strong>A</strong> is uniquely defined by a triangle, i.e. the upper or the lower triangular matrix. It is assumed, that the input to describe <strong>A</strong> is either a Cholesky factor or
part of matrix <strong>A</strong> itself. This is defined by the user with the boolean inputs <em>isCholesky</em> and <em>upper</em> which is true when <strong>A</strong> is already Cholesky factor and when <strong>A</strong> is upper triangular respectively.
<p>
Considering the Cholesky decomposition
<blockquote><pre>
         T
  <strong>A</strong> = <strong>L</strong>*<strong>L</strong>

</pre></blockquote>
with lower triangular matrix <strong>L</strong> the equation above could be rewritten as
<blockquote><pre>
       T
  <strong>X</strong>*<strong>L</strong>*<strong>L</strong> = <strong>B</strong>

</pre></blockquote>
which is solved with BLAS function <em>dtrmm</em> applied to a upper triangular matrix and subsequently to a lower triangular matrix.
<p>
In contrast to function <em>solveSymRight</em> this function is implemented in C-code
</p>
</html>"));

end solveSymRight_C;

within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zgetri
  import Modelica_LinearSystems2.Math.Complex;

  input Real LU_real[ :, :]
    "LU factorization of zgetrf of a square matrix, real part";
  input Real LU_imag[size(LU_real,1), size(LU_real,2)]
    "LU factorization of zgetrf of a square matrix, imaginary part";
  input Integer pivots[size(LU_real, 1)] "Pivot vector of zgetrf";
  output Real inv_real[size(LU_real, 1), size(LU_real, 2)]=LU_real
    "Inverse of matrix P*L*U, real part";
  output Real inv_imag[size(LU_real, 1), size(LU_real, 2)]=LU_imag
    "Inverse of matrix P*L*U, imaginary part";
  output Integer info;

  external "FORTRAN 77" c_inter_zgetri(size(LU_real, 1), inv_real, inv_imag, size(LU_real,1), pivots, info);
  annotation (Include="
  #include<f2c.h>
extern  int zgetri_(integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);

int c_inter_zgetri_(integer *n, doublereal *a_real, doublereal *a_imag, integer *lda, integer *pivots,  integer *info)
{
   doublecomplex *a;
   doublecomplex *work;

  integer nn=*n;
  integer mm=*lda;
  integer lwork=-1;
  integer i;

   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));

   for(i=0;i<nn*mm;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }


   work = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex));
   zgetri_(n, a, lda, pivots, work, &lwork, info);
   lwork=(int)(work[0].r);
   free(work);
   work = (doublecomplex *) malloc((lwork+1)*sizeof(doublecomplex));
   zgetri_(n, a, lda, pivots, work, &lwork, info);
   for(i=0;i<nn*mm;i++)
   {
     a_real[i]=a[i].r;
     a_imag[i]=a[i].i;
   }

   free(a);
   free(work);
  return 0;
}", Library={"zlapack"});

end wrapper_zgetri;

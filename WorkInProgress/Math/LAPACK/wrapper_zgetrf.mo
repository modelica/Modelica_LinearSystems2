within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zgetrf
  import Modelica_LinearSystems2.Math.Complex;

//  input Modelica_LinearSystems2.Math.Complex A[ :, :] "Square or rectangular matrix";
  input Real A_real[ :, :] "Square or rectangular matrix";
  input Real A_imag[size(A_real,1), size(A_real,2)]
    "Square or rectangular matrix";
  output Real LU_real[size(A_real, 1), size(A_real, 2)]=A_real
    "LU factorization in packed format, real part";
  output Real LU_imag[size(A_real, 1), size(A_real, 2)]=A_imag
    "LU factorization in packed format, imaginary part";
  output Integer pivots[min(size(A_real, 1), size(A_real, 2))] "Pivot vector";
//  output Real inv_real[size(LU_real, 1), size(LU_real, 2)]=LU_real "Real part of inverse of matrix P*L*U";
//  output Real inv_imag[size(LU_imag, 1), size(LU_imag, 2)]=LU_imag "Real part of inverse of matrix P*L*U";

  output Integer info;

  external "FORTRAN 77" c_inter_zgetrf(size(A_real, 1), size(A_real, 2), LU_real, LU_imag, size(A_real,1), pivots, info);
  annotation (Include="
  #include<f2c.h>
extern  int zgetrf_(integer *, integer *, doublecomplex *, integer *, integer *, integer *);

int c_inter_zgetrf_(integer *m, integer *n, doublereal *a_real, doublereal *a_imag, integer *lda, integer *pivots, integer *info) 
{
   doublecomplex *a;
     
  integer nn=*n;
  integer mm=*m;
  integer i;
  
   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));      
        
   for(i=0;i<nn*mm;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }
   
   zgetrf_(m, n, a, lda, pivots, info);
  
   for(i=0;i<nn*mm;i++)
   {
     a_real[i]=a[i].r;
     a_imag[i]=a[i].i;
   }
     
   free(a);
  return 0;
}", Library={"zlapack"});

end wrapper_zgetrf;

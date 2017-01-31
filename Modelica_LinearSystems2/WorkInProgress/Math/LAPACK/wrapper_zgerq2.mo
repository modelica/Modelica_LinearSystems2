within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zgerq2 "Wrapper for lapack function zgerq2"
  import Modelica_LinearSystems2.Math.Complex;

//  input Modelica_LinearSystems2.Math.Complex A[ :, :] "Square or rectangular matrix";
  input Real A_real[ :, :] "Square or rectangular matrix";
  input Real A_imag[size(A_real,1), size(A_real,2)]
    "Square or rectangular matrix";
  output Real RQ_real[size(A_real, 1), size(A_real, 2)]=A_real
    "RQ factorization in packed format, real part";
  output Real RQ_imag[size(A_real, 1), size(A_real, 2)]=A_imag
    "RQ factorization in packed format, imaginary part";
  output Real tau_real[min(size(A_real, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  output Real tau_imag[min(size(A_real, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  output Integer info;

protected
  Integer n=size(A_real, 1) "Row dimension of A";
  Integer m=size(A_real, 2) "Column dimension of A";
  external "FORTRAN 77" c_inter_zgerq2(n, m, RQ_real, RQ_imag, n, tau_real, tau_imag, info);
  annotation (Include="
  #include<f2c.h>
extern  int zgerq2_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);

int c_inter_zgerq2_(integer *n, integer *m, doublereal *a_real, doublereal *a_imag, integer *lda, doublereal *tau_real, doublereal *tau_imag,  integer *info)
{
   doublecomplex *a;
   doublecomplex *tau;
   doublecomplex *work;

  integer nn=*n;
  integer mm=*m;
  integer mn=min(mm,nn);
  integer i;

   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));
   tau = (doublecomplex *) malloc((mn+1)*sizeof(doublecomplex));
   work = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex));

   for(i=0;i<nn*mm;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }

   zgerq2_(n, m, a, lda, tau, work, info);

   for(i=0;i<nn*mm;i++)
   {
     a_real[i]=a[i].r;
     a_imag[i]=a[i].i;
   }
   for(i=0;i<mn;i++)
   {
     tau_real[i]=tau[i].r;
     tau_imag[i]=tau[i].i;
   }

   free(a);
   free(tau);
   free(work);
  return 0;
}", Library={"zlapack"});

end wrapper_zgerq2;

within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zgeqrf
  import Modelica_LinearSystems2.Math.Complex;

//  input Modelica_LinearSystems2.Math.Complex A[ :, :] "Square or rectangular matrix";
  input Real A_real[ :, :] "Square or rectangular matrix";
  input Real A_imag[size(A_real,1), size(A_real,2)]
    "Square or rectangular matrix";
  output Real QR_real[size(A_real, 1), size(A_real, 2)]=A_real
    "QR factorization in packed format, real part";
  output Real QR_imag[size(A_real, 1), size(A_real, 2)]=A_imag
    "QR factorization in packed format, imaginary part";
  output Real tau_real[min(size(A_real, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  output Real tau_imag[min(size(A_real, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  output Integer info;

protected
  Integer ncol=size(A_real, 2) "Column dimension of A";
  Integer lwork[2*ncol];
  Real work_real[ncol] "work array real";
  Real work_imag[ncol] "work array imaginary";
  external "FORTRAN 77" c_inter_zgeqrf(size(A_real, 1), ncol, QR_real, QR_imag, size(A_real,1), tau_real, tau_imag, info);
  annotation (Include="
  #include<f2c.h>
extern  int zgeqrf_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, integer *);

int c_inter_zgeqrf_(integer *m, integer *n, doublereal *a_real, doublereal *a_imag, integer *lda, doublereal *tau_real, doublereal *tau_imag,  integer *info)
{
   doublecomplex *a;
   doublecomplex *tau;
   doublecomplex *work;

  integer nn=*n;
  integer mm=*m;
  integer lwork=-1;
  integer mn=min(mm,nn);
  integer i;

   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));
   tau = (doublecomplex *) malloc((mn+1)*sizeof(doublecomplex));

   for(i=0;i<nn*mm;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }

   work = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex));
   zgeqrf_(m, n, a, lda, tau, work, &lwork, info);
   lwork=(int)(work[0].r);
   free(work);
   work = (doublecomplex *) malloc((lwork+1)*sizeof(doublecomplex));
   zgeqrf_(m, n, a, lda, tau, work, &lwork, info);

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

end wrapper_zgeqrf;

within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zungrq
  import Modelica_LinearSystems2.Math.Complex;

  input Real RQ_real[ :, :] "Square or rectangular matrix";
  input Real RQ_imag[size(RQ_real,1), size(RQ_real,2)]
    "Square or rectangular matrix";
  input Real tau_real[min(size(RQ_real, 1), size(RQ_real, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  input Real tau_imag[min(size(RQ_real, 1), size(RQ_real, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  output Real Q_real[size(RQ_real, 1),size(RQ_real, 2)]=RQ_real;
  output Real Q_imag[size(RQ_real, 1),size(RQ_real, 2)]=RQ_imag;
  output Integer info;
//  Integer info;
  external "FORTRAN 77" c_inter_zungrq(size(RQ_real, 1), size(RQ_real, 2), size(tau_real,1),Q_real, Q_imag, size(RQ_real,1), tau_real, tau_imag, info);
  annotation (Include="
  #include<f2c.h>
extern  int zungrq_(integer *, integer *, integer *, doublecomplex *, integer *,doublecomplex *, doublecomplex *, integer *, integer *);

int c_inter_zungrq_(integer *m, integer *n, integer *k, doublereal *rq_real, doublereal *rq_imag, integer *lda, doublereal *tau_real, doublereal *tau_imag, integer *info) 
{
   doublecomplex *a;
   doublecomplex *tau;
   doublecomplex *work;
     
   integer nn=*n;
   integer mm=*m;
   integer mn=min(mm,nn);
   integer i;
   integer lwork=-1;
  
   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));      
   tau = (doublecomplex *) malloc((mn+1)*sizeof(doublecomplex));      

        
   for(i=0;i<nn*mm;i++)
   {
     a[i].r=rq_real[i];
     a[i].i=rq_imag[i];
   }
    for(i=0;i<mn;i++)
    {
      tau[i].r=tau_real[i];
      tau[i].i=tau_imag[i];
    }
  
   work = (doublecomplex *) malloc(nn*sizeof(doublecomplex));
   zungrq_(m, n, k, a, lda, tau, work, &lwork, info);
   lwork=(int)(work[0].r);
   free(work);
   work = (doublecomplex *) malloc((lwork+1)*sizeof(doublecomplex));
   zungrq_(m, n, k, a, lda, tau, work, &lwork, info);

   for(i=0;i<nn*mm;i++)
   {
     rq_real[i]=a[i].r;
     rq_imag[i]=a[i].i;
   }
  
   free(a);
   free(tau);
   free(work);
  return 0;
}", Library={"zlapack"});

end wrapper_zungrq;

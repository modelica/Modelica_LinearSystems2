within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function wrapper_zunmrq "Wrapper for lapack function zunmrq"

  input Real RQ_real[ :, :] "Square or rectangular matrix";
  input Real RQ_imag[size(RQ_real,1), size(RQ_real,2)]
    "Square or rectangular matrix";
  input Real tau_real[min(size(RQ_real, 1), size(RQ_real, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  input Real tau_imag[min(size(RQ_real, 1), size(RQ_real, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  input Real C_real[:,:];
  input Real C_imag[size(C_real,1), size(C_real,2)];
  input Boolean left=true;
  input Boolean trans=false;

  output Real QC_real[size(C_real, 1),size(C_real, 2)]=C_real;
  output Real QC_imag[size(C_real, 1),size(C_real, 2)]=C_imag;
  output Integer info;

protected
  Integer lwork = if left then 2*size(C_real,2) else 2*size(C_real,1);
  Integer dim2=if left then size(C_real,1) else size(C_real,2);
  String side = if left then "L" else "R";
  String herm = if trans then "C" else "N";

  external "FORTRAN 77" c_inter_zunmrq(side, herm,  size(QC_real,1),size(QC_real,2), size(tau_real,1), dim2, RQ_real, RQ_imag,  size(tau_real,1), tau_real, tau_imag, QC_real, QC_imag, size(C_real,1), lwork, info);
  annotation (Include="
  #include<f2c.h>
extern  int zunmrq_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);

 #include <stdio.h>

int c_inter_zunmrq_(char *side, char *trans, integer *m, integer *n, integer *k, integer *dim2, doublereal *rq_real, doublereal *rq_imag, integer *lda, doublereal *tau_real, doublereal *tau_imag,
        doublereal *c_real, doublereal *c_imag, integer *ldc, integer *lwork, integer *info)
{
   doublecomplex *rq;
   doublecomplex *c;
   doublecomplex *tau;
   doublecomplex *work;

   integer nn=*n;
   integer mm=*m;
   integer kk=*k;
   integer ddim2=*dim2;
   integer llwork=*lwork;
   integer i;
     FILE *fileptr;

   fileptr = fopen(\"test.txt\",\"w\");

   rq = (doublecomplex *) malloc((kk*ddim2+1)*sizeof(doublecomplex));
   c = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));
   tau = (doublecomplex *) malloc((kk+1)*sizeof(doublecomplex));
   work = (doublecomplex *) malloc((llwork+1)*sizeof(doublecomplex));

   for(i=0;i<kk*ddim2;i++)
   {
     rq[i].r=rq_real[i];
     rq[i].i=rq_imag[i];
   }
   for(i=0;i<nn*mm;i++)
   {
     c[i].r=c_real[i];
     c[i].i=c_imag[i];
   }
    for(i=0;i<kk;i++)
    {
      tau[i].r=tau_real[i];
      tau[i].i=tau_imag[i];
    }

 zunmrq_(side, trans, m, n, k, rq, lda, tau, c, ldc, work, lwork, info);
 fprintf(fileptr,\"info = %d\\n\",*info);

   for(i=0;i<nn*mm;i++)
   {
     c_real[i]=c[i].r;
     c_imag[i]=c[i].i;
   }

   free(rq);
   free(c);
   free(tau);
   free(work);
   fclose(fileptr);
  return 0;
}", Library={"zlapack"});

end wrapper_zunmrq;

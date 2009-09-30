within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function wrapper_zgesvd

  input Real A_real[:,:] "Square or rectangular matrix";
  input Real A_imag[size(A_real, 1),size(A_real, 2)]
    "Square or rectangular matrix";
  output Real sigma[min(size(A_real, 1), size(A_real, 2))] "singular values";
  output Real U_real[size(A_real, 1),size(A_real, 1)];//=identity(size(A_real, 1)) "Left orthogonal matrix, real part";
  output Real U_imag[size(A_real, 1),size(A_real, 1)];//=zeros(size(A_real, 1), size(A_real, 1)) "Left orthogonal matrix, imaginary part";
  output Real VT_real[size(A_real, 2),size(A_real, 2)];//=identity(size(A_real, 2)) "Left orthogonal matrix, real part";
  output Real VT_imag[size(A_real, 2),size(A_real, 2)];//=zeros(size(A_real, 2),  size(A_real, 2)) "Left orthogonal matrix, imaginary part";
  output Integer info;

protected
  Integer nrow=size(A_real, 1) "Column dimension of A";
  Integer ncol=size(A_real, 2) "Column dimension of A";
  Real rwork[5*min(nrow, ncol)];
external"FORTRAN 77" c_inter_zgesvd(
    "A",
    "A",
    size(A_real, 1),
    size(A_real, 2),
    A_real,
    A_imag,
    size(A_real, 1),
    sigma,
    U_real,
    U_imag,
    size(A_real, 1),
    VT_real,
    VT_imag,
    size(A_real, 2),
    rwork,
    info);
  annotation (Include="
  #include<f2c.h>
extern  int zgesvd_(char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *, 
        integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, integer *);

int c_inter_zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
        doublereal *a_real ,doublereal *a_imag, integer *lda, doublereal *s, doublereal *u_real, doublereal *u_imag, 
        integer *ldu, doublereal *vt_real, doublereal *vt_imag, integer *ldvt, doublereal *rwork, integer *info) 
{
   doublecomplex *a;
   doublecomplex *u;
   doublecomplex *vt;
   doublecomplex *work;
     
  integer nn=*n;
  integer mm=*m;
  integer lwork=-1;
  integer minmn=min(mm,nn);
  integer maxmn=max(mm,nn);
  integer i;
  
   a = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));      
   u = (doublecomplex *) malloc((mm*mm+1)*sizeof(doublecomplex));      
   vt = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
        
   for(i=0;i<nn*mm;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }
   
   work = (doublecomplex *) malloc((2*minmn+maxmn+1)*sizeof(doublecomplex));
   zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &lwork, rwork, info);
   lwork=(int)(work[0].r);
   free(work);
   work = (doublecomplex *) malloc((lwork+1)*sizeof(doublecomplex));
   zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &lwork, rwork, info);
  
 
 
   for(i=0;i<mm*mm;i++)
   {
     u_real[i]=u[i].r;
     u_imag[i]=u[i].i;
   }
   for(i=0;i<nn*nn;i++)
   {
     vt_real[i]=vt[i].r;
     vt_imag[i]=vt[i].i;
   }
   
   free(a);
   free(u);
   free(vt);
   free(work);
  return 0;
}", Library={"zlapack"});

end wrapper_zgesvd;

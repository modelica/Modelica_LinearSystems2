within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function wrapper_zgesv

  input Real A_real[:,size(A_real,1)] "Square or rectangular matrix";
  input Real A_imag[size(A_real, 1),size(A_real, 2)]
    "Square or rectangular matrix";
  input Real B_real[size(A_real,1),:] "Square or rectangular matrix";
  input Real B_imag[size(B_real, 1),size(B_real, 2)]
    "Square or rectangular matrix";
  output Real X_real[size(A_real, 1), size(B_real, 2)]=B_real;
  output Real X_imag[size(A_real, 1), size(B_real, 2)]=B_imag;
  output Integer info;

protected
  Real Awork_real[size(A_real, 1), size(A_real, 1)]=A_real;
  Real Awork_imag[size(A_imag, 1), size(A_imag, 1)]=A_imag;
  Integer ipiv[size(A_real, 1)];

external"FORTRAN 77" c_inter_zgesv(size(A_real, 1), size(B_real, 2), Awork_real, Awork_imag, size(A_real, 1), ipiv,  X_real,  X_imag, size(A_real, 1), info);
  annotation (Include="
  #include<f2c.h>
extern  int zgesv_(integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);

int c_inter_zgesv_(integer *n, integer *nrhs, doublereal *a_real, doublereal *a_imag, integer *lda, integer *ipiv, doublereal *b_real, doublereal *b_imag, integer *ldb, integer *info) 
{
   doublecomplex *a;
   doublecomplex *b;
     
   integer nn=*n;
   integer mm=*nrhs;
   integer i;
  
   a = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
   b = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));      
        
   for(i=0;i<nn*nn;i++)
   {
     a[i].r=a_real[i];
     a[i].i=a_imag[i];
   }
   for(i=0;i<nn*mm;i++)
   {
     b[i].r=b_real[i];
     b[i].i=b_imag[i];
   }
   
   zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
  
   for(i=0;i<nn*mm;i++)
   {
     b_real[i]=b[i].r;
     b_imag[i]=b[i].i;
   }
   
   free(a);
   free(b);
  return 0;
}", Library={"zlapack"});

end wrapper_zgesv;

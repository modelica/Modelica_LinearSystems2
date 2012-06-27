within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function wrapper_calcK
  "Computes the feedback matrix from the assigned eigenvalues, closed loop eigenvectors and the B matrix factorization"

  input Real A[:,size(A, 1)] "Real square system matrix";
  input Real U0[size(A, 1),:] "U0 and Z are the decompositions of B";
  input Real Z[size(U0, 2),size(U0, 2)] "Z and U0 are the decompositions of B";
  input Real gamma_real[size(A,1)] "Eigenvalue vector, real part";
  input Real gamma_imag[size(A,1)] "Eigenvalue vector, imaginary part";
  input Real X_real[n,n] "Eigenvectors, real part";
  input Real X_imag[n,n] "Eigenvectors, imaginary part";
  input Integer nre "Number of real eigenvalues";

  output Real K[size(U0,2), size(A,1)] "Feedback matrix";

protected
  Integer n=size(A,1);
  Integer m=size(U0,2);

external"FORTRAN 77" c_inter_calcK(n, m, A, U0, Z, gamma_real, gamma_imag, X_real, X_imag, nre, K);

  annotation (Include="
  #include<f2c.h>
  
  
extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
extern  int zgesv_(integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);

int c_inter_calcK_(integer *n, integer *m, doublereal *a, doublereal *u0, doublereal *z, doublereal *gamma_real, doublereal *gamma_imag, doublereal *x_real, doublereal *x_imag, integer *nre, doublereal *k) 
{
   static integer c1 = 1;
   
   doublecomplex *gamma;
   doublecomplex *xh;
   doublecomplex *x2;
   doublereal *ma;
   doublereal *mp;
   
   doublecomplex nc={0.0,0.0};
   doublecomplex ic={1.0,0.0};
   doublereal ndd=0.0;
   doublereal idd=1.0;
     
    
   integer nn=*n;
   integer mm=*m;
   integer nnre=*nre;// number of real eigenvalues
   integer nncp=(nn-nnre)/2;// number of conjugated complex eigenvalues
   integer i;
   integer ii;
   integer info;
   integer *ipiv;
   
   char *all=\"A\";
   char *conj=\"C\";
   char *trans=\"T\";
   char *no=\"N\";
   char *fro=\"F\";
   
  
   ma = (doublereal *) malloc((nn*nn+1)*sizeof(doublereal));  
   mp = (doublereal *) malloc((nn*mm+1)*sizeof(doublereal));         
   if(nn>mm)  
   {
      gamma = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex)); 

   for(i=0;i<nn;i++)
   {
     gamma[i].r=gamma_real[i];
     gamma[i].i=gamma_imag[i];
   }


   xh = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
   x2 = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
     
   ipiv = (integer *) malloc((nn+1)*sizeof(integer));      
        
   
     
   for(i=0;i<nn;i++)    
     for(ii=0;ii<nn;ii++)    
     {
       x2[ii*nn+i].r = gamma[i].r*x_real[i*nn+ii]-gamma[i].i*x_imag[i*nn+ii];
       x2[ii*nn+i].i = -gamma[i].r*x_imag[i*nn+ii]-gamma[i].i*x_real[i*nn+ii];
       xh[i*nn+ii].r = x_real[ii*nn+i];
       xh[i*nn+ii].i = -x_imag[ii*nn+i];
     }
   
   
   zgesv_(n, n, xh, n, ipiv, x2, n, &info); 
   for(i=0;i<nn;i++)    
     for(ii=0;ii<nn;ii++)
       ma[i*nn+ii] = -x2[ii*nn+i].r+a[i*nn+ii];
   
   dgemm_(no, trans, m, n, m, &idd, z, m, u0, n, &ndd, mp, m);//mp=z*u0'
   dgemm_(no, no, m, n, n, &idd, mp, m, ma, n, &ndd, k, m);//k=mp*ma'
   free(xh);
   free(x2);
   free(ipiv); 
   free(gamma);
   }// if nn>mm
   else
   {
     for(i=0;i<nn*nn;i++)
     {
       ma[i]=a[i];
     }
     for(i=0;i<nn;i++)
     {
       ma[nn*i+i]=ma[nn*i+i]-gamma_real[i];
     }
     for(i=0;i<nncp;i++)
     {
//         ma[nn*(nnre+2*i)+nnre+2*i] = ma[nn*(nnre+2*i)+nnre+2*i]-gamma_real[nnre + 2*i];
//         ma[nn*(nnre+2*i+1)+nnre+2*i+1] = ma[nn*(nnre+2*i+1)+nnre+2*i+1]-gamma_real[nnre + 2*i];
         ma[nn*(nnre+2*i+1)+nnre+2*i] = ma[nn*(nnre+2*i+1)+nnre+2*i]+gamma_imag[nnre + 2*i];
         ma[nn*(nnre+2*i)+nnre+2*i+1] = ma[nn*(nnre+2*i)+nnre+2*i+1]-gamma_imag[nnre + 2*i];
     }
     
      dgemm_(no, trans, m, n, m, &idd, z, m, u0, n, &ndd, mp, m);//mp=z*u0'; z'=inv(z) with z resulting from svd
      dgemm_(no, no, m, n, n, &idd, mp, m, ma, n, &ndd, k, m);//k=mp*ma
   }
   free(mp);
   free(ma);

  return 0;
}", Library={"zlapack"});

end wrapper_calcK;

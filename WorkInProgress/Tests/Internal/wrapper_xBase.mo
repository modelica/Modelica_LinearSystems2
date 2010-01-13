within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function wrapper_xBase
  "Compute the eigenvector bases according to Kautsky algorithm"

  input Real A[:,size(A, 1)] "Real square system matrix";
  input Real B[size(A,1),:] "Real input matrix";
  input Real gamma_real[size(A,1)] "Eigenvalue vector, real part";
  input Real gamma_imag[size(A,1)] "Eigenvalue vector, imaginary part";
  input Integer ncp "Number of complex pairs of eigenvalues";

  output Real U0[size(A,1),size(B,2)] "U0";
  output Real Z[size(B,2),size(B,2)] "Z";
  output Real S_real[size(A,1),size(B,2)*(size(A,1)-ncp)]
    "Eigenvector bases, real part";
  output Real S_imag[size(A,1),size(B,2)*(size(A,1)-ncp)]
    "Eigenvector bases, imaginary part";
  output Integer rankB "Rank of matrix B";

protected
  Integer n=size(A,1);
  Integer m=size(B,2);

external"FORTRAN 77" c_inter_xBase(A, B, n, m, gamma_real, gamma_imag, ncp, U0, Z, S_real, S_imag, rankB);
  annotation (Include="
  #include<f2c.h>
  
extern int zgesvd_(char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, integer *);
//extern int nullSpace_(doublecomplex *, integer *, integer *, doublecomplex *, integer *);

int nullSpace_(doublecomplex *a, integer *n, integer *m, doublecomplex *v, integer *nullity)
{
  integer nn=*n;
  integer mm=*m;
  integer minnm=min(nn,mm);
  integer lworkzsvd=-1;
  integer info;
  integer rank;
  integer i;
  integer ii;
 
    
  doublecomplex *workzsvd;
  doublecomplex *u;
  doublecomplex *vt;
  doublecomplex *aa;
  doublereal *sigma;
  doublereal sigmamax;
  
  doublereal eps;
  doublereal *rwork;
   
  char *all=\"A\";
   

  sigma=(doublereal *) malloc((minnm+1)*sizeof(doublereal)); 
  aa=(doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex)); 
  u=(doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex)); 
  vt=(doublecomplex *) malloc((mm*mm+1)*sizeof(doublecomplex)); 
  rwork = (doublereal *) malloc((5*minnm+1)*sizeof(doublereal)); 
  
  for(i=0;i<nn*mm;i++)
    aa[i] = a[i];

  workzsvd = (doublecomplex *) malloc(max(1,2*minnm+max(mm,nn))*sizeof(doublecomplex));
  zgesvd_(all, all, n, m, aa, n, sigma, u, n, vt, m, workzsvd, &lworkzsvd, rwork, &info);
  lworkzsvd=(int)(workzsvd[0].r);
  free(workzsvd);
  workzsvd = (doublecomplex *) malloc((lworkzsvd+1)*sizeof(doublecomplex));
  zgesvd_(all, all, n, m, aa, n, sigma, u, n, vt, m, workzsvd, &lworkzsvd, rwork, &info);
  
  sigmamax=0.0;
  for(i=0;i<minnm;i++)
    if(sigma[i]>sigmamax)
      sigmamax=sigma[i];
      
  eps=max(nn,mm)*sigmamax*1e-12;
  rank=0;
  i=nn;
  while(i>0)
  {
     if(sigma[i-1]>eps)
     {
        rank=i;
        i=0;
     }
      i=i-1;
  }
  *nullity=mm-rank;
  
  for(i=0;i<mm;i++)
    for(ii=0;ii<mm;ii++)
    {
      v[i*mm+ii].r = vt[ii*mm+i].r;
      v[i*mm+ii].i = -vt[ii*mm+i].i;
    }

  free(workzsvd);
  free(u);
  free(vt);
  free(sigma);
  free(rwork);
  free(aa);
return 0;
      
};

 #define VOID void
 typedef char CHAR;
 typedef short SHORT;
 typedef long LONG;
 typedef unsigned char   u_char;
 typedef unsigned short  u_short;
 typedef unsigned int    u_int;
 typedef unsigned long   u_long;
 typedef unsigned __int64 u_int64;
 #include <winsock2.h>
 #include <windows.h>  
 #include <stdio.h> 
  
extern int dgesvd_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
//extern int zgesvd_(char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, integer *);
extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
extern int nullSpace_(doublecomplex *, integer *, integer *, doublecomplex *, integer *);


 
int c_inter_xBase_(doublereal *a, doublereal *b, integer *n, integer *m, doublereal *gamma_real, doublereal *gamma_imag, integer *ncp, doublereal *u0, doublereal *z, doublereal *s_real, doublereal *s_imag, integer *rankB) 
{
   static integer c1 = 1;
   
   doublecomplex *gamma;
   doublecomplex *s;
   doublecomplex *c2;
   doublecomplex *ac;
   doublecomplex *sr;
   doublereal *workdsvd;
//   doublecomplex *workzsvd;
   doublecomplex *u;
   doublecomplex *v;
   doublecomplex *sigma;
   
   
   doublereal *sigmaB;
   doublereal *uB;
   doublecomplex *u1T;
   doublecomplex *u1T_mul_A;
   doublecomplex *u1T_mod;//=[u1T;u1T_mul_A]
   doublereal *vBt;
//   doublereal *rwork;
   
   doublecomplex nc={0.0,0.0};
   doublecomplex ic={1.0,0.0};
   doublecomplex normc={0.0,0.0};
     
   doublereal *work;
     
   integer nn=*n;
   integer mm=*m;
   integer nm=min(nn,mm);
   integer nncp=*ncp;
   integer nre=nn-2*nncp;
   integer i;
   integer ii;
   integer iii;
   integer iv;
   integer rowc2;   
   integer k;
   integer info;
   integer lworkdsvd=-1;
   integer lworkzsvd=-1;
   integer rrankB;
   integer rows_u1T;
   integer minr;
//   integer trowuB2;
   integer rank;
   integer nullity;
   integer nullitysr;
   
   char *all=\"A\";
   char *conj=\"C\";
   char *no=\"N\";
   char *fro=\"F\";
   char mess[200];
   
   FILE *fileptr;  
    
   fileptr = fopen(\"test.txt\",\"w\"); 
    
  
   gamma = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex));      
   s = (doublecomplex *) malloc((nn*(nn-nncp)*mm+1)*sizeof(doublecomplex));      
   ac = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
   
   sigmaB = (doublereal *) malloc((nm+1)*sizeof(doublereal));      
   uB = (doublereal *) malloc((nn*nn+1)*sizeof(doublereal));      
   vBt = (doublereal *) malloc((mm*mm+1)*sizeof(doublereal));      

        
   for(i=0;i<nn;i++)
   {
     gamma[i].r=gamma_real[i];
     gamma[i].i=gamma_imag[i];
   }
   for(i=0;i<mm*mm;i++)
   {
     z[i] = 0.0;
   }        
   for(i=0;i<nn*nn;i++)
   {
     ac[i].r = -a[i];
     ac[i].i= 0.0;
   }        

//begin decompostion of B   
   workdsvd = (doublereal *) malloc((max(3*nm+max(mm,nn),5*min(mm,nn)-4)+1)*sizeof(doublereal));
   dgesvd_(all, all, n, m, b, n, sigmaB, uB, n, vBt, m, workdsvd,  &lworkdsvd, &info);
   lworkdsvd=(int)(workdsvd[0]);
   free(workdsvd);
   workdsvd = (doublereal *) malloc((lworkdsvd+1)*sizeof(doublereal));
   dgesvd_(all, all, n, m, b, n, sigmaB, uB, n, vBt, m, workdsvd,  &lworkdsvd, &info);
   rrankB = 0;
   i = nm;
   while(i>0)
   {
    if(sigmaB[i-1]>1e-10)
    {
      rrankB = i;
      i = 0;
    }
     i = i - 1;
   }//end while;
   *rankB = rrankB;
   
   
   rows_u1T=nn-rrankB;
   u1T = (doublecomplex *) malloc((nn*rows_u1T+1)*sizeof(doublecomplex)); 
  //   c = (doublecomplex *) malloc((nn*rows_u1T+1)*sizeof(doublecomplex)); 
   u1T_mul_A = (doublecomplex *) malloc((nn*rows_u1T+1)*sizeof(doublecomplex)); 
   u1T_mod = (doublecomplex *) malloc((2*nn*rows_u1T+1)*sizeof(doublecomplex)); 

   for(i=0; i<rrankB; i++)
     for(ii=0; ii<rrankB; ii++)
       z[i*rrankB+ii] = vBt[ii*rrankB+i]/sigmaB[i];


   for(i=0;i<nn*rrankB;i++)
     u0[i] = uB[i];
   for(i=0;i<rows_u1T;i++)
     for(ii=0;ii<nn;ii++)
     {
       u1T[ii*rows_u1T+i].r = uB[(i+rrankB)*nn+ii];
       u1T[ii*rows_u1T+i].i = 0.0;
     }
//end decompostion of B   

fprintf(fileptr,\"uB= \\n\");
for(i=0;i<nn;i++)
{
  for(ii=0;ii<nn;ii++)
    fprintf(fileptr,\"%f,  \", uB[ii*nn+i]);
  fprintf(fileptr,\"\\n\");
}

// sprintf(mess,\"u1T_0=%f, u1T_1=%f, u1T_2=%f, u1T_3=%f, u1T_4=%f, u1T_5=%f \\n\",u1T[0], u1T[1], u1T[2], u1T[3], u1T[4], u1T[5]);
// MessageBoxA(NULL,mess,\"u1T\",MB_OK); 
 
 
// begin calculation of Sr      
  if(2*rrankB-nn>0 && nncp>0)
  {
    zgemm_(no, no, &rows_u1T, n, n, &ic, u1T, &rows_u1T, ac, n, &nc, u1T_mul_A, n);
    for(i=0;i<nn;i++)
    {
      for(ii=0;ii<rows_u1T;ii++)
        u1T_mod[ii*rows_u1T+i] = u1T[ii*rows_u1T+i];
      for(ii=0;ii<rows_u1T;ii++)
        u1T_mod[rows_u1T*nn+ii*rows_u1T+i] = u1T_mul_A[ii*rows_u1T+i];
    }  
    minr=min(nn,2*rows_u1T);
//    trowuB2t=2*rows_u1T;
//    sigma=(doublecomplex *) malloc((minr+1)*sizeof(doublecomplex)); 
    v=(doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));
  
    nullSpace_(u1T_mod, &rows_u1T, n, v, &nullity);
    nullitysr=nullity;
    rank = nn-nullitysr;
   
    sr = (doublecomplex *) malloc((nn*max(0,nullitysr*nn+1))*sizeof(doublecomplex)); 
  
    for(i=0;i<nullitysr;i++)
      for(ii=0;ii<nn;ii++)
        sr[i*nn+ii] = v[(i+rank)*nn+ii];
    
    c2=(doublecomplex *) malloc((nn*(nullitysr+rows_u1T)+1)*sizeof(doublecomplex));    
    
//  free(sigma);
//  free(u);
  free(v);
//  free(rwork);
//  free(workdzsvd);
}          
// end calculation of Sr


 fprintf(fileptr,\"\\n\\n sr= \\n\");
 for(i=0;i<nn;i++)
 {
   for(ii=0;ii<nullitysr;ii++)
     fprintf(fileptr,\"%f,  \", sr[ii*nn+i].r);
   fprintf(fileptr,\"\\n\");
 }


//begin computation of the nullspaces, i.e. the bases of the eigenvectors
//   v=(doublecomplex *) malloc((nn*(rows_u1T+nullitysr)+1)*sizeof(doublecomplex)); 
   v=(doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex)); 
   for(i=0;i<nn - nncp;i++)
   {

      for(ii=0;ii<nn;ii++)
      {
        ac[ii*nn+ii].r = -a[ii*nn+ii] + gamma[i].r;
        ac[ii*nn+ii].i = gamma[i].i;
      }
// fprintf(fileptr,\"\\n\\n ac= \\n\");
// for(ii=0;ii<nn;ii++)
// {
//   for(iii=0;iii<nn;iii++)
//     fprintf(fileptr,\"%f,  \", ac[iii*nn+ii].r);
//   fprintf(fileptr,\"\\n\");
// }
//       
      
      zgemm_(no, no, &rows_u1T, n, n, &ic, u1T, &rows_u1T, ac, n, &nc, u1T_mul_A, &rows_u1T);
      
//  fprintf(fileptr,\"\\n\\n u1T_mul_A= \\n\");
//  for(ii=0;ii<rows_u1T;ii++)
//  {
//    for(iii=0;iii<nn;iii++)
//      fprintf(fileptr,\"%f,  \", u1T_mul_A[iii*rows_u1T+ii].r);
//    fprintf(fileptr,\"\\n\");
//  }
    

        if(i>=nre && 2*rrankB-nn>0)
        {
          for(ii=0;ii<rows_u1T;ii++)
            for(iii=0;iii<nn;iii++)
              c2[ii*(rows_u1T+nullitysr)+iii] = u1T_mul_A[ii*rows_u1T+iii];
          for(ii=0;ii<nullitysr;ii++)
            for(iii=0;iii<nn;iii++)
            {
              c2[ii*(nullitysr)+rows_u1T+iii] = sr[iii*nn+ii];
//              c2[ii*(nullitysr)+rows_u1T+iii].i = -sr[iii*nn+ii].i;
            }
          rowc2=rows_u1T;
          nullSpace_(c2, &rowc2, n, v, &nullity);
        }
        else
          nullSpace_(u1T_mul_A, &rows_u1T, n, v, &nullity);
          
        for(ii=0;ii<rrankB;ii++)
          for(iii=0;iii<nn;iii++)
          {
            s[ii*nn+iii+i*nn*rrankB].r = v[ii*nn+iii+nn*(nn-rrankB)].r;
            s[ii*nn+iii+i*nn*rrankB].i = v[ii*nn+iii+nn*(nn-rrankB)].i;
          }
          
fprintf(fileptr,\"\\n\\n v= \\n\");
for(ii=0;ii<nn;ii++)
{
  for(iii=0;iii<nn;iii++)
    fprintf(fileptr,\"%f,  \", v[iii*nn+ii].r);
  fprintf(fileptr,\"\\n\");
}          
          
          
   }// end for i;
   
fprintf(fileptr,\"\\n\\n S= \\n\");
for(i=0;i<nn;i++)
{
  for(ii=0;ii<(nn-nncp)*rrankB;ii++)
    fprintf(fileptr,\"%f,  \", s[ii*nn+i].r);
  fprintf(fileptr,\"\\n\");
}

  fprintf(fileptr,\"\\nvor free v\");

free(v);   
//end computation of the nullspaces, i.e. the bases of the eigenvectors
   fprintf(fileptr,\"\\nnach free v\");
    for(i=0;i<nn*(nn-nncp)*mm;i++)
    {
      s_real[i]=s[i].r;
      s_imag[i]=s[i].i;
    }
  fprintf(fileptr,\"\\nvor free sr\");
   if(2*rrankB-nn>0 && nncp>0)
   {
     free(sr);
     free(c2); 
   }
   
  fprintf(fileptr,\"\\nnach free sr\");
   
   fclose(fileptr);
   free(gamma);
   free(s);
   free(ac);
  fprintf(fileptr,\"\\nnach free ac\");
   
   free(sigmaB);
   free(uB);
   free(u1T);
//   free(c);
   free(u1T_mod);
   free(u1T_mul_A);
   free(vBt);
   free(workdsvd);
//   free(workzsvd);
//   free(rwork);
  return 0;
}", Library={"zlapack"});

end wrapper_xBase;

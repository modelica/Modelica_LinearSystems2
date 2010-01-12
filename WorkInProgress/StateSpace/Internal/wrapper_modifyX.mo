within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function wrapper_modifyX
  "Contains a C sub routine of robust pole assignment to modify the eigenvector matrix X according to Kautsky algorithm"
  input Real X_real[:,size(X_real, 1)] "Eigenvector matrix, real part";
  input Real X_imag[size(X_real, 1),size(X_real, 2)]
    "Eigenvector matrix, imaginary part";
  input Integer n "Order of X";
  input Real S_real[size(X_real, 1),:] "Eigenvector bases, real part";
  input Real S_imag[size(S_real, 1),size(S_real, 2)]
    "Eigenvector bases, imaginary part";
  input Integer m
    "Rank of the system input matrix B; S_real and S_imag must have n*m columns";
  input Integer ncp "number of complex pairs";
  input Integer steps "Number of iterations";
  input Boolean IniX=false "Initial values of X are provided";

  output Real Xm_real[size(X_real, 1),size(X_real, 2)]=X_real;
  output Real Xm_imag[size(X_imag, 1),size(X_imag, 2)]=X_imag;

protected
  Integer inix=if IniX then 1 else 0;

external"FORTRAN 77" c_inter_modifyX2(
    Xm_real,
    Xm_imag,
    n,
    S_real,
    S_imag,
    m,
    ncp,
    steps,
    inix);

  annotation (Include="
  #include<f2c.h>
 

extern  int zgeqrf_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, integer *);
extern  int zungqr_(integer *, integer *, integer *, doublecomplex *, integer *,doublecomplex *, doublecomplex *, integer *, integer *);
extern  int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *, doublereal *);
extern void zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
extern int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);


int c_inter_modifyX2_(doublereal *x_real, doublereal *x_imag, integer *n, doublereal *s_real, doublereal *s_imag, integer *m, integer *ncp, integer *steps, integer *inix) 
{
   static integer c1 = 1;
   doublecomplex *x;
   doublecomplex *xj;
   doublecomplex *s;
   doublecomplex *ss;
   doublecomplex *tauqrf;
   doublecomplex *workqrf;
   doublecomplex *workgqr;
   doublecomplex *y;
   doublecomplex *qx;
   
   doublecomplex nc={0.0,0.0};
   doublecomplex ic={1.0,0.0};
   doublecomplex normc={0.0,0.0};
     
   doublereal norm=2;
   doublereal *work;
     
   integer nn=*n;
   integer mm=*m;
   integer nncp=*ncp;
   integer nre=nn-2*nncp;
   integer i;
   integer ii;
   integer iii;
   integer k;
   integer *info;
   integer info2;
   integer lworkqrf=-1;
   integer lworkgqr=-1;
   
   char *conj=\"C\";
   char *no=\"N\";
   char *fro=\"F\";
   
   if(nn>mm)
   {
   x = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
   xj = (doublecomplex *) malloc((nn*nn+1)*sizeof(doublecomplex));      
   s = (doublecomplex *) malloc((nn*nn*mm+1)*sizeof(doublecomplex));      
   ss = (doublecomplex *) malloc((nn*mm+1)*sizeof(doublecomplex));      
   y = (doublecomplex *) malloc((mm+1)*sizeof(doublecomplex));      
   qx = (doublecomplex *) malloc((nn+1)*sizeof(doublecomplex));      
   tauqrf = (doublecomplex *) malloc((nn)*sizeof(doublecomplex));      

   work = (doublereal *) malloc((2)*sizeof(doublereal));
   workqrf = (doublecomplex *) malloc((2*nn+1)*sizeof(doublecomplex));
   
   zgeqrf_(n, n, xj, n, tauqrf, workqrf, &lworkqrf, &info2);      
 
   lworkqrf=(int)(workqrf[0].r);
   free(workqrf);
   workqrf = (doublecomplex *) malloc((lworkqrf+1)*sizeof(doublecomplex));

   workgqr = (doublecomplex *) malloc((2*nn+1)*sizeof(doublecomplex));
   zungqr_(n, n, n, xj, n, tauqrf, workgqr, &lworkgqr, &info2);
   lworkgqr=(int)(workgqr[0].r);
   free(workgqr);
   workgqr = (doublecomplex *) malloc((lworkgqr+1)*sizeof(doublecomplex));
       
    for(i=0;i<nn*nn*mm;i++)
    {
      s[i].r=s_real[i];
      s[i].i=s_imag[i];
    }     

   if(*inix==1)
   { 
     for(i=0;i<nn*nn;i++)
     {
       x[i].r=x_real[i];
       x[i].i=x_imag[i];
     }
   }
   else
   {
     for(i=0;i<nn - nncp;i++)
       zcopy_(n, &s[mm*i*nn], &c1, &x[i*nn], &c1);  
     for(i=0;i<nncp;i++)
       for(ii=0;ii<nn;ii++)
       {
          x[ii+ (nre + nncp + i)*nn].r = x[ii+ (nre + i)*nn].r;
          x[ii+ (nre + nncp + i)*nn].i = -x[ii+ (nre + i)*nn].i;
       }
   }
   
   if(mm>1)
   {
    for(i=0;i<nn*nn;i++)
      xj[i]=nc;   
  
  k = 0;
  while(k < *steps)
  {
    k = k + 1;
 

    for(i=1; i<=nn-nncp;i++)
    {
      if(i==1)
      {
        for(ii=0;ii<nn*(nn-1);ii++)
          xj[ii] = x[ii+nn];
      }
      else
      {
        for(ii=0;ii<(i-1)*nn;ii++)
          xj[ii] = x[ii];
        for(ii=i*nn;ii<nn*nn;ii++)
          xj[ii-nn] = x[ii];
      }//end if
      
      for(ii=0;ii<nn;ii++)
        xj[nn*(nn-1)+ii]=nc;   
        
      zgeqrf_(n, n, xj, n, tauqrf, workqrf, &lworkqrf, &info2);  
      zungqr_(n, n, n, xj, n, tauqrf, workgqr, &lworkgqr, &info2);
      
      for(ii=0;ii<mm*nn;ii++)
        ss[ii] = s[(i-1)*mm*nn+ii];
      for(ii=0;ii<nn;ii++)
        qx[ii] = xj[nn*(nn-1)+ii];

      zgemv_(conj, n, m, &ic, ss, n, qx, &c1, &nc, y, &c1);
      norm=zlange_(fro, m, &c1, y, n, work);

      normc.r=1/norm;
      normc.i=0.0;
      zgemv_(no, n, m, &normc, ss, n, y, &c1, &nc, qx, &c1);//y=qx
      
//       zdotc_(&normc, m, qx, &c1, qx, &c1);//qx^H*qx
//       if(i>nncp && normc.r*normc.r+normc.i*normc.i>0.81)
//       {
//            idx = 1 + rem(k, mm-nullitySr);
//            qx = (qx + S[:, i*nn+idx])/sqrt(2);
//        }

     for(ii=0;ii<nn;ii++)
       x[(i-1)*nn+ii] = qx[ii];
       
      if(i>nre)
      {
        for(ii=0;ii<nn;ii++)
        {
          x[(i+nncp-1)*nn+ii].r = qx[ii].r;
          x[(i+nncp-1)*nn+ii].i = -qx[ii].i;
        }
      }//end if;

     }// end for i
}//end while;
}// end if mm>1

   for(i=0;i<nn*nn;i++)
   {
     x_real[i]=x[i].r;
     x_imag[i]=x[i].i;
   }


  
   free(x);
   free(xj);
   free(s);
   free(ss);
   free(y);
   free(qx);
   free(workqrf);
   free(tauqrf);
   free(workgqr);
   free(work);
   }//end if nn>mm
   else
   {
     for(i=0;i<nn*nn;i++)
     {
       x_real[i]=0.0;
       x_imag[i]=0.0;
     }
     for(i=0;i<nre;i++)
     {
       x_real[i*nn+i]=1.0;
     }
     for(i=0;i<nncp;i++)
     {
         x_real[nn*(nre+2*i)+nre+2*i] = 0.5;
         x_real[nn*(nre+2*i)+nre+2*i+1] = 0.5;
         x_imag[nn*(nre+2*i+1)+nre+2*i+1] = 0.5;
         x_imag[nn*(nre+2*i+1)+nre+2*i] = -0.5;
     }

     
   }
   
  return 0;
}", Library={"zlapack"});

end wrapper_modifyX;

within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function multiplyRealMatrices3cCall "Some case studies to complex matrices"

  input Integer n=50;
  input Integer nr=20;
  input Real C1[n,n]=fill(1,n,n);
  output Real C2[n,n]=C1;

external"FORTRAN 77" c_dgemm(C1, C2, n, nr);
  annotation (Include="
  #include<f2c.h>
  #include <stdio.h>


extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);

int c_dgemm_(doublereal *c1, doublereal *c2, integer *n, integer *nr)
{

   static doublereal ndd=0.0;
   static doublereal idd=1.0;
   static integer incx=1;
   doublereal *c3;

   integer nn=*n;
   integer nxn=nn*nn;
   integer nnr=*nr;
   integer i;
   integer ii;

   char *no=\"N\";

 //   FILE *fileptr;
 //   fileptr = fopen(\"test.txt\",\"w\");

   c3 = (doublereal *) malloc((nn*nn+1)*sizeof(doublereal));


   for(i=0;i<nnr;i++)
   {
     dgemm_(no, no, n, n, n, &idd, c1, n, c2, n, &ndd, c3, n);
     dcopy_(&nxn, c3, &incx, c2, &incx);
   }


//    for(i=0;i<nn;i++)
//    {
//      for(ii=0;ii<nn;ii++)
//        fprintf(fileptr,\"%f,  \", c3[ii*nn+i]);
//      fprintf(fileptr,\"\\n\");
//    }
//   fclose(fileptr);

   free(c3);
   return 0;
}", Library={"lapack"});

end multiplyRealMatrices3cCall;

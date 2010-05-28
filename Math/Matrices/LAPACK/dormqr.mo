within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dormqr
  "overwrites the general real M-by-N matrix C with Q * C or C * Q or Q' * C or C * Q', where Q is an orthogonal matrix of a QR factorization as returned by dgeqrf"

  input Real C[:,:];
  input Real A[:,:];
  input Real tau[:];
  input String side="L";
  input String trans="N";

  output Real Cout[size(C, 1),size(C, 2)]=C
    "contains Q*C or Q**T*C or C*Q**T or C*Q";

  output Integer info;
protected
  Integer m=size(C, 1);
  Integer n=size(C, 2);
  Integer k=if side == "L" then m else n;
  Integer lda=if side == "L" then max(1, m) else max(1, n);
  Integer ldc=max(1, m);
  Integer lwork=if side == "L" then max(1, n) else max(1, m);
  Real work[lwork];

external "Fortran 77" dormqr(
    side,
    trans,
    m,
    n,
    k,
    A,
    lda,
    tau,
    Cout,
    ldc,
    work,
    lwork,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="   Purpose  
   =======  
 
   DORMQR overwrites the general real M-by-N matrix C with  
 
                   SIDE = 'L'     SIDE = 'R'  
   TRANS = 'N':      Q * C          C * Q  
   TRANS = 'T':      Q**T * C       C * Q**T  
 
   where Q is a real orthogonal matrix defined as the product of k  
   elementary reflectors  
 
         Q = H(1) H(2) . . . H(k)  
 
   as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N  
   if SIDE = 'R'.  
 
   Arguments  
   =========  
 
   SIDE    (input) CHARACTER*1  
           = 'L': apply Q or Q**T from the Left;  
           = 'R': apply Q or Q**T from the Right.  
 
   TRANS   (input) CHARACTER*1  
           = 'N':  No transpose, apply Q;  
           = 'T':  Transpose, apply Q**T.  
 
   M       (input) INTEGER  
           The number of rows of the matrix C. M >= 0.  
 
   N       (input) INTEGER  
           The number of columns of the matrix C. N >= 0.  
 
   K       (input) INTEGER  
           The number of elementary reflectors whose product defines  
           the matrix Q.  
           If SIDE = 'L', M >= K >= 0;  
           if SIDE = 'R', N >= K >= 0.  
 
   A       (input) DOUBLE PRECISION array, dimension (LDA,K)  
           The i-th column must contain the vector which defines the  
           elementary reflector H(i), for i = 1,2,...,k, as returned by  
           DGEQRF in the first k columns of its array argument A.  
           A is modified by the routine but restored on exit.  
 
   LDA     (input) INTEGER  
           The leading dimension of the array A.  
           If SIDE = 'L', LDA >= max(1,M);  
           if SIDE = 'R', LDA >= max(1,N).  
 
   TAU     (input) DOUBLE PRECISION array, dimension (K)  
           TAU(i) must contain the scalar factor of the elementary  
           reflector H(i), as returned by DGEQRF.  
 
   C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)  
           On entry, the M-by-N matrix C.  
           On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.  
 
   LDC     (input) INTEGER  
           The leading dimension of the array C. LDC >= max(1,M).  
 
   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)  
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.  
 
   LWORK   (input) INTEGER  
           The dimension of the array WORK.  
           If SIDE = 'L', LWORK >= max(1,N);  
           if SIDE = 'R', LWORK >= max(1,M).  
           For optimum performance LWORK >= N*NB if SIDE = 'L', and  
           LWORK >= M*NB if SIDE = 'R', where NB is the optimal  
           blocksize.  
 
           If LWORK = -1, then a workspace query is assumed; the routine  
           only calculates the optimal size of the WORK array, returns  
           this value as the first entry of the WORK array, and no error  
           message related to LWORK is issued by XERBLA.  
 
   INFO    (output) INTEGER  
           = 0:  successful exit  
           < 0:  if INFO = -i, the i-th argument had an illegal value  
 
   =====================================================================  "));
end dormqr;

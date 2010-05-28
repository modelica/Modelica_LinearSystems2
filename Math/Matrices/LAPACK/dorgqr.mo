within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dorgqr
  "generates a real orthogonal matrix Q which is defined as the product of elementary reflectors, as returned by DGEQRF"

  input Real Q[:,:] "Orthogonal matrix as a result of elementary reflectors";
  input Real tau[:] "scalar factors of the elementary reflectors";

  output Real Qout[size(Q,1),size(Q,2)]=Q;
  output Integer info;
protected
  Integer m=size(Q, 1);
  Integer n=size(Q, 2);
  Integer k=size(tau, 1);
  Integer lda=max(1, m);
  Integer lwork=3*max(1, m);
  Real work[lwork];

external "Fortran 77" dorgqr(
    m,
    n,
    k,
    Qout,
    lda,
    tau,
    work,
    lwork,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="   Purpose  
   =======  
 
   DORGQR generates an M-by-N real matrix Q with orthonormal columns,  
   which is defined as the first N columns of a product of K elementary  
   reflectors of order M  
 
         Q  =  H(1) H(2) . . . H(k)  
 
   as returned by DGEQRF.  
 
   Arguments  
   =========  
 
   M       (input) INTEGER  
           The number of rows of the matrix Q. M >= 0.  
 
   N       (input) INTEGER  
           The number of columns of the matrix Q. M >= N >= 0.  
 
   K       (input) INTEGER  
           The number of elementary reflectors whose product defines the  
           matrix Q. N >= K >= 0.  
 
   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)  
           On entry, the i-th column must contain the vector which  
           defines the elementary reflector H(i), for i = 1,2,...,k, as  
           returned by DGEQRF in the first k columns of its array  
           argument A.  
           On exit, the M-by-N matrix Q.  
 
   LDA     (input) INTEGER  
           The first dimension of the array A. LDA >= max(1,M).  
 
   TAU     (input) DOUBLE PRECISION array, dimension (K)  
           TAU(i) must contain the scalar factor of the elementary  
           reflector H(i), as returned by DGEQRF.  
 
   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)  
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.  
 
   LWORK   (input) INTEGER  
           The dimension of the array WORK. LWORK >= max(1,N).  
           For optimum performance LWORK >= N*NB, where NB is the  
           optimal blocksize.  
 
           If LWORK = -1, then a workspace query is assumed; the routine  
           only calculates the optimal size of the WORK array, returns  
           this value as the first entry of the WORK array, and no error  
           message related to LWORK is issued by XERBLA.  
 
   INFO    (output) INTEGER  
           = 0:  successful exit  
           < 0:  if INFO = -i, the i-th argument has an illegal value  
 
   =====================================================================  
"));
end dorgqr;

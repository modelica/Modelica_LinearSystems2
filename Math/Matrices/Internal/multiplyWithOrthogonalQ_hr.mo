within Modelica_LinearSystems2.Math.Matrices.Internal;
function multiplyWithOrthogonalQ_hr
  "Overwrites the general real M-by-N matrix C with Q * C or C * Q or Q' * C or C * Q', depending on inputs trans and side"

  input Real C[:,:];
  input Real A[:,:];
  input Real tau[size(A, 2) - 1];
  input String side="L";
  input String trans="N";
  input Integer ilo=1
    "lowest index where the original matrix had been Hessenbergform";
  input Integer ihi=size(A, 2)
    "highest index where the original matrix had been Hessenbergform";
  output Real Cout[size(C, 1),size(C, 2)]=C
    "contains the Hessenberg form in the upper triangle and the first subdiagonal and below the first subdiagonal it contains the elementary reflectors which represents (with array tau) as a product the orthogonal matrix Q";

  output Integer info;
protected
  Integer m=size(C, 1);
  Integer n=size(C, 2);
  Integer lda=max(1, size(A, 2));
  Integer ldc=max(1, m);
  Integer lwork=2*size(A, 2);
  Real work[lwork];

external "Fortran 77" dormhr(
      side,
      trans,
      m,
      n,
      ilo,
      ihi,
      A,
      lda,
      tau,
      Cout,
      ldc,
      work,
      lwork,
      info)
          annotation(Library = {"lapack"});

  annotation (Documentation(info="   Purpose  
    -- LAPACK routine (version 3.0) --  
      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,  
      Courant Institute, Argonne National Lab, and Rice University  
      June 30, 1999  
 
      .. Scalar Arguments ..  
      ..  
      .. Array Arguments ..  
      ..  
 
   Purpose  
   =======  
 
     DORMHR overwrites the general real M-by-N matrix C with  
 
                   SIDE = 'L'     SIDE = 'R'  
   TRANS = 'N':      Q * C          C * Q  
   TRANS = 'T':      Q**T * C       C * Q**T  
 
   where Q is a real orthogonal matrix of order nq, with nq = m if  
   SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of  
   IHI-ILO elementary reflectors, as returned by DGEHRD:  
 
   Q = H(ilo) H(ilo+1) . . . H(ihi-1).  
 
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
 
   ILO     (input) INTEGER  
   IHI     (input) INTEGER  
           ILO and IHI must have the same values as in the previous call  
           of DGEHRD. Q is equal to the unit matrix except in the  
           submatrix Q(ilo+1:ihi,ilo+1:ihi).  
           If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and  
           ILO = 1 and IHI = 0, if M = 0;  
           if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and  
           ILO = 1 and IHI = 0, if N = 0.  
 
   A       (input) DOUBLE PRECISION array, dimension  
                                (LDA,M) if SIDE = 'L'  
                                (LDA,N) if SIDE = 'R'  
           The vectors which define the elementary reflectors, as  
           returned by DGEHRD.  
 
   LDA     (input) INTEGER  
           The leading dimension of the array A.  
           LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.  
 
   TAU     (input) DOUBLE PRECISION array, dimension  
                                (M-1) if SIDE = 'L'  
                                (N-1) if SIDE = 'R'  
           TAU(i) must contain the scalar factor of the elementary  
           reflector H(i), as returned by DGEHRD.  
 
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
end multiplyWithOrthogonalQ_hr;

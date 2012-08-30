within Modelica_LinearSystems2.Math.Matrices;
function orthogonalQ
  "Generates a real orthogonal matrix Q defined as the product of IHI-ILO elementary reflectors"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real A[:,size(A, 1)];
  input Real tau[size(A, 1) - 1] "Scalar factors of the elementary reflectors";
  input Integer ilo=1
    "Lowest index where the original matrix had been Hessenbergform - ilo must have the same value as in the previous call of DGEHRD";
  input Integer ihi
    "Highest index where the original matrix had been Hessenbergform  - ihi must have the same value as in the previous call of DGEHRD";

  output Real Q[size(A, 1),size(A, 2)]
    "Orthogonal matrix as a result of elementary reflectors";
  output Integer info;

algorithm
  (Q,info) := LAPACK.dorghr(
    A,
    ilo,
    ihi,
    tau);
  annotation (Documentation(info="<html>
<p>
This function generates a real orthogonal matrix Q which is defined as the product
of IHI-ILO elementary reflectors of order N, as returned by DGEHRD.
</p>

<pre>
Lapack documentation:

   Purpose
   =======

   DORGHR generates a real orthogonal matrix Q which is defined as the
   product of IHI-ILO elementary reflectors of order N, as returned by
   DGEHRD:

   Q = H(ilo) H(ilo+1) . . . H(ihi-1).

   Arguments
   =========

   N       (input) INTEGER
           The order of the matrix Q. N &gt;= 0.

   ILO     (input) INTEGER
   IHI     (input) INTEGER
           ILO and IHI must have the same values as in the previous call
           of DGEHRD. Q is equal to the unit matrix except in the
           submatrix Q(ilo+1:ihi,ilo+1:ihi).
           1 &lt;= ILO &lt;= IHI &lt;= N, if N &gt; 0; ILO=1 and IHI=0, if N=0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the vectors which define the elementary reflectors,
           as returned by DGEHRD.
           On exit, the N-by-N orthogonal matrix Q.

   LDA     (input) INTEGER
           The leading dimension of the array A. LDA &gt;= max(1,N).

   TAU     (input) DOUBLE PRECISION array, dimension (N-1)
           TAU(i) must contain the scalar factor of the elementary
   ), as returned by DGEHRD.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK. LWORK &gt;= IHI-ILO.
           For optimum performance LWORK &gt;= (IHI-ILO)*NB, where NB is
           the optimal blocksize.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           &lt; 0:  if INFO = -i, the i-th argument had an illegal value

   =====================================================================

</pre>
</html>"));
end orthogonalQ;

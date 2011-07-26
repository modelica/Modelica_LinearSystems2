within Modelica_LinearSystems2.WorkInProgress.Optimizer;
function constrainedLeastSquares
  "Solve linearly constrained least squares problem with equality and inequality constraints (interface to Fortran function DLSEI)"

  extends Modelica.Icons.Function;
  input Real W[:, :] "System matrix = [E,f; A,b; G,h]";
  input Integer n_E=0
    "Number of equations to be exactly satisfied (number of rows of matrix E)";
  input Integer n_G=0
    "Number of inequality equations (number of rows of matrix G)";
  input Real options[:]={1.0} "Option vector (see info)";
  output Real x[size(W, 2) - 1] "Solution vector if mode = 0 or 1";
  output Integer mode
    "Status after completion (mode = 0/1 is successful, 2/3/4 is error)";
  output Real residueNorm_A=0 "Length of b-A*x, if mode = 0 or 1";
  output Real residueNorm_E=0 "Length of f-E*x, if mode = 0 or 1";
  output Real W_out[size(W, 1), size(W, 2)]=W
    "Covariance matrix, if requested by options vector";
protected
  Integer n_x=size(W, 2) - 1 "Dimension of solution vector x";
  Integer n_A=size(W, 1) - n_E - n_G
    "Number of equations to be approximately satisfied";
  Integer IP1=2*(n_E + n_x) + (n_G + 2)*(n_x + 7) + max(n_A + n_G, n_x)
    "Size of real work array";
  Integer IP2=n_G + 2*n_x + 2 "Size of integer work array";
  Real WS[IP1] "Real work array";
  Integer IP[IP2] "Integer work array";
external "Fortran 77" dlsei(W_out, size(W, 1), n_E, n_A, n_G, n_x, options, x,
     residueNorm_E, residueNorm_A, mode, WS, IP, IP1, IP2);

  annotation (
    Include="#include <dlsei.h>",
    Library="dlsei",
    Documentation(info="<HTML>
<p>
This function solves a least squares problem with 
linear equality and linear inequality constraints.
The problem to be solved is mathematically formulated as:
</p>
<pre>
   min |<b>A</b>*<b>x</b> - <b>b</b>|^2
   subject to   <b>E</b>*<b>x</b>  = <b>f</b>;
                <b>G</b>*<b>x</b> >= <b>h</b>;
</pre>
<p>
Note, that the inequality relation is treated element wise.
Remarks:
</p>
<ul>
<li> If there are infinitely many solutions for <b>x</b> (because
     <b>A</b> does not have full rank), the solution <b>x</b> with the
     smallest norm is selected that minimizes the problem above.</li>
<li> If the equality constraints cannot be satisfied, a
     generalized inverse solution is computed such that the
     residual vector length <b>f</b>-<b>E</b>*<b>x</b> is
     minimized.</li>
<li> The number of rows of matrices <b>A</b>, <b>E</b>, and/or <b>G</b> 
     are arbitrary.
     Especially, these matrices may have zero rows.</li>
</ul>
<p>
As input argument <b>W</b> the assembled matrices have to be
given:
</p>
<pre>
         | <b>E  f</b> |
    <b>W</b> =  | <b>A  b</b> |
         | <b>G  h</b> |
</pre>
<p>
With output argument <b>mode</b> the type of solution is
characterized:
</p>
<pre>
  mode = 0: A solution has been computed. Both equality and 
            inequality constraints are compatible and have been satisfied.
       = 1: A solution has been computed. Equality constraints are 
            contradictory. A generalized inverse solution of E*x=f was used
            to minimize the residual vector length f-E*x.
            In this sense, the solution is still meaningful.
       > 1: Error, no solution has been computed
       = 2: Inequality constraints are contradictory.
       = 3: Both equality and inequality constraints are contradictory.
       = 4: Usage error (input arguments are wrong).
</pre>


<p><b>Copyright &copy; 2004-2005, DLR Institute of Robotics and Mechatronics </b></p>

<p>
The Fortran interface of DLSEI is defined as:
</p>
<pre>
     SUBROUTINE DLSEI (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME,
    +                 RNORML, MODE, WS, IP)
***BEGIN PROLOGUE  DLSEI
***PURPOSE  Solve a linearly constrained least squares problem with
            equality and inequality constraints, and optionally compute
            a covariance matrix.
***LIBRARY   SLATEC
***CATEGORY  K1A2A, D9
***TYPE      DOUBLE PRECISION (LSEI-S, DLSEI-D)
***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
             QUADRATIC PROGRAMMING
***AUTHOR  Hanson, R. J., (SNLA)
           Haskell, K. H., (SNLA)
***DESCRIPTION
     Abstract
     This subprogram solves a linearly constrained least squares
     problem with both equality and inequality constraints, and, if the
     user requests, obtains a covariance matrix of the solution
     parameters.
     Suppose there are given matrices E, A and G of respective
     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
     respective lengths ME, MA and MG.  This subroutine solves the
     linearly constrained least squares problem
                   EX = F, (E ME by N) (equations to be exactly
                                       satisfied)
                   AX = B, (A MA by N) (equations to be
                                       approximately satisfied,
                                       least squares sense)
                   GX .GE. H,(G MG by N) (inequality constraints)
     The inequalities GX .GE. H mean that every component of the
     product GX must be .GE. the corresponding component of H.
     In case the equality constraints cannot be satisfied, a
     generalized inverse solution residual vector length is obtained
     for F-EX.  This is the minimal length possible for F-EX.
     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
     rank of the matrix E is estimated during the computation.  We call
     this value KRANKE.  It is an output parameter in IP(1) defined
     below.  Using a generalized inverse solution of EX=F, a reduced
     least squares problem with inequality constraints is obtained.
     The tolerances used in these tests for determining the rank
     of E and the rank of the reduced least squares problem are
     given in Sandia Tech. Rept. SAND-78-1290.  They can be
     modified by the user if new values are provided in
     the option list of the array PRGOPT(*).
     The user must dimension all arrays appearing in the call list..
     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
     where K=MAX(MA+MG,N).  This allows for a solution of a range of
     problems in the given working space.  The dimension of WS(*)
     given is a necessary overestimate.  Once a particular problem
     has been run, the output parameter IP(3) gives the actual
     dimension required for that problem.
     The parameters for DLSEI( ) are
     Input.. All TYPE REAL variables are DOUBLE PRECISION
     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
     ME,MA,MG,N    first dimensioning parameter equal to MDW.
                   For this discussion let us call M = ME+MA+MG.  Then
                   MDW must satisfy MDW .GE. M.  The condition
                   MDW .LT. M is an error.
                   The array W(*,*) contains the matrices and vectors
                                  (E  F)
                                  (A  B)
                                  (G  H)
                   in rows and columns 1,...,M and 1,...,N+1
                   respectively.
                   The integers ME, MA, and MG are the
                   respective matrix row dimensions
                   of E, A and G.  Each matrix has N columns.
     PRGOPT(*)    This real-valued array is the option vector.
                  If the user is satisfied with the nominal
                  subprogram features set
                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
                  Otherwise PRGOPT(*) is a linked list consisting of
                  groups of data of the following form
                  LINK
                  KEY
                  DATA SET
                  The parameters LINK and KEY are each one word.
                  The DATA SET can be comprised of several words.
                  The number of items depends on the value of KEY.
                  The value of LINK points to the first
                  entry of the next group of data within
                  PRGOPT(*).  The exception is when there are
                  no more options to change.  In that
                  case, LINK=1 and the values KEY and DATA SET
                  are not referenced.  The general layout of
                  PRGOPT(*) is as follows.
               ...PRGOPT(1) = LINK1 (link to first entry of next group)
               .  PRGOPT(2) = KEY1 (key to the option change)
               .  PRGOPT(3) = data value (data value for this change)
               .       .
               .       .
               .       .
               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
               .                       next group)
               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
               .  PRGOPT(LINK1+2) = data value
               ...     .
               .       .
               .       .
               ...PRGOPT(LINK) = 1 (no more options to change)
                  Values of LINK that are nonpositive are errors.
                  A value of LINK .GT. NLINK=100000 is also an error.
                  This helps prevent using invalid but positive
                  values of LINK that will probably extend
                  beyond the program limits of PRGOPT(*).
                  Unrecognized values of KEY are ignored.  The
                  order of the options is arbitrary and any number
                  of options can be changed with the following
                  restriction.  To prevent cycling in the
                  processing of the option array, a count of the
                  number of options changed is maintained.
                  Whenever this count exceeds NOPT=1000, an error
                  message is printed and the subprogram returns.
                  Options..
                  KEY=1
                         Compute in W(*,*) the N by N
                  covariance matrix of the solution variables
                  as an output parameter.  Nominally the
                  covariance matrix will not be computed.
                  (This requires no user input.)
                  The data set for this option is a single value.
                  It must be nonzero when the covariance matrix
                  is desired.  If it is zero, the covariance
                  matrix is not computed.  When the covariance matrix
                  is computed, the first dimensioning parameter
                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
                  KEY=10
                         Suppress scaling of the inverse of the
                  normal matrix by the scale factor RNORM**2/
                  MAX(1, no. of degrees of freedom).  This option
                  only applies when the option for computing the
                  covariance matrix (KEY=1) is used.  With KEY=1 and
                  KEY=10 used as options the unscaled inverse of the
                  normal matrix is returned in W(*,*).
                  The data set for this option is a single value.
                  When it is nonzero no scaling is done.  When it is
                  zero scaling is done.  The nominal case is to do
                  scaling so if option (KEY=1) is used alone, the
                  matrix will be scaled on output.
                  KEY=2
                         Scale the nonzero columns of the
                         entire data matrix.
                  (E)
                  (A)
                  (G)
                  to have length one.  The data set for this
                  option is a single value.  It must be
                  nonzero if unit length column scaling
                  is desired.
                  KEY=3
                         Scale columns of the entire data matrix
                  (E)
                  (A)
                  (G)
                  with a user-provided diagonal matrix.
                  The data set for this option consists
                  of the N diagonal scaling factors, one for
                  each matrix column.
                  KEY=4
                         Change the rank determination tolerance for
                  the equality constraint equations from
                  the nominal value of SQRT(DRELPR).  This quantity can
                  be no smaller than DRELPR, the arithmetic-
                  storage precision.  The quantity DRELPR is the
                  largest positive number such that T=1.+DRELPR
                  satisfies T .EQ. 1.  The quantity used
                  here is internally restricted to be at
                  least DRELPR.  The data set for this option
                  is the new tolerance.
                  KEY=5
                         Change the rank determination tolerance for
                  the reduced least squares equations from
                  the nominal value of SQRT(DRELPR).  This quantity can
                  be no smaller than DRELPR, the arithmetic-
                  storage precision.  The quantity used
                  here is internally restricted to be at
                  least DRELPR.  The data set for this option
                  is the new tolerance.
                  For example, suppose we want to change
                  the tolerance for the reduced least squares
                  problem, compute the covariance matrix of
                  the solution parameters, and provide
                  column scaling for the data matrix.  For
                  these options the dimension of PRGOPT(*)
                  must be at least N+9.  The Fortran statements
                  defining these options would be as follows:
                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
                  PRGOPT(2)=1 (covariance matrix key)
                  PRGOPT(3)=1 (covariance matrix wanted)
                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
                  PRGOPT(5)=5 (least squares equas.  tolerance key)
                  PRGOPT(6)=... (new value of the tolerance)
                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
                  PRGOPT(8)=3 (user-provided column scaling key)
                  CALL DCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
                    scaling factors from the user array D(*)
                    to PRGOPT(9)-PRGOPT(N+8))
                  PRGOPT(N+9)=1 (no more options to change)
                  The contents of PRGOPT(*) are not modified
                  by the subprogram.
                  The options for WNNLS( ) can also be included
                  in this array.  The values of KEY recognized
                  by WNNLS( ) are 6, 7 and 8.  Their functions
                  are documented in the usage instructions for
                  subroutine WNNLS( ).  Normally these options
                  do not need to be modified when using DLSEI( ).
     IP(1),       The amounts of working storage actually
     IP(2)        allocated for the working arrays WS(*) and
                  IP(*), respectively.  These quantities are
                  compared with the actual amounts of storage
                  needed by DLSEI( ).  Insufficient storage
                  allocated for either WS(*) or IP(*) is an
                  error.  This feature was included in DLSEI( )
                  because miscalculating the storage formulas
                  for WS(*) and IP(*) might very well lead to
                  subtle and hard-to-find execution errors.
                  The length of WS(*) must be at least
                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
                  where K = max(MA+MG,N)
                  This test will not be made if IP(1).LE.0.
                  The length of IP(*) must be at least
                  LIP = MG+2*N+2
                  This test will not be made if IP(2).LE.0.
     Output.. All TYPE REAL variables are DOUBLE PRECISION
     X(*),RNORME,  The array X(*) contains the solution parameters
     RNORML        if the integer output flag MODE = 0 or 1.
                   The definition of MODE is given directly below.
                   When MODE = 0 or 1, RNORME and RNORML
                   respectively contain the residual vector
                   Euclidean lengths of F - EX and B - AX.  When
                   MODE=1 the equality constraint equations EX=F
                   are contradictory, so RNORME .NE. 0.  The residual
                   vector F-EX has minimal Euclidean length.  For
                   MODE .GE. 2, none of these parameters is defined.
     MODE          Integer flag that indicates the subprogram
                   status after completion.  If MODE .GE. 2, no
                   solution has been computed.
                   MODE =
                   0  Both equality and inequality constraints
                      are compatible and have been satisfied.
                   1  Equality constraints are contradictory.
                      A generalized inverse solution of EX=F was used
                      to minimize the residual vector length F-EX.
                      In this sense, the solution is still meaningful.
                   2  Inequality constraints are contradictory.
                   3  Both equality and inequality constraints
                      are contradictory.
                   The following interpretation of
                   MODE=1,2 or 3 must be made.  The
                   sets consisting of all solutions
                   of the equality constraints EX=F
                   and all vectors satisfying GX .GE. H
                   have no points in common.  (In
                   particular this does not say that
                   each individual set has no points
                   at all, although this could be the
                   case.)
                   4  Usage error occurred.  The value
                      of MDW is .LT. ME+MA+MG, MDW is
                      .LT. N and a covariance matrix is
                      requested, or the option vector
                      PRGOPT(*) is not properly defined,
                      or the lengths of the working arrays
                      WS(*) and IP(*), when specified in
                      IP(1) and IP(2) respectively, are not
                      long enough.
     W(*,*)        The array W(*,*) contains the N by N symmetric
                   covariance matrix of the solution parameters,
                   provided this was requested on input with
                   the option vector PRGOPT(*) and the output
                   flag is returned with MODE = 0 or 1.
     IP(*)         The integer working array has three entries
                   that provide rank and working array length
                   information after completion.
                      IP(1) = rank of equality constraint
                              matrix.  Define this quantity
                              as KRANKE.
                      IP(2) = rank of reduced least squares
                              problem.
                      IP(3) = the amount of storage in the
                              working array WS(*) that was
                              actually used by the subprogram.
                              The formula given above for the length
                              of WS(*) is a necessary overestimate.
                              If exactly the same problem matrices
                              are used in subsequent executions,
                              the declared dimension of WS(*) can
                              be reduced to this output value.
     User Designated
     Working Arrays..
     WS(*),IP(*)              These are respectively type real
                              and type integer working arrays.
                              Their required minimal lengths are
                              given above.
***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
                 linear least squares problems with equality and
                 nonnegativity constraints, Report SAND77-0552, Sandia
                 Laboratories, June 1978.
               K. H. Haskell and R. J. Hanson, Selected algorithms for
                 the linearly constrained least squares problem - a
                 users guide, Report SAND78-1290, Sandia Laboratories,
                 August 1979.
               K. H. Haskell and R. J. Hanson, An algorithm for
                 linear least squares problems with equality and
                 nonnegativity constraints, Mathematical Programming
                 21 (1981), pp. 98-118.
               R. J. Hanson and K. H. Haskell, Two algorithms for the
                 linearly constrained least squares problem, ACM
                 Transactions on Mathematical Software, September 1982.
***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DLSI,
                    DNRM2, DSCAL, DSWAP, XERMSG
***REVISION HISTORY  (YYMMDD)
   790701  DATE WRITTEN
   890531  Changed all specific intrinsics to generic.  (WRB)
   890618  Completely restructured and extensively revised (WRB & RWC)
   890831  REVISION DATE from Version 3.2
   891214  Prologue converted to Version 4.0 format.  (BAB)
   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
   900604  DP version created from SP version.  (RWC)
   920501  Reformatted the REFERENCES section.  (WRB)
</pre>
</HTML>
", revisions="<html>
<p><b>Release Notes:</b></p>
<ul>
<li><i>Sept. 3, 2004</i>
       Martin Otter (DLR): implemented
</li>
</ul>
</html>"),
    uses(Modelica(version="2.2.2 development")));
end constrainedLeastSquares;

within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dtrsyl "DTRSYL solves the real Sylvester matrix equation"

  input Real A[:,:];
  input Real B[:,:];
  input Real C[if tranA then size(A,1) else size(A, 2),if tranB then size(B,1) else size(B, 2)];

  input Boolean tranA=false;
  input Boolean tranB=false;
  input Integer isgn=1;
  output Real X[size(C,1),size(C,2)]=C;
  output Real scale;
  output Integer info;
protected
  Integer m=if tranA then size(A,1) else size(A, 2);
  Integer n=if tranB then size(B,1) else size(B, 2);
  String trana = if tranA then "T" else "N";
  String tranb = if tranB then "T" else "N";

external "Fortran 77" dtrsyl(trana, tranb, isgn, m, n, A, m, B, n, X, m, scale, info) annotation(Library = {"lapack"});

  annotation (Documentation(info="   
   Purpose   
    =======   

    DTRSYL solves the real Sylvester matrix equation:   

       op(A)*X + X*op(B) = scale*C or   
       op(A)*X - X*op(B) = scale*C,   

    where op(A) = A or A**T, and  A and B are both upper quasi-   
    triangular. A is M-by-M and B is N-by-N; the right hand side C and   
    the solution X are M-by-N; and scale is an output scale factor, set   
    <= 1 to avoid overflow in X.   

    A and B must be in Schur canonical form (as returned by DHSEQR), that   
    is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;   
    each 2-by-2 diagonal block has its diagonal elements equal and its   
    off-diagonal elements of opposite sign.   

    Arguments   
    =========   

    TRANA   (input) CHARACTER*1   
            Specifies the option op(A):   
            = 'N': op(A) = A    (No transpose)   
            = 'T': op(A) = A**T (Transpose)   
            = 'C': op(A) = A**H (Conjugate transpose = Transpose)   

    TRANB   (input) CHARACTER*1   
            Specifies the option op(B):   
            = 'N': op(B) = B    (No transpose)   
            = 'T': op(B) = B**T (Transpose)   
            = 'C': op(B) = B**H (Conjugate transpose = Transpose)   

    ISGN    (input) INTEGER   
            Specifies the sign in the equation:   
            = +1: solve op(A)*X + X*op(B) = scale*C   
            = -1: solve op(A)*X - X*op(B) = scale*C   

    M       (input) INTEGER   
            The order of the matrix A, and the number of rows in the   
            matrices X and C. M >= 0.   

    N       (input) INTEGER   
            The order of the matrix B, and the number of columns in the   
            matrices X and C. N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,M)   
            The upper quasi-triangular matrix A, in Schur canonical form.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,N)   
            The upper quasi-triangular matrix B, in Schur canonical form.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,N).   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N right hand side matrix C.   
            On exit, C is overwritten by the solution matrix X.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M)   

    SCALE   (output) DOUBLE PRECISION   
            The scale factor, scale, set <= 1 to avoid overflow in X.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1: A and B have common or very close eigenvalues; perturbed   
                 values were used to solve the equation (but the matrices   
                 A and B are unchanged).   

    =====================================================================   


 
   DTRSEN first collects the selected eigenvalues by computing an  
   orthogonal transformation Z to move them to the top left corner of T.  
   In other words, the selected eigenvalues are the eigenvalues of T11  
   in:  
 
                 Z'*T*Z = ( T11 T12 ) n1  
                          (  0  T22 ) n2  
                             n1  n2  
 
   where N = n1+n2 and Z' means the transpose of Z. The first n1 columns  
   of Z span the specified invariant subspace of T.  
 
   If T has been obtained from the real Schur factorization of a matrix  
   A = Q*T*Q', then the reordered real Schur factorization of A is given  
   by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span  
   the corresponding invariant subspace of A.  
 
   The reciprocal condition number of the average of the eigenvalues of  
   T11 may be returned in S. S lies between 0 (very badly conditioned)  
   and 1 (very well conditioned). It is computed as follows. First we  
   compute R so that  
 
                          P = ( I  R ) n1  
                              ( 0  0 ) n2  
                                n1 n2  
 
   is the projector on the invariant subspace associated with T11.  
   R is the solution of the Sylvester equation:  
 
                         T11*R - R*T22 = T12.  
 
   Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote  
   the two-norm of M. Then S is computed as the lower bound  
 
                       (1 + F-norm(R)**2)**(-1/2)  
 
   on the reciprocal of 2-norm(P), the true reciprocal condition number.  
   S cannot underestimate 1 / 2-norm(P) by more than a factor of  
   sqrt(N).  
 
   An approximate error bound for the computed average of the  
   eigenvalues of T11 is  
 
                          EPS * norm(T) / S  
 
   where EPS is the machine precision.  
 
   The reciprocal condition number of the right invariant subspace  
   spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.  
   SEP is defined as the separation of T11 and T22:  
 
                      sep( T11, T22 ) = sigma-min( C )  
 
   where sigma-min(C) is the smallest singular value of the  
   n1*n2-by-n1*n2 matrix  
 
      C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )  
 
   I(m) is an m by m identity matrix, and kprod denotes the Kronecker  
   product. We estimate sigma-min(C) by the reciprocal of an estimate of  
   the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)  
   cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).  
 
   When SEP is small, small changes in T can cause large changes in  
   the invariant subspace. An approximate bound on the maximum angular  
   error in the computed right invariant subspace is  
 
                       EPS * norm(T) / SEP  
 
   =====================================================================  "));
end dtrsyl;

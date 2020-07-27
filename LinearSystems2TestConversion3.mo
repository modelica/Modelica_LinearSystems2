within ;
package LinearSystems2TestConversion3
  package Math
    function complexNumerics
      import Modelica_LinearSystems2.Math.Complex;

    protected
      Complex j = Modelica_LinearSystems2.Math.Complex.j();
      Complex c1=2+3*j;
      Complex c2=3+4*j;
      Complex c3;
      Complex cv[3] "Vector";
      Complex cm[3,2] "Matrix";
    algorithm
      c3 := Complex(2);
      cv := {Complex(2), Complex(1,7), Complex(3,-3)};
      cm := [cv, Modelica_LinearSystems2.Math.Complex.Vectors.reverse(cv)];
      c3 := Modelica_LinearSystems2.Math.Complex.'constructor'(9, -4);
      Complex.'-'.negate(c1);
      Modelica_LinearSystems2.Math.Complex.'-'.subtract(c1, c2);
      Complex.'+'(c1, c2);
      Modelica_LinearSystems2.Math.Complex.'*'(c1, c2);
      Complex.'/'(c1, c2);
      Modelica_LinearSystems2.Math.Complex.'=='(c1, c2);
      Complex.'String'(c3);
      Complex.'abs'(c3);
      Complex.'sqrt'(c3);
      Modelica_LinearSystems2.Math.Complex.'max'(cv);
      Complex.exp(c1);
      Complex.log(c1);
      Complex.sin(c1);
      Complex.cos(c1);
      Complex.arg(c1);
      Complex.conj(c1);
      Modelica_LinearSystems2.Math.Complex.real(c1);
      Modelica_LinearSystems2.Math.Complex.imag(c1);
      Modelica_LinearSystems2.Math.Complex.eigenValues(diagonal({2,3,6}));
      Complex.eigenVectors(diagonal({2,3,6}));
      Modelica_LinearSystems2.Math.Complex.frequency(c1);
      Complex.Vectors.print("c1", c=cv);
      Modelica_LinearSystems2.Math.Complex.Vectors.printHTML(cv);
      Modelica_LinearSystems2.Math.Complex.Vectors.length(cv);
      Complex.Vectors.norm(cv);
      Complex.Vectors.normalize(cv);
      Complex.Vectors.sortComplex(cv);
      Complex.Vectors.multiply(cv,cv);
      Complex.Vectors.reverse(cv);
      Modelica_LinearSystems2.Math.Complex.Matrices.print(cm);
      Complex.Matrices.matMatMul(cm, cm);
      Complex.Matrices.matVecMul(cm, cv);
    end complexNumerics;

    function readMatrixGainTest

    protected
      Real K0[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(m=1,n=1);
      Real K1[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain("filename.mat",m=1,n=1);
      Real K2[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(matrixName="testK",m=1,n=1);
    algorithm
    end readMatrixGainTest;

    function callAllLAPACK "Call all functions from Modelica_LinearSystems2.Math.Matrices.LAPACK"
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgecon
      input Real LU_of_A[:,:] = [1, 2, 3; 3, 4, 5; 3, 2, 3] "LU factroization of a real matrix A";
      input Boolean inf=false "Is true if infinity norm is used and false for 1-norm";
    protected
       Real anorm = 3 "input; norm of A";
    public
      output Real rcond "Reciprocal condition number of A";
      output Integer info "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgees
      input Real A[:,:] = [1, 2, 3; 3, 4, 5; 3, 2, 3] "Square matrix";
      output Real T[:,:] "Real Schur form with A = Z*T*Z'";
      output Real Z[:,:] "orthogonal matrix Z of Schur vectors";
      output Real eval_real[:] "real part of the eigenvectors of A";
      output Real eval_imag[:] "imaginary part of the eigenvectors of A";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev
      output Real alphaReal[:] "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
      output Real alphaImag[:] "Imaginary part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
      output Real lEigenVectors[:,:] "left eigenvectors of matrix A";
      output Real rEigenVectors[:,:] "right eigenvectors of matrix A";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev_eigenValues
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx
      output Real AS[:,:] "AS iss the real Schur form of the balanced version of the input matrix A";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx_eigenValues
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgegv
      input Real B[:,:] = [1, 2, 3; 3, 2, 3; 3, 4, 5] "";
      output Real beta[:] "Denominator of eigenvalue";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgehrd
      input Integer ilo=1 "lowest index where the original matrix had been Hessenbergform";
      input Integer ihi=size(A, 1) "highest index where the original matrix had been Hessenbergform";
      output Real Aout[:,:] "contains the Hessenberg form in the upper triangle and the first subdiagonal and below the first subdiagonal it contains the elementary reflectors which represents (with array tau) as a product the orthogonal matrix Q";
      output Real tau[:] "scalar factors of the elementary reflectors";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqp3
      input Integer lwork1=3*size(A, 2) + 1 "size of work array; should be optimized with Modelica_LinearSystems2.Math.Matrices.Internal.dgeqp3_workdim";
      output Integer jpvt[:] "pivoting indices";
      output Real work[:] "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqrf
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesdd
      output Real sigma[:] "";
      output Real U[:,:] "";
      output Real VT[:,:] "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvd
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvx
      input Boolean transposed=true "True, if matrix A is transformed on input, i.e. system is A**T * X = B";
      output Real X[:,:] "Matrix X[n,nrhs]";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgetrs
      input Real LU[:,:] = LU_of_A "LU factorization of dgetrf of a square matrix";
      input Integer pivots[:] = {3, 4, 5} "Pivot vector of dgetrf";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev
      input Integer nA=size(A, 1) "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev_eigenValues
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dggevx
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dhgeqz
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dhseqr
      input Real H[:,:] = A "";
      input Integer lwork=max(1, size(H, 1)) "";
      input Boolean eigenValuesOnly=true "";
      input String compz="N" "";
      output Real Ho[:,:] "";
      output Real Zo[:,:] "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dlange
      input String norm="1" "specifies the norm, i.e. 1, I, F, M";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dlansy
      input Boolean upper=true "Specifies whether the upper or lower triangular part of A is referenced";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dorghr
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr
    protected
       Real Q[:,:] = A "Orthogonal matrix of elementary reflectors";
    public
      output Real Qout[:,:] "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr_x
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dormhr
      input Real C[:,:] = A "";
      input String side="L" "";
      input String trans="N" "";
      output Real Cout[:,:] "contains the Hessenberg form in the upper triangle and the first subdiagonal and below the first subdiagonal it contains the elementary reflectors which represents (with array tau) as a product the orthogonal matrix Q";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dormqr
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrevc
      input String howmny="B" "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsen
      input String job="N" "";
      input String compq="V" "";
      input Boolean select[:] = {true, false} "";
      output Real To[:,:] "";
      output Real Qo[:,:] "";
      output Real wr[:] "";
      output Real wi[:] "";
      output Integer m "";
      output Real s "";
      output Real sep "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsyl
      input Boolean tranA=false "";
      input Boolean tranB=false "";
      input Integer isgn=1 "";
      output Real scale "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelsx
      output Integer rank "Effective rank of A";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgemm
      input Real a=1 "Factor a";
      input Real b=0 "Factor b";
      input Boolean transA=false "True if transformed A is used";
      input Boolean transB=false "True if transformed B is used";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dpotrf
      output Real Acholesky[:,:] "Cholesky factor";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrmm
      input Real alpha=1 "Factor alpha";
      input Boolean right=true "True if A is right multiplication";
      input Boolean transBout=false "RENAMED HERE, ORIGINALLY: trans; True if op(A) means transposed(A)";
      input Boolean unitTriangular=false "True if A is unit triangular, i.e. all diagonal elements of A are equal to 1";
      output Real Bout[:,:] "Matrix Bout=alpha*op( A )*B,   or   B := alpha*B*op( A )";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsm
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.drot
    protected
       Real x[:] = {5, 4, 2} "input; ";
       Real y[:] = x "input; ";
       Real c = 2 "input; ";
    public
      input Integer incx=1 "";
      input Integer incy=1 "";
      output Real xr[:] "";
      output Real yr[:] "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.drotg
      output Real r "";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsv
      input Real bx[:] = {7, 5, 2} "RENAMED HERE, ORIGINALLY: b; Input vector b";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dposv
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dpocon
      input Real cholA[:,:] = A "Cholesky factor of matrix A";
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelqf
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dorglq
      // Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrtri
      output Real invA[:,:] "Inverse of A";

    algorithm
      anorm :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dlange(A, norm);
      anorm :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dlansy( A, norm, upper);
      (rcond,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgecon( LU_of_A, inf, anorm);
      (T,Z,eval_real,eval_imag,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgees(A);
      (alphaReal,alphaImag,lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev(A);
      (alphaReal,alphaImag,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev_eigenValues(A);
      (alphaReal,alphaImag,lEigenVectors,rEigenVectors,AS,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx(A);
      (alphaReal,alphaImag,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx_eigenValues(A);
      (alphaReal,alphaImag,beta,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgegv(A, B);
      (Aout,tau,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgehrd( A, ilo, ihi);
      (Aout,jpvt,tau,info,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqp3(A, lwork1);
      (Aout,tau,info,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqrf(A, lwork1);
      (sigma,U,VT,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesdd(A);
      (sigma,U,VT,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvd(A);
      (X,info,rcond) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvx( A, B, transposed);
      X :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgetrs( LU, pivots, B);
      (alphaReal,alphaImag,beta,lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev( A, B, nA);
      (alphaReal,alphaImag,beta,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev_eigenValues(A, B);
      (alphaReal,alphaImag,beta,lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggevx(A, B);
      (alphaReal,alphaImag,beta,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dhgeqz(A, B);
      (alphaReal,alphaImag,info,Ho,Zo,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dhseqr( H, lwork, eigenValuesOnly, compz, Z);
      (Aout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorghr( A, ilo, ihi, tau);
      (Q,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorglq(A, tau);
      (Qout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr(Q, tau);
      (Aout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr_x(Q, tau);
      (Cout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dormhr( C, A, tau, side, trans, ilo, ihi);
      (Cout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dormqr( C, A, tau, side, trans);
      (lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrevc( T, side, howmny, Q);
      (To,Qo,wr,wi,m,s,sep,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsen( job, compq, select, T, Q);
      (X,scale,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsyl( A, B, C, tranA, tranB, isgn);
      (X,info,rank) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelsx( A, B, rcond);
      Cout :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgemm( A, B, C, a, b, transA, transB);
      (Acholesky,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dpotrf(A, upper);
      Bout :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrmm( A, B, alpha, right, upper, transBout, unitTriangular);
      X :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsm( A, B, alpha, right, upper, transBout, unitTriangular);
      (c,s,r) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.drotg(a, b);
      x :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrsv( A, bx, upper, transBout, unitTriangular);
      (xr,yr) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.drot( x, y, c, s, incx, incy);
      (X,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dposv( A, B, upper);
      (rcond,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dpocon( cholA, anorm, upper);
      (Aout,tau,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelqf(A);
      (invA,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrtri(A, upper);
    end callAllLAPACK;
  end Math;

  package Types

    model Issue13 "Conversion test concerning github issue #13"
      import Modelica_LinearSystems2.Types;
      parameter Modelica_LinearSystems2.Types.Grid grid = Modelica_LinearSystems2.Types.Grid.Equidistant;
      parameter Types.AnalogFilter analogFilter = Modelica_LinearSystems2.Types.AnalogFilter.Chebyshev;
      parameter Modelica_LinearSystems2.Types.FilterType filterType = Types.FilterType.HighPass;
      parameter Modelica_LinearSystems2.Types.Method method = Modelica_LinearSystems2.Types.Method.Trapezoidal;
      parameter Modelica_LinearSystems2.Types.StaircaseMethod staircaseMethod = Modelica_LinearSystems2.Types.StaircaseMethod.QR;
      parameter Types.TimeResponse timeResponse = Types.TimeResponse.Impulse;
      parameter Modelica_LinearSystems2.Types.Window window = Modelica_LinearSystems2.Types.Window.Bartlett;

    end Issue13;
  end Types;

  package Streams

    function readSystemDimensionTest

    protected
      Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension();
      Integer xuy2[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension2();
      Integer xuy3[3] = Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension();
    algorithm

    end readSystemDimensionTest;
  end Streams;
  annotation (uses(Modelica_LinearSystems2(version="2.4.0")));
end LinearSystems2TestConversion3;

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
      Modelica_LinearSystems2.Math.Complex.Internal.C_transpose(cm);
    end complexNumerics;

    function readMatrixGainTest

    protected
      Real K0[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(m=1,n=1);
      Real K1[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain("filename.mat",m=1,n=1);
      Real K2[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(matrixName="testK",m=1,n=1);
    algorithm
    end readMatrixGainTest;

    function printVectorTest
      output Boolean identical "= true, if strings are identical";
    protected
      String s2 = Modelica_LinearSystems2.Math.Vectors.printVector( {3,33,7}, 2, "vec");
      String s3 = Modelica_LinearSystems2.Math.Vectors.printVector( v={3,33,7}, significantDigits=2);
      String sm = Modelica.Math.Vectors.toString({3, 33, 7}, "vec", 2);
    algorithm
      identical := Modelica.Utilities.Strings.isEqual(s2, sm);

    end printVectorTest;

    function printMatrixTest
      output Boolean identical "= true, if strings are identical";
    protected
      Real A[4,3] = [1, 0,  0; 6, 5,  0; 1, -2,  2; 0, 4,  44];

      String m2 = Modelica_LinearSystems2.Math.Matrices.printMatrix(A, 2, "mtx");
      String m3 = Modelica_LinearSystems2.Math.Matrices.printMatrix(M=A, significantDigits=3);
      String mm = Modelica.Math.Matrices.toString(A, "mtx", 2);
    algorithm
      identical := Modelica.Utilities.Strings.isEqual(m2, mm);

    end printMatrixTest;

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
      //
      (Aout,tau,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgehrd( A, ilo, ihi);
      (Aout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorghr( A, ilo, ihi, tau);
      (Cout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dormhr( C, A, tau, side, trans, ilo, ihi);
      //
      (Aout,tau,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelqf(A);
      (Q,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorglq(Aout, tau);
      //
      (Aout,jpvt,tau,info,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqp3(A, lwork1);
      //
      (Aout,tau,info,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeqrf(A, lwork1);
      (Qout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr(Q, tau);
      //
      (sigma,U,VT,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesdd(A);
      (sigma,U,VT,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvd(A);
      (X,info,rcond) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgesvx( A, B, transposed);

      //(LU, pivots, info) := Modelica.Math.Matrices.LAPACK.dgetrf(A);
      X :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dgetrs(LU, pivots, B);
      (alphaReal,alphaImag,beta,lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev( A, B, nA);
      (alphaReal,alphaImag,beta,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev_eigenValues(A, B);
      (alphaReal,alphaImag,beta,lEigenVectors,rEigenVectors,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dggevx(A, B);
      (alphaReal,alphaImag,beta,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dhgeqz(A, B);
      (alphaReal,alphaImag,info,Ho,Zo,work) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dhseqr( H, lwork, eigenValuesOnly, compz, Z);
      (Aout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr_x(Q, tau);

      (Cout,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dormqr( C, A, tau, side, trans);
      (T, Z, alphaReal, alphaImag) :=Modelica_LinearSystems2.Math.Matrices.rsf(A);
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

      (invA,info) :=Modelica_LinearSystems2.Math.Matrices.LAPACK.dtrtri(A, upper);
    end callAllLAPACK;

    function callVectors
      output Integer index;
      output Real length;

    protected
      Integer v[3] = {1, 2, 3};
      Integer e1 = 2;

    algorithm
      index := Modelica_LinearSystems2.Math.Vectors.find(e1,v);
      length := Modelica_LinearSystems2.Math.Vectors.length(v);

    end callVectors;

    function callMatrices
      output Real r;
      output Real det;
      output Real n "p-norm of matrix A";
      output Real x[size(b1, 1)] "Vector x such that A*x = b1";
      output Real xr[size(A, 2)]
        "Vector x such that min|A*x-b|^2 if size(A,1) >= size(A,2) or min|x|^2 and A*x=b, if size(A,1) < size(A,2)";
      output Real xls[size(Als, 1)] "State vector";
      output Real t "Trace of A";
      output Real HC[size(A, 1),size(A, 2)];
      output Real HH[size(A, 1),size(A, 2)] "Upper Hessenberg form";
      output Real H[size(A, 1),size(A, 2)];
      output Real U[size(A, 1),size(A, 2)];
      output Real Aflip[size(A, 1), size(A, 2)];
      output Real Afud[size(A, 1), size(A, 2)];
      output Real X2[size(A, 2), size(B, 2)];
      output Real LU[size(A, 1),size(A, 2)];
      output Integer pivots[3]; //min(size(A, 1), size(A, 2))];
      output Real x1[size(A, 1)];
      output Real X1[size(B, 1), size(B, 2)];
      output Real X[size(B, 1),size(B, 2)] "Matrix X such that A*X = B";
      output Real Z[size(A, 2), :] "Orthonormal nullspace of matrix A";
      output Real S[size(A, 1), size(A, 2)] "Real Schur form of A";
      output Real rcond "Reciprocal condition number of A";

    protected
      Real A[3,3] = [1, 0,  0; 6, 5,  0; 1, -2,  2] "Square matrix";
      Real SA[size(A, 1),size(A, 2)] = A*transpose(A);
      Real B[size(A, 1),1] = [1;0;1];
      Real b1[size(A, 1)] = {7,13,10};

      Real QZ[size(A, 1),size(A, 2)];
      Real alphaReal[size(A, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
      Real alphaImag[size(A, 1)] "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

      Real Als[:,:] = [-1.0, 0.0, 0.0; 0.0, -2.0, 0.0; 0.0, 0.0, -3.0];
      Real Bls[:,:] = [1.0; 1.0; 0.0];
      Real Cls[:,:] = [1.0, 1.0, 1.0];
      Real Dls[:,:] = [0.0];
      Integer nx = size(Als,2);
      Integer nu = size(Bls,2);
      Integer ny = size(Cls,2);

    algorithm
      HC := Modelica_LinearSystems2.Math.Matrices.cholesky(SA, true);
      r := Modelica_LinearSystems2.Math.Matrices.conditionNumber(A);
      det := Modelica_LinearSystems2.Math.Matrices.det(A);
      Aflip := Modelica_LinearSystems2.Math.Matrices.fliplr(A);
      Afud := Modelica_LinearSystems2.Math.Matrices.flipud(A);
      (H, U) := Modelica_LinearSystems2.Math.Matrices.hessenberg(A);
      xr := Modelica_LinearSystems2.Math.Matrices.leastSquares(A,b1);
      X2 := Modelica_LinearSystems2.Math.Matrices.leastSquares2(A,B);
      xls := Modelica_LinearSystems2.Math.Matrices.equalityLeastSquares(
        Als, -Bls*fill(1,nu), Cls, Dls*fill(1,nu));
      (LU, pivots) := Modelica_LinearSystems2.Math.Matrices.LU(A);
      x1 := Modelica_LinearSystems2.Math.Matrices.LU_solve(LU, pivots, b1);
      X1 := Modelica_LinearSystems2.Math.Matrices.LU_solve2(LU, pivots, B);
      HH := Modelica_LinearSystems2.Math.Matrices.toUpperHessenberg(A);
      n := Modelica_LinearSystems2.Math.Matrices.norm(A);
      Z := Modelica_LinearSystems2.Math.Matrices.nullspace(A);
      (S, QZ, alphaReal, alphaImag) := Modelica_LinearSystems2.Math.Matrices.rsf2(A);
      rcond := Modelica_LinearSystems2.Math.Matrices.rcond(A);
      x := Modelica_LinearSystems2.Math.Matrices.solve(A, b1);
      X := Modelica_LinearSystems2.Math.Matrices.solve2(A, B);
      t := Modelica_LinearSystems2.Math.Matrices.trace(A);

    end callMatrices;

    package Polynomials

      function polynomialDegree
        output Integer degree;
        output Integer degree2;
      protected
        Modelica_LinearSystems2.Math.Polynomial p = Modelica_LinearSystems2.Math.Polynomial({0,0,4,0,1});
      algorithm
        degree := Modelica_LinearSystems2.Math.Polynomial.degree(p);
        // shall be = 2
        degree2 := Modelica_LinearSystems2.Math.Polynomial.degree2(p);
        // shall be = 2

      end polynomialDegree;
    end Polynomials;
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

      output Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension(
        Modelica_LinearSystems2.DataDir + "abcd_siso.mat");
      output Integer xuy2[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension2(
        Modelica_LinearSystems2.DataDir + "abcd_siso.mat");
      output Integer xuy3[3] = Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension(
        Modelica_LinearSystems2.DataDir + "abcd_siso.mat");
    algorithm

    end readSystemDimensionTest;

    function readMatrixXTest
      import Modelica_LinearSystems2.Internal.Streams;

      input String fileName = Modelica_LinearSystems2.DataDir + "abcd.mat";
      input String mName = "ABCD";

      output Integer dim[2] = Modelica_LinearSystems2.Internal.Streams.readMatrixOnFileSize(fileName, mName);
      output Real A[:,:] = Streams.ReadMatrixA(fileName, mName);
      output Real B[:,:] = Streams.ReadMatrixB(fileName, mName);
      output Real C[:,:] = Modelica_LinearSystems2.Internal.Streams.ReadMatrixC(fileName, mName);
      output Real D[:,:] = Modelica_LinearSystems2.Internal.Streams.ReadMatrixD(fileName, mName);
    protected
      Integer xuy[3] = Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension(
        fileName, mName);
      Integer nx = xuy[1];
      Integer nu = xuy[2];
      Integer ny = xuy[3];

      Real A2[nx,nx] = Streams.ReadMatrixA2(fileName, mName, nx=nx);
      Real B2[nx,nu] = Streams.ReadMatrixB2(fileName, mName, nx=nx, nu=nu);
      Real C2[ny,nx] = Modelica_LinearSystems2.Internal.Streams.ReadMatrixC2(
        fileName, mName, nx=nx, ny=ny);
      Real D2[ny,nu] = Modelica_LinearSystems2.Internal.Streams.ReadMatrixD2(
        fileName, mName, nx=nx, nu=nu, ny=ny);
    algorithm

      annotation();
    end readMatrixXTest;

    function otherClassesTest

    protected
      Modelica_LinearSystems2.Internal.Streams.AnalyseOptions ao;
      Modelica_LinearSystems2.StateSpace ssi = Modelica_LinearSystems2.StateSpace(
        A=[1,1;3,0],
        B=[1;1],
        C=[1,0],
        D=[0],
        xNames={"x1","x2"},
        uNames={"u1"}, yNames={"y1"});

      String s = Modelica_LinearSystems2.Internal.Streams.stateSpaceString_html(ssi);
    algorithm

      annotation();
    end otherClassesTest;
  end Streams;

  package Controllers
    extends Modelica.Icons.ExamplesPackage;

    model FirstExample
      extends Modelica.Icons.Example;

      import Modelica_LinearSystems2;

      parameter Modelica.Units.SI.AngularFrequency w=10
        "Undamped natural frequency";
      parameter Real D=0.1 "Damping ratio";

      Modelica.Blocks.Sources.Step step(
        startTime=0.5,
        height=1.2,
        offset=0.2) annotation (
          Placement(transformation(extent={{-80,0},{-60,20}})));
      Modelica_LinearSystems2.Controller.StateSpace stateSpace(
        x_start={0.1,0},
        initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.InitialState,
        system=Modelica_LinearSystems2.StateSpace(
          A=[0,1; -w*w,-2*w*D],
          B=[0; w*w],
          C=[1,0],
          D=[0]),
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
        annotation(Placement(transformation(extent={{-20,40},{0,60}})));
      Modelica_LinearSystems2.Controller.TransferFunction transferFunction(system(n={1,2}, d={1,2,3}), blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption) annotation (Placement(transformation(extent={{-20,0},{0,20}})));
      Modelica_LinearSystems2.Controller.ZerosAndPoles zerosAndPoles(system(
          n1={1},
          n2=fill(0, 0, 2),
          d1=fill(0, 0),
          d2=[1,1; 1,1]), blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption) annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
      inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(sampleTime=0.1, blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Continuous) annotation (Placement(transformation(extent={{60,60},{80,80}})));
    equation
      connect(step.y, stateSpace.u[1]) annotation (Line(
          points={{-59,10},{-40,10},{-40,50},{-22,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(transferFunction.u, step.y) annotation (Line(
          points={{-22,10},{-59,10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(zerosAndPoles.u, step.y) annotation (Line(
          points={{-22,-30},{-40,-30},{-40,10},{-59,10}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (    experiment(StopTime=5));
    end FirstExample;

    model Discretization1
      extends Modelica.Icons.Example;

      import Modelica_LinearSystems2.Controller;

      parameter Real w=20 "Angular frequency";
      parameter Real D=0.1 "Damping";

      Controller.SecondOrder continuous(w=w, D=D)
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(sampleTime=
            0.01)
        annotation (Placement(transformation(extent={{60,60},{80,80}})));
      Modelica.Blocks.Sources.Step step(
        height=1.2,
        offset=0.2,
        startTime=0.1) annotation (
          Placement(transformation(extent={{-80,40},{-60,60}})));
      Modelica_LinearSystems2.Controller.SecondOrder explicitEuler(
        w=w,
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ExplicitEuler)
        annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

      Modelica_LinearSystems2.Controller.SecondOrder implicitEuler(
        w=w,
        D=D,
        blockType=Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ImplicitEuler)
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));

      Modelica_LinearSystems2.Controller.SecondOrder trapezoid(
        w=w,
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.Trapezoidal)
        annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));

      Controller.SecondOrder impulseExact(
        w=w,
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Controller.Types.MethodWithGlobalDefault.ImpulseExact)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));

      Modelica_LinearSystems2.Controller.SecondOrder stepExact(
        w=w,
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.StepExact)
        annotation (Placement(transformation(extent={{20,-20},{40,0}})));

      Modelica_LinearSystems2.Controller.SecondOrder rampExact(
        w=w,
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.RampExact)
        annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

    equation
      connect(step.y, continuous.u) annotation (Line(
          points={{-59,50},{-42,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, explicitEuler.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,10},{-42,10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, implicitEuler.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,-30},{-42,-30}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, trapezoid.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,-70},{-42,-70}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, impulseExact.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,80},{0,80},{0,30},{18,30}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, stepExact.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,80},{0,80},{0,-10},{18,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step.y, rampExact.u) annotation (Line(
          points={{-59,50},{-52,50},{-52,80},{0,80},{0,-50},{18,-50}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        experiment(Tolerance=1e-006),
        Documentation(info=""));
    end Discretization1;

    model Discretization2
      extends Modelica.Icons.Example;

      import Modelica_LinearSystems2.Controller;
      import Modelica_LinearSystems2.Controller.SecondOrder;

      parameter Real w=20 "Angular frequency";
      parameter Real D=0.1 "Damping";

      inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(
        sampleTime=0.01)
        annotation (Placement(transformation(extent={{60,60},{80,80}})));

      SecondOrder impulseExact(
        D=D,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ImpulseExact,
        w=w)
        annotation (Placement(transformation(extent={{0,-30},{20,-10}})));

      SecondOrder continuous(D=D, w=w)
        annotation (Placement(transformation(extent={{0,10},{20,30}})));
      Controller.Derivative derivative(T=1e-8) annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
      Modelica.Blocks.Sources.Pulse pulse(
        startTime=0.1,
        period=1,
        width=sampleClock.sampleTime*100)
        annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
      Modelica.Blocks.Sources.Step step1(
        startTime=0.1,
        height=1,
        offset=0) annotation (Placement(transformation(extent={{-80,10},{-60,30}})));

    equation
      connect(pulse.y, impulseExact.u) annotation (Line(
          points={{-59,-20},{-2,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(continuous.u, derivative.y) annotation (Line(
          points={{-2,20},{-19,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step1.y, derivative.u) annotation (Line(
          points={{-59,20},{-42,20}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        experiment(Tolerance=1e-006),
        Documentation(info=""));
    end Discretization2;

    model InverseDoublePendulumWithObserver
      extends Modelica.Icons.Example;

      extends Modelica_LinearSystems2.Controller.Templates.SimpleObserverStateSpaceControl(
        redeclare Modelica_LinearSystems2.Controller.Examples.Components.InverseDoublePendulum3 plant(
          additionalMeasurableOutputs=true,
          m_trolley=1,
          n=6,
          phi2_start=0,
          length=1,
          cartDisturbance=true,
          bodyDisturbance=true,
          l=2,
          secondAngle=false,
          m_load=1,
          phi1_start=1.5707963267949),
        preFilter(
          matrixName="M_pa",
          fileName=Modelica_LinearSystems2.Controller.DataDir + "inverseDoublePendulumController.mat",
          matrixOnFile=true),
        feedbackMatrix(
          matrixOnFile=true,
          matrixName="K_pa",
          fileName=Modelica_LinearSystems2.Controller.DataDir + "inverseDoublePendulumController.mat"),
        sampleClock(sampleTime=0.002, blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Continuous),
        observer(
          systemName="stateSpace",
          matrixOnFile=true,
          initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.InitialState,
          methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.StepExact,
          x_start={0,0,0,0,0,0},
          observerMatrixName="K_ob2",
          blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption,
          withDelay=true,
          fileName=Modelica_LinearSystems2.Controller.DataDir + "inverseDoublePendulumControllerO.mat"));

      Modelica.Blocks.Sources.Pulse pulse(
        offset=0,
        startTime=1,
        width=50,
        period=30,
        amplitude=5)
        annotation (Placement(transformation(extent={{-140,-10},{-120,10}})));
      Modelica_LinearSystems2.Controller.Examples.Components.AccelerationLimiter accelerationLimiter(
        v_limit=20,
        velocityLimitation=false,
        withDelay2=false,
        a_limit=1) annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      Modelica_LinearSystems2.Controller.Noise noise(
        firstSeed={43,123,162},
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        y_min=-0.005,
        y_max=0.005,
        sampleFactor=200) annotation (Placement(transformation(extent={{130,40},{110,60}})));
      Modelica_LinearSystems2.Controller.Noise noise1(
        sampleFactor=100,
        blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
        y_min=-0.025,
        y_max=0.025) annotation (Placement(transformation(extent={{50,40},{70,60}})));

    initial equation
      //feedback.y = {0.0};
      // plant.u = {0.0};

    equation
      connect(pulse.y, accelerationLimiter.u) annotation (Line(
          points={{-119,0},{-112,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(noise1.y, plant.dist) annotation (Line(
          points={{71,50},{85.2,50},{85.2,5.6}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(noise.y, plant.dist2) annotation (Line(
          points={{109,50},{95.4,50},{95.4,5.6}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(preFilter.u[1], accelerationLimiter.s) annotation (Line(
          points={{-72,0},{-80,0},{-80,6},{-89,6}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-140,-100},{140,
                100}}), graphics={Text(
              extent={{44,76},{134,64}},
              lineColor={0,0,0},
              textString="disturbance"), Rectangle(extent={{-76,28},{108,-60}},
                lineColor={255,0,0})}),
        experiment(
          StopTime=60,
          __Dymola_NumberOfIntervals=2000,
          Tolerance=1e-005),
        Documentation(info=""));
    end InverseDoublePendulumWithObserver;

    model MixingUnit
      extends Modelica.Icons.Example;

      extends Modelica_LinearSystems2.Controller.Templates.TwoDOFinverseModelController(
        redeclare Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit plant_inv(mixingUnit(
            c(start=c_start, fixed=true),
            T_c(start=T_c_start, fixed=true),
            T(start=T_start, fixed=true))),
        redeclare Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit plant(mixingUnit(c(start=c_start, fixed=true), T(start=T_start, fixed=true))),
        filter(
          order=3,
          normalized=false,
          f_cut=freq,
          initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.NoInit),
        redeclare Modelica_LinearSystems2.Controller.PI controller(
          k=10,
          T=10,
          initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.InitialState));

      import Modelica.Units.SI;
      parameter Real x10 = 0.42
        "Initial value of state x1 (related concentration of substance A in tank)";
      parameter Real x10_inv = 0.6 "Initial value of state x1 of inverted model";
      parameter Real x20 = 0.01
        "Initial value of state x2 (related temperature in tank)";
      parameter Real u0 = -0.0224
        "Initial related temperature of cooling medium [-]";
      parameter SI.Frequency freq = 1/300 "Critical frequency of filter";

      final parameter Real c0 = 0.848
        "Nominal concentration of substance A on intake";
      final parameter SI.Temperature T0 = 308.5
        "Nominal temperature of substance A on intake";
      final parameter Real c_start(unit="mol/l") = c0*(1-x10)
        "Initial concentration of substance A in tank";
      final parameter Real c_inv_start(unit="mol/l") = c0*(1-x10_inv)
        "Initial concentration of substance A in tank";
      final parameter SI.Temperature T_start = T0*(1+x20)
        "Initial temperature in tank";
      final parameter Real c_high_start(unit="mol/l") = c0*(1-0.72)
        "Concentration change height";
      final parameter SI.Temperature T_c_start = T0*(1+u0)
        "Initial temperature of cooling medium";

      Modelica.Blocks.Sources.Step step1(
        height=c_high_start - c_start,
        offset=c_start,
        startTime=25)
        annotation (Placement(transformation(extent={{-120,10},{-100,30}}, rotation=0)));
      inner Modelica_LinearSystems2.Controller.SampleClock sampleClock annotation (Placement(transformation(extent={{80,80},{100,100}})));
    equation
      connect(step1.y, filter.u) annotation (Line(
          points={{-99,20},{-92,20}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=500),
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-100},{120,
                100}}), graphics));
    end MixingUnit;

    model ExamplesUtilities
      extends Modelica.Icons.Example;

      import Modelica_LinearSystems2.Controller;
      import Modelica_LinearSystems2.Controller.Examples.Components;
      import Modelica_LinearSystems2.Controller.Examples.Components.DoublePendulum2;

      Modelica_LinearSystems2.Controller.Examples.Components.AccelerationLimiter accelerationLimiter annotation (Placement(transformation(extent={{-50,80},{-30,100}})));
      Controller.Examples.Components.DoublePendulum doublePendulum annotation (Placement(transformation(extent={{59,78},{89,98}})));
      DoublePendulum2 doublePendulum2 annotation (Placement(transformation(extent={{60,0},{80,20}})));
      Components.InverseDoublePendulum inverseDoublePendulum annotation (Placement(transformation(extent={{60,40},{80,60}})));
      Modelica_LinearSystems2.Controller.Examples.Components.InverseDoublePendulum2 inverseDoublePendulum2 annotation (Placement(transformation(extent={{60,-40},{80,-20}})));
      Controller.Examples.Components.InverseDoublePendulum3 inverseDoublePendulum3 annotation (Placement(transformation(extent={{60,-80},{80,-60}})));
      Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit mixingUnit annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
      Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit1 mixingUnit1 annotation (Placement(transformation(extent={{-50,-42},{-30,-22}})));
      Components.SeriesConnection seriesConnection annotation (Placement(transformation(extent={{-50,-80},{-30,-60}})));
      Modelica_LinearSystems2.Controller.Examples.Components.TwoPoint twoPoint annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
      Modelica.Blocks.Sources.Constant const(k=0.1) annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
      Modelica.Blocks.Sources.Constant const1(k=0.1) annotation (Placement(transformation(extent={{0,60},{20,80}})));
      Modelica.Blocks.Sources.Constant const2(k=0.1) annotation (Placement(transformation(extent={{0,-20},{20,0}})));
    equation
      connect(const.y, accelerationLimiter.u) annotation (Line(points={{-79,30},{-70,30},{-70,90},{-52,90}}, color={0,0,127}));
      connect(const.y, twoPoint.u) annotation (Line(points={{-79,30},{-70,30},{-70,50},{-52,50}}, color={0,0,127}));
      connect(const.y, mixingUnit.u) annotation (Line(points={{-79,30},{-70,30},{-70,10},{-52,10}}, color={0,0,127}));
      connect(const.y, mixingUnit1.T_c) annotation (Line(points={{-79,30},{-70,30},{-70,-32},{-52,-32}}, color={0,0,127}));
      connect(const1.y, doublePendulum.u) annotation (Line(points={{21,70},{40.5,70},{40.5,88},{57,88}}, color={0,0,127}));
      connect(const1.y, inverseDoublePendulum.u) annotation (Line(points={{21,70},{40,70},{40,50},{58,50}}, color={0,0,127}));
      connect(const2.y, doublePendulum2.u[1]) annotation (Line(points={{21,-10},{40,-10},{40,10},{58,10}}, color={0,0,127}));
      connect(const2.y, inverseDoublePendulum2.u[1]) annotation (Line(points={{21,-10},{40,-10},{40,-30},{58,-30}}, color={0,0,127}));
      connect(const2.y, inverseDoublePendulum3.u[1]) annotation (Line(points={{21,-10},{40,-10},{40,-70},{58,-70}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end ExamplesUtilities;

    package Templates
      model DoublePendulum "Crane trolley controlled by a state feedback controller"
        extends Modelica.Icons.Example;
        extends Modelica_LinearSystems2.Controller.Templates.SimpleStateSpaceControl
                                                 (
          redeclare Modelica_LinearSystems2.Controller.Examples.Components.DoublePendulum2 plant(
            additionalMeasurableOutputs=true,
            m_trolley=5,
            m_load=20,
            length=2,
            n=6,
            l=6,
            phi1_start=-0.69813170079773,
            phi2_start=-0.34906585039887),
          preFilter(
            matrixName="M_pa",
            fileName=Modelica_LinearSystems2.DataDir + "doublePendulumController.mat",
            matrixOnFile=true),
          feedbackMatrix(
            matrixOnFile=true,
            matrixName="K_pa",
            fileName=Modelica_LinearSystems2.DataDir + "doublePendulumController.mat"),
          sampleClock(sampleTime=0.01, blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Continuous));

        Modelica.Blocks.Sources.Pulse pulse(
          offset=0,
          amplitude=3,
          width=50,
          startTime=5,
          period=10)
          annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
        Modelica_LinearSystems2.Controller.FirstOrder firstOrder(T=0.25) annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      equation
        connect(firstOrder.u, pulse.y) annotation (Line(
            points={{-72,0},{-76,0},{-76,0},{-79,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(firstOrder.y, preFilter.u[1]) annotation (Line(
            points={{-49,0},{-46,0},{-46,0},{-42,0}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (
          experiment(
            StopTime=40,
            __Dymola_NumberOfIntervals=2000,
            Tolerance=1e-005));
      end DoublePendulum;

      model PartialPlantMIMOExtend
        extends Modelica_LinearSystems2.Controller.Templates.PartialPlantMIMO(
          l=1,
          additionalMeasurableOutputs=true,
          m=1,
          n=1);
        Modelica.Blocks.Math.Gain gain(k=1) annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      equation
        connect(u[1], gain.u) annotation (Line(points={{-120,0},{-62,0}}, color={0,0,127}));
        connect(gain.y, y[1]) annotation (Line(points={{-39,0},{110,0}}, color={0,0,127}));
        connect(gain.y, ym[1]) annotation (Line(points={{-39,0},{0,0},{0,-110}}, color={0,0,127}));
      end PartialPlantMIMOExtend;

      model PartialPlantSISOExtend
        extends Modelica_LinearSystems2.Controller.Templates.PartialPlantSISO;
        Modelica.Blocks.Math.Gain gain(k=1) annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      equation
        connect(u, gain.u) annotation (Line(points={{-120,0},{-92,0},{-92,0},{-62,0}}, color={0,0,127}));
        connect(gain.y, y) annotation (Line(points={{-39,0},{34,0},{34,0},{110,0}}, color={0,0,127}));
        connect(gain.y, ym) annotation (Line(points={{-39,0},{0,0},{0,-110}}, color={0,0,127}));
      end PartialPlantSISOExtend;

      model PlantTemplate_SISOExtend
        extends Modelica_LinearSystems2.Controller.Templates.PlantTemplate_SISO(l=2, additionalMeasurableOutputs=true);
        Modelica.Blocks.Math.Gain gain(k=1) annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Modelica.Blocks.Sources.Constant const(k=0) annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
      equation
        connect(u, gain.u) annotation (Line(points={{-120,0},{-62,0}}, color={0,0,127}));
        connect(gain.y, y) annotation (Line(points={{-39,0},{110,0}}, color={0,0,127}));
        connect(gain.y, ym[1]) annotation (Line(points={{-39,0},{0,0},{0,-112.5}},
                                                                                 color={0,0,127}));
        connect(const.y, ym[2]) annotation (Line(points={{-39,-50},{-2,-50},{-2,-102},{0,-102},{0,-107.5}}, color={0,0,127}));
      end PlantTemplate_SISOExtend;

      model PlantTemplate_SISOIntantiate
        replaceable Modelica_LinearSystems2.Controller.Templates.PlantTemplate_SISO plantTemplate_SISO(additionalMeasurableOutputs=true, l=2) annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
        replaceable Modelica_LinearSystems2.Controller.Templates.PartialPlantSISO plant annotation (Placement(transformation(extent={{-20,20},{0,40}})));
        Modelica.Blocks.Sources.Constant const(k=0) annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Modelica.Blocks.Math.Gain gain(k=1) annotation (Placement(transformation(extent={{40,0},{60,20}})));
        Modelica.Blocks.Routing.Multiplex2 multiplex2 annotation (Placement(transformation(extent={{40,-70},{60,-50}})));
      equation
        connect(const.y, plantTemplate_SISO.u) annotation (Line(points={{-59,0},{-40,0},{-40,-30},{-22,-30}}, color={0,0,127}));
        connect(const.y, plant.u) annotation (Line(points={{-59,0},{-40,0},{-40,30},{-22,30}}, color={0,0,127}));
        connect(plant.ym, gain.u) annotation (Line(points={{-10,19},{-10,10},{38,10}}, color={0,0,127}));
        connect(plantTemplate_SISO.ym[1], multiplex2.u1[1]) annotation (Line(points={{-10,-41.25},{-10,-54},{38,-54}}, color={0,0,127}));
        connect(plantTemplate_SISO.ym[2], multiplex2.u2[1]) annotation (Line(points={{-10,-40.75},{-10,-66},{38,-66}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end PlantTemplate_SISOIntantiate;
    end Templates;

    model DoublePendulum
      Modelica_LinearSystems2.Controller.Examples.Components.DoublePendulum doublePendulum annotation (Placement(transformation(extent={{1,20},{31,40}})));
      Modelica_LinearSystems2.Controller.Examples.Components.InverseDoublePendulum inverseDoublePendulum annotation (Placement(transformation(extent={{0,-40},{20,-20}})));
      Modelica.Blocks.Sources.Constant const(k=0.1) annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
    equation
      connect(const.y, doublePendulum.u) annotation (Line(points={{-39,0},{-20,0},{-20,30},{-1,30}}, color={0,0,127}));
      connect(const.y, inverseDoublePendulum.u) annotation (Line(points={{-39,0},{-20,0},{-20,-30},{-2,-30}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end DoublePendulum;
  end Controllers;
  annotation (uses(Modelica_LinearSystems2(version="2.4.0"), Modelica(version="4.0.0")));
end LinearSystems2TestConversion3;

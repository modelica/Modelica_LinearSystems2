within Modelica_LinearSystems2.Math.Matrices;
package Examples
  function exampleHessenberg
    import Modelica_LinearSystems2.Math.Matrices;
    input String fileName=DataDir + "m.mat"
                                annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="matrix file")));
    input String matrixName="A" "Name of the matrix";
  protected
    Real M[:,:]=Matrices.fromFile(fileName, matrixName);
    Real tau[size(M, 1) - 1];
    Real Q[size(M, 1),size(M, 2)];

    Real V[size(M, 1),size(M, 2)];
    Real H[size(M, 1),size(M, 2)];
    Integer info;
  algorithm
    H := Modelica_LinearSystems2.Math.Matrices.Internal.hessenberg2(M);
    Matrices.printMatrix(
        M,
        6,
        "M");
    Math.Matrices.printMatrix(
        H,
        6,
        "Hu");
    H := Modelica_LinearSystems2.Math.Matrices.Internal.hessenberg2(M, "l");
    Math.Matrices.printMatrix(
        H,
        6,
        "Hl");
    (H,V,tau) := Modelica_LinearSystems2.Math.Matrices.toUpperHessenberg(
        M,
        1,
        size(M, 1));
    Math.Matrices.printMatrix(
        H,
        6,
        "H_dgehrd");
    Math.Matrices.printMatrix(
        V,
        6,
        "V_dgehrd");
    Math.Vectors.printVector(
        tau,
        6,
        "tau");
    Q := Modelica_LinearSystems2.Math.Matrices.orthogonalQ(
        V,
        tau,
        1,
        size(V, 1));
    Math.Matrices.printMatrix(
        Q*H*transpose(Q),
        6,
        "Q*H*Q'");
    H := Modelica_LinearSystems2.Math.Matrices.hessenberg(M);
    Math.Matrices.printMatrix(
        H,
        6,
        "H_hess_lapack");

  end exampleHessenberg;

  function exampleQR
    "Example for the usage of QR2-function, QR factorization with colomns pivoting"
    import Modelica_LinearSystems2.Math.Matrices;
    input String fileName=DataDir + "m.mat"
                                annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="matrix file")));
    input String matrixName="A" "Name of the matrix";
  protected
    Real M[:,:]=Matrices.fromFile(fileName, matrixName);
    Real tau[min(size(M, 1), size(M, 2))];
    Integer p[min(size(M, 1), size(M, 2))];
    Real Q[size(M, 1),size(M, 2)];
    Real R[size(M, 2),size(M, 2)];
    Real P[size(M, 2),size(M, 2)]=fill(
          0,
          size(M, 2),
          size(M, 2));
    Real M2[size(M, 1),size(M, 2)];
    Real QR[size(M, 1),size(M, 2)];
    Real QR2[size(M, 1),size(M, 2)];

    Integer info;
  algorithm

    (Q,R,p,tau) := Modelica_LinearSystems2.Math.Matrices.Internal.QR2(
                                M);
    QR := Q*R;
    for i in 1:size(M, 2) loop
      P[p[i], i] := 1;
      M2[:, p[i]] := M[:, i];
      QR2[:, i] := QR[:, p[i]];

    end for;
    Modelica.Utilities.Streams.print(
      "Show results of QR2 - QR factorization with pivoting:\n-----------------------------------------------------");
    Matrices.printMatrix(
        M,
        6,
        "M");
    Matrices.printMatrix(
        Q,
        6,
        "Q");
    Matrices.printMatrix(
        R,
        6,
        "R");
    Vectors.printVector(
        p,
        6,
        "p");
    Matrices.printMatrix(
        QR,
        6,
        "QR");
    Matrices.printMatrix(
        QR2,
        6,
        "QR2");
    Matrices.printMatrix(
        M*P,
        6,
        "M*P");
    Matrices.printMatrix(
        M2,
        6,
        "M2");
    QR2 := QR*transpose(P);
    Matrices.printMatrix(
        QR2,
        6,
        "QR2");

    Modelica.Utilities.Streams.print(
      "Show results of QR factorization without pivoting:\n-----------------------------------------------------");
    (Q,R,tau,QR2) := Modelica_LinearSystems2.Math.Matrices.QR(
                              M);
    QR := Q*R;
    Matrices.printMatrix(
        Q,
        6,
        "Q");
    Matrices.printMatrix(
        R,
        6,
        "R");
    Matrices.printMatrix(
        QR,
        6,
        "QR");
    Matrices.printMatrix(
        QR2,
        6,
        "QR2");

  end exampleQR;

  function exampleSVD
    "Example for the usage of dgesdd and dgesvd lapack routines"
    import Modelica_LinearSystems2.Math.Matrices;
    import Modelica_LinearSystems2.Math.Vectors;
    input String fileName=DataDir + "m.mat"
                                annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="matrix file")));

  protected
    Real M[:,:]=Matrices.fromFile(fileName, "A");
    Real U1[size(M,1),size(M,1)];
    Real VT1[size(M,2),size(M,2)];
    Real sigma1[size(M,1)];

    Integer info;
  algorithm

    (sigma1,U1,VT1):=Matrices.LAPACK.dgesvd(M);

    Matrices.printMatrix(U1, 6, "U1");
    Matrices.printMatrix(VT1, 6, "VT1");
    Vectors.printVector(sigma1, 6, "sigma1");

  end exampleSVD;

  function care "solve a continuous algebraic Riccati equation"
    extends Modelica.Icons.Function;
    import Modelica_LinearSystems2.Math.Matrices;

  protected
    Real A[2,2]=[4,3; -9/2,-7/2];
    Real B[2,1]=[1; -1];
    Real R[1,1]=[1];
    Real Q[2,2]=[9,6; 6,4];
  public
    output Real X1[2,2]=Matrices.care(A, B, R, Q, false);
    output Real X2[2,2]=Matrices.care(A, B, R, Q, true);
    output Real X3[2,2]=[9*(1 + sqrt(2)),6*(1 + sqrt(2)); 6*(1 + sqrt(2)),4*(1 +
        sqrt(2))];

  algorithm
     Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
     Matrices.printMatrix(X1, 16, "X1");
     Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
     Matrices.printMatrix(X2, 16, "X2");
     Modelica.Utilities.Streams.print("Exact solution X3");
     Matrices.printMatrix(X3, 16, "X3");
  end care;

end Examples;

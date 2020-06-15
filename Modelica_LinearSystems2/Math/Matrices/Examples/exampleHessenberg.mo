within Modelica_LinearSystems2.Math.Matrices.Examples;
function exampleHessenberg
  "Example for the transformation of a matrix to upper Hessenberg form"
  import Modelica_LinearSystems2.Math.Matrices;
  input String fileName=DataDir + "m.mat"
    "Name of file where the matrix is saved"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="Open file with matrix")));
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
  (H,V,tau) := Modelica.Math.Matrices.Utilities.toUpperHessenberg(
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
  H := Modelica.Math.Matrices.hessenberg(M);
  Math.Matrices.printMatrix(
      H,
      6,
      "H_hess_lapack");

end exampleHessenberg;

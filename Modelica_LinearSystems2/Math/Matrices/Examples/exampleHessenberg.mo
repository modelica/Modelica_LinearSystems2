within Modelica_LinearSystems2.Math.Matrices.Examples;
function exampleHessenberg
  "Example for the transformation of a matrix to upper Hessenberg form"

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
  Modelica.Math.Matrices.toString(M, "M", 6);
  Modelica.Math.Matrices.toString(H, "Hu", 6);

  H := Modelica_LinearSystems2.Math.Matrices.Internal.hessenberg2(M, "l");
  Modelica.Math.Matrices.toString(H, "Hl", 6);

  (H,V,tau) := Modelica.Math.Matrices.Utilities.toUpperHessenberg(M, 1, size(M, 1));
  Modelica.Math.Matrices.toString(H, "H_dgehrd", 6);
  Modelica.Math.Matrices.toString(V, "V_dgehrd", 6);
  Modelica.Math.Vectors.toString(tau, "tau", 6);

  Q := Modelica_LinearSystems2.Math.Matrices.orthogonalQ(V, tau, 1, size(V, 1));
  Modelica.Math.Matrices.toString(Q*H*transpose(Q), "Q*H*Q'", 6);

  H := Modelica.Math.Matrices.hessenberg(M);
  Modelica.Math.Matrices.toString(H, "H_hess_lapack", 6);

end exampleHessenberg;

within Modelica_LinearSystems2.Math.Matrices.Examples;
function exampleQR
  "Example for the usage of QR2-function, QR factorization with columns pivoting"

  input String fileName=DataDir + "m.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="matrix file")));
  input String matrixName="A" "Name of the matrix";
protected
  Real M[:,:]=Modelica_LinearSystems2.Math.Matrices.fromFile(fileName, matrixName);
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
  Modelica.Math.Matrices.toString(M,"M",6);
  Modelica.Math.Matrices.toString(Q,"Q",6);
  Modelica.Math.Matrices.toString(R,"R",6);
  Modelica.Math.Vectors.toString(p, "p", 6);
  Modelica.Math.Matrices.toString(QR,"QR",6);
  Modelica.Math.Matrices.toString(QR2,"QR2",6);
  Modelica.Math.Matrices.toString(M*P,"M*P",6);
  Modelica.Math.Matrices.toString(M2,"M2",6);
  QR2 := QR*transpose(P);
  Modelica.Math.Matrices.toString(QR2,"QR2",6);

  Modelica.Utilities.Streams.print(
    "Show results of QR factorization without pivoting:\n-----------------------------------------------------");
  (Q,R,tau,QR2) := Modelica_LinearSystems2.Math.Matrices.QR(M);
  QR := Q*R;
  Modelica.Math.Matrices.toString(Q,"Q",6);
  Modelica.Math.Matrices.toString(R,"R",6);
  Modelica.Math.Matrices.toString(QR,"QR",6);
  Modelica.Math.Matrices.toString(QR2,"QR2",6);

end exampleQR;

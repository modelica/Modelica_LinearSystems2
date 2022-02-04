within Modelica_LinearSystems2.Math.Matrices.Examples;
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

  (sigma1,U1,VT1):=Modelica.Math.Matrices.LAPACK.dgesvd(M);

  Matrices.printMatrix(U1, 6, "U1");
  Matrices.printMatrix(VT1, 6, "VT1");
  Vectors.printVector(sigma1, 6, "sigma1");

end exampleSVD;

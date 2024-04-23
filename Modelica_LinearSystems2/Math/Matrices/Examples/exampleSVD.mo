within Modelica_LinearSystems2.Math.Matrices.Examples;
function exampleSVD
  "Example for the usage of dgesdd and dgesvd lapack routines"

  input String fileName=DataDir + "m.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="matrix file")));

protected
  Real M[:,:]=Modelica_LinearSystems2.Math.Matrices.fromFile(fileName, "A");
  Real U1[size(M,1),size(M,1)];
  Real VT1[size(M,2),size(M,2)];
  Real sigma1[size(M,1)];

  Integer info;
algorithm

  (sigma1,U1,VT1):=Modelica.Math.Matrices.LAPACK.dgesvd(M);

  Modelica.Math.Matrices.toString(U1, "U1", 6);
  Modelica.Math.Matrices.toString(VT1, "VT1", 6);
  Modelica.Math.Vectors.toString(sigma1, "sigma1", 6);

end exampleSVD;

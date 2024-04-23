within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
encapsulated function schur
  "Pole assignment design algorithm for multi input systems"
  import Modelica;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:]=[1,2,3,4;0,2,5,6;0,0,-2,-1;0,0,0,-3];
  input Real alpha=0;
  output Real T[:,:];
  output Real Z[:,:];
  output Real alphaReal[size(A, 1)]
    "Real part of eigenvalue=alphaReal+i*alphaImag";
  output Real alphaImag[size(A, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

algorithm
  (T,Z,alphaReal,alphaImag) := MatricesMSL.realSchur(A);
  Matrices.printMatrix(T,6,"T1");
  Modelica.Math.Vectors.toString(alphaReal, "ar1", 6);

// reorder real Schur form according to alpha
   (T,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSFc(
       T,
       identity(size(A, 1)),
       alphaReal,
       alphaImag,
       alpha);
  Matrices.printMatrix(T,6,"T2");
  Modelica.Math.Vectors.toString(alphaReal, "ar1", 6);
end schur;

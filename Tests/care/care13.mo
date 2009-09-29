within Modelica_LinearSystems2.Tests.care;
function care13 "Example 13  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e-6;
  Real A[4,4]=[0.0,  0.4,        0.0,       0.0;
               0.0,  0.0,        0.345,     0.0;
               0.0, -0.524/eps, -0.465/eps, 0.262/eps;
               0.0,  0.0,        0.0,      -1/eps];

  Real B[4,1]=[0; 0; 0; 1/eps];
  Real R[1,1]=[1];
  Real Q[4,4]=[1, 0, 0, 0; 0, 0, 0, 0; 0, 0, 1, 0; 0, 0, 0, 0];

public
  output Real X1[4,4];
  output Real X2[4,4];
  output Real X3[4,4]=[7.3840246663292781e+000,   5.9047620549742454e+000,   3.9930823439894581e-006,   9.9999999994873053e-007;
                       5.9047620549742454e+000,   7.1516063009772344e+000,   3.7996988588428911e-006,   8.6123471821676468e-007;
                       3.9930823439894581e-006,   3.7996988588428911e-006,   1.0402935802354235e-006,   1.8035961902063092e-007;
                       9.9999999994873053e-007,   8.6123471821676468e-007,   1.8035961902063092e-007,   4.6187574179129038e-008];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);

   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("MATLAB solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care13;

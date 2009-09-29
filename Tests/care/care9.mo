within Modelica_LinearSystems2.Tests.care;
function care9 "Example 9  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e-6;
  Real h=sqrt(1+2*eps);
  Real A[2,2]=[0.0, eps; 0.0, 0.0];
  Real B[2,1]=[0.0; 1.0];
  Real R[1,1]=[1];
  Real Q[2,2]=identity(2);

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=[h/eps, 1; 1, h];
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care9;

within Modelica_LinearSystems2.Tests.care;
function care10 "Example 10  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e-7;
  Real h=1+eps;
  Real A[2,2]=[h, 1.0; 1.0, h];
  Real B[2,2]=identity(2);
  Real R[2,2]=identity(2);
  Real Q[2,2]=eps^2*identity(2);
  Real x11=0.5*(2*h+sqrt(2*h^2+2)+eps*sqrt(2));
  Real x12=x11/(x11-h);

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=[x11, x12; x12, x11];
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care10;

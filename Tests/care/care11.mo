within Modelica_LinearSystems2.Tests.care;
function care11 "Example 11  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=0;
  Real A[2,2]=[3-eps, 1.0; 4.0, 2-eps];
  Real B[2,1]=[1.0; 1.0];
  Real R[1,1]=[1];
  Real Q[2,2]=[4*eps-11, 2*eps-5; 2*eps-5, 2*eps-2];

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=[2, 1; 1, 1];
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care11;

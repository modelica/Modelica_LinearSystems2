within Modelica_LinearSystems2.Math.Matrices.Examples;
function care "Solve a continuous algebraic Riccati equation"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real A[2,2]=[4,3; -9/2,-7/2];
  Real B[2,1]=[1; -1];
  Real R[1,1]=[1];
  Real Q[2,2]=[9,6; 6,4];
public
  output Real X1[2,2]=Matrices.care(A, B, R, Q, false);
  output Real X2[2,2]=Matrices.care(A, B, R, Q, true);
  output Real X3[2,2]=[9*(1 + sqrt(2)),6*(1 + sqrt(2)); 6*(1 + sqrt(2)),4*(1 +
      sqrt(2))];

algorithm
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");
end care;

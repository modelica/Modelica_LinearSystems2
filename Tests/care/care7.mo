within Modelica_LinearSystems2.Tests.care;
function care7 "Example 7  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1.0;
  Real A[2,2]=[1.0, 0.0; 0.0, -2.0];
  Real B[2,1]=[eps; 0.0];
  Real R[1,1]=[1];
  Real C[1,2]=[1, 1];
  Real Qtilde[1,1]=[1];
  Real Q[2,2]=transpose(C)*Qtilde*C;
  Real h;
public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   h := sqrt(1+eps^2);
   X3:=[(1+h)/(eps^2), 1/(2+h); 1/(2+h), 0.25*(1-(eps^2)/((2+h)^2))];
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");
   eps := 1e-6;
   B:=[eps; 0.0];
   h := sqrt(1+eps^2);
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=[(1+h)/(eps^2), 1/(2+h); 1/(2+h), 0.25*(1-(eps^2)/(2+h))];
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care7;

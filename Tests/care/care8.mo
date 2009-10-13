within Modelica_LinearSystems2.Tests.care;
function care8 "Example 8  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e-8;
  Real A[2,2]=[-0.1, 0.0; 0.0, -0.02];
  Real B[2,2]=[0.1, 0.0; 0.001, 0.01];
  Real R[2,2]=[1+eps, 1.0; 1.0, 1.0];
  Real C[1,2]=[10, 100];
  Real Qtilde[1,1]=[1];
  Real Q[2,2]=transpose(C)*Qtilde*C;
  Real h;

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2]=[7.4700062938350172e+001,   8.2995600931349986e+002;
                       8.2995600931349986e+002,   9.2213602958269603e+003];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   h := sqrt(1+eps^2);
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("MATLAB solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care8;

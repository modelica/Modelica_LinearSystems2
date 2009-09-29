within Modelica_LinearSystems2.Tests.care;
function care12 "Example 12  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e6;
  Real V[3,3]=[1/3, -2/3, -2/3; -2/3, 1/3, -2/3; -2/3, -2/3, 1/3];
  Real A0[3,3]=eps*[1, 0, 0; 0, 2, 0; 0, 0, 3];
  Real Q0[3,3]=[1/eps, 0, 0; 0, 1, 0; 0, 0, eps];
  Real A[3,3]=V*A0*V;
  Real B[3,3]=identity(3);
  Real R[3,3]=eps*identity(3);
  Real Q[3,3]=V*Q0*V;
  Real x1=eps^2+sqrt(eps^4+1);
  Real x2=2*eps^2+sqrt(4*eps^4+eps);
  Real x3=3*eps^2+sqrt(9*eps^4+eps^2);
  Real D[3,3]=[x1, 0, 0; 0, x2, 0; 0, 0, x3];

public
  output Real X1[3,3];
  output Real X2[3,3];
  output Real X3[3,3];

algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=V*D*V;
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("Exact solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care12;

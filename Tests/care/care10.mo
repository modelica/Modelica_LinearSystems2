within Modelica_LinearSystems2.Tests.care;
function care10 "Example 10  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;
  input String outputFile = "";

protected
  Real eps=1e-7;
  Real h=1+eps;
  Real A[2,2]=[h, 1.0; 1.0, h];
  Real B[2,2]=identity(2);
  Real R[2,2]=identity(2);
  Real G[2,2]=B*transpose(B);
  Real Q[2,2]=eps^2*identity(2);
  Real x11=0.5*(2*h+sqrt(2*h^2+2)+eps*sqrt(2));
  Real x12=x11/(x11-h);
  Real Qr1[2,2];
  Real Qr2[2,2];
  Real Qr3[2,2];
  Real deltaQ1;
  Real deltaQ2;
  Real deltaQ3;
  Real H[4,4]=[A,-G; -Q,-transpose(A)];
  Real condH=Modelica_LinearSystems2.Math.Matrices.conditionNumber(
                                               H);
  Real normH=Matrices.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];
  output Real ku1;
  output Real ku2;
  output Real ku3;

algorithm
  X1:=Matrices.care(A, B, R, Q, false);
  X2:=Matrices.care(A, B, R, Q, true);
  X3:=[x11, x12; x12, x11];
  Qr1 := X1*G*X1-transpose(A)*X1-X1*A;
  Qr2 := X2*G*X2-transpose(A)*X2-X2*A;
  Qr3 := X3*G*X3-transpose(A)*X3-X3*A;
  deltaQ1 := Modelica.Math.Matrices.norm(Q-Qr1)/Modelica.Math.Matrices.norm(Q);
  deltaQ2 := Modelica.Math.Matrices.norm(Q-Qr2)/Modelica.Math.Matrices.norm(Q);
  deltaQ3 := Modelica.Math.Matrices.norm(Q-Qr3)/Modelica.Math.Matrices.norm(Q);

  ku1:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X1);
  ku2:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X2);
  ku3:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X3);

  condX1 := Modelica_LinearSystems2.Math.Matrices.conditionNumber(X1);
  normX1 := Matrices.norm(X1, 2);
  condX2 := Modelica_LinearSystems2.Math.Matrices.conditionNumber(X2);
  normX2 := Matrices.norm(X2, 2);
  condX3 := Modelica_LinearSystems2.Math.Matrices.conditionNumber(X3);
  normX3 := Matrices.norm(X3, 2);
  Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement",outputFile);
  Modelica.Utilities.Streams.print(Matrices.printMatrix(X1, 16, "X1"),outputFile);
  Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement",outputFile);
  Modelica.Utilities.Streams.print(Matrices.printMatrix(X2, 16, "X2"),outputFile);
  Modelica.Utilities.Streams.print("Exact solution X3",outputFile);
  Modelica.Utilities.Streams.print(Matrices.printMatrix(X3, 16, "X3"),outputFile);
  Modelica.Utilities.Streams.print("\n normH = " + String(normH),outputFile);
  Modelica.Utilities.Streams.print("\n condH = " + String(condH),outputFile);
  Modelica.Utilities.Streams.print("\n normX1 = " + String(normX1),outputFile);
  Modelica.Utilities.Streams.print("\n condX1 = " + String(condX1),outputFile);
  Modelica.Utilities.Streams.print("\n ku1 = " + String(ku1),outputFile);
  Modelica.Utilities.Streams.print("\n normX2 = " + String(normX2),outputFile);
  Modelica.Utilities.Streams.print("\n condX2 = " + String(condX2),outputFile);
  Modelica.Utilities.Streams.print("\n ku2 = " + String(ku2),outputFile);
  Modelica.Utilities.Streams.print("\n normX3 = " + String(normX3),outputFile);
  Modelica.Utilities.Streams.print("\n condX3 = " + String(condX3),outputFile);
  Modelica.Utilities.Streams.print("\n ku3 = " + String(ku3),outputFile);
  Modelica.Utilities.Streams.print("\n deltaQ1 = " + String(deltaQ1),outputFile);
  Modelica.Utilities.Streams.print("\n deltaQ2 = " + String(deltaQ2),outputFile);
  Modelica.Utilities.Streams.print("\n deltaQ3 = " + String(deltaQ3),outputFile);
end care10;

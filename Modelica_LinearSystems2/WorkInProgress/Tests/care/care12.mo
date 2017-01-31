within Modelica_LinearSystems2.WorkInProgress.Tests.care;
function care12 "Example 12  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;
  input String outputFile = "";

protected
  Real eps=1e6;
  Real V[3,3]=[1/3, -2/3, -2/3; -2/3, 1/3, -2/3; -2/3, -2/3, 1/3];
  Real A0[3,3]=eps*[1, 0, 0; 0, 2, 0; 0, 0, 3];
  Real Q0[3,3]=[1/eps, 0, 0; 0, 1, 0; 0, 0, eps];
  Real A[3,3]=V*A0*V;
  Real B[3,3]=identity(3);
  Real R[3,3]=eps*identity(3);
  Real Q[3,3]=V*Q0*V;
//  Real G[3,3]=B*Modelica.Math.Matrices.solve2(R,B);
  Real G[3,3]=1/eps*identity(3);
  Real x1=eps^2+sqrt(eps^4+1);
  Real x2=2*eps^2+sqrt(4*eps^4+eps);
  Real x3=3*eps^2+sqrt(9*eps^4+eps^2);
  Real D[3,3]=[x1, 0, 0; 0, x2, 0; 0, 0, x3];
  Real Qr1[3,3];
  Real Qr2[3,3];
  Real Qr3[3,3];
  Real deltaQ1;
  Real deltaQ2;
  Real deltaQ3;
  Real H[6,6]=[A,-G; -Q,-transpose(A)];
  Real condH=Modelica_LinearSystems2.Math.Matrices.conditionNumber(
                                               H);
  Real normH=Matrices.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;
  Real resX1;
  Real resX2;
public
  output Real X1[3,3];
  output Real X2[3,3];
  output Real X3[3,3];
  output Real ku1;
  output Real ku2;
  output Real ku3;

algorithm
  X1:=Matrices.care(A, B, R, Q, false);
  X2:=Matrices.care(A, B, R, Q, true);
  X3:=V*D*V;
  Qr1 := X1*G*X1-transpose(A)*X1-X1*A;
  Qr2 := X2*G*X2-transpose(A)*X2-X2*A;
  Qr3 := X3*G*X3-transpose(A)*X3-X3*A;
  deltaQ1 := Modelica.Math.Matrices.norm(Q-Qr1)/Modelica.Math.Matrices.norm(Q);
  deltaQ2 := Modelica.Math.Matrices.norm(Q-Qr2)/Modelica.Math.Matrices.norm(Q);
  deltaQ3 := Modelica.Math.Matrices.norm(Q-Qr3)/Modelica.Math.Matrices.norm(Q);
  resX1 := Modelica.Math.Matrices.norm(X1-X3)/Modelica.Math.Matrices.norm(X3);
  resX2 := Modelica.Math.Matrices.norm(X2-X3)/Modelica.Math.Matrices.norm(X3);

 ku1:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(
    A,
    Q,
    G,
    X1);
 ku2:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(
    A,
    Q,
    G,
    X2);
 ku3:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(
    A,
    Q,
    G,
    X3);

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
  Modelica.Utilities.Streams.print("Residum of solution X1: resX1 = "+String(resX1),outputFile);
  Modelica.Utilities.Streams.print("Residum of solution X2: resX2 = "+String(resX2),outputFile);

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
end care12;

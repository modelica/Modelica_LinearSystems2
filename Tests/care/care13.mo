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
  Real G[4,4]=B*transpose(B);
  Real Q[4,4]=[1, 0, 0, 0; 0, 0, 0, 0; 0, 0, 1, 0; 0, 0, 0, 0];
  Real Qr1[4,4];
  Real Qr2[4,4];
  Real deltaQ1;
  Real deltaQ2;
  Real H[8,8]=[A,-G; -Q,-transpose(A)];
  Real condH=Matrices.Internal.conditionNumber(H);
  Real normH=Matrices.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;

public
  output Real X1[4,4];
  output Real X2[4,4];
  output Real X3[4,4]=[7.3840246663292781e+000,   5.9047620549742454e+000,   3.9930823439894581e-006,   9.9999999994873053e-007;
                       5.9047620549742454e+000,   7.1516063009772344e+000,   3.7996988588428911e-006,   8.6123471821676468e-007;
                       3.9930823439894581e-006,   3.7996988588428911e-006,   1.0402935802354235e-006,   1.8035961902063092e-007;
                       9.9999999994873053e-007,   8.6123471821676468e-007,   1.8035961902063092e-007,   4.6187574179129038e-008];
  output Real ku1;
  output Real ku2;
  output Real ku3=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X3);
algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   ku1:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X1);
   ku2:=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X2);
  Qr1 := X1*G*X1-transpose(A)*X1-X1*A;
    Qr2 := X2*G*X2-transpose(A)*X2-X2*A;
    deltaQ1 := Modelica.Math.Matrices.norm(Q-Qr1)/Modelica.Math.Matrices.norm(Q);
    deltaQ2 := Modelica.Math.Matrices.norm(Q-Qr2)/Modelica.Math.Matrices.norm(Q);
   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("MATLAB solution X3");
   Matrices.printMatrix(X3, 16, "X3");
  condX1 := Matrices.Internal.conditionNumber(X1);
  normX1 := Matrices.norm(X1, 2);
  condX2 := Matrices.Internal.conditionNumber(X2);
  normX2 := Matrices.norm(X2, 2);
  condX3 := Matrices.Internal.conditionNumber(X3);
  normX3 := Matrices.norm(X3, 2);

  Modelica.Utilities.Streams.print("\n normH = " + String(normH));
  Modelica.Utilities.Streams.print("\n condH = " + String(condH));
  Modelica.Utilities.Streams.print("\n normX1 = " + String(normX1));
  Modelica.Utilities.Streams.print("\n condX1 = " + String(condX1));
  Modelica.Utilities.Streams.print("\n ku1 = " + String(ku1));
  Modelica.Utilities.Streams.print("\n normX2 = " + String(normX2));
  Modelica.Utilities.Streams.print("\n condX2 = " + String(condX2));
  Modelica.Utilities.Streams.print("\n ku2 = " + String(ku2));
  Modelica.Utilities.Streams.print("\n normX3 = " + String(normX3));
  Modelica.Utilities.Streams.print("\n condX3 = " + String(condX3));
  Modelica.Utilities.Streams.print("\n ku3 = " + String(ku3));
  Modelica.Utilities.Streams.print("\n deltaQ1 = " + String(deltaQ1));
  Modelica.Utilities.Streams.print("\n deltaQ2 = " + String(deltaQ2));

end care13;

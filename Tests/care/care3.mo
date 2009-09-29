within Modelica_LinearSystems2.Tests.care;
function care3 "Example 3 from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real A[4,4]=[0.0,   1.0,   0.0,   0.0;
               0.0,  -1.890,   3.900e-01,  -5.530;
               0.0,  -3.400e-02,  -2.980,   2.430;
               3.400e-02,  -1.100e-03,  -9.900e-01,  -2.100e-01];
  Real B[4,2]=[ 0.0,   0.0;
                3.600e-01,  -1.60;
               -9.500e-01,  -3.200e-02;
                3.000e-02,   0.0];
  Real R[2,2]=[1, 0; 0, 1];
  Real Q[4,4]=[2.313,   2.727,   6.880e-01,   2.300e-02;
               2.727,   4.271,  1.148,   3.230e-01;
               6.880e-01,   1.148,   3.130e-01,   1.020e-01;
               2.300e-02,   3.230e-01,   1.020e-01,   8.300e-02];
  Real G[4,4]=B*transpose(B);
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
  output Real X1[4,4]=Matrices.care(A, B, R, Q, false);
  output Real X2[4,4]=Matrices.care(A, B, R, Q, true);
  output Real X3[4,4]=[
  1.3238595718183990,               9.0153284952164081e-001,   5.4663403916715425e-001,  -1.7672385587639674;
  9.0153284952164081e-001,   9.6068122262991273e-001,   4.3342816873410484e-001,  -1.1989126854651069;
  5.4663403916715425e-001,   4.3342816873410484e-001,   4.6054882548934967e-001,  -1.3632873589876660;
 -1.7672385587639674,              -1.1989126854651069,     -1.3632873589876660,              4.4611816254580852];

  output Real ku1=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X1);
  output Real ku2=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X2);
  output Real ku3=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X3);

algorithm
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

end care3;

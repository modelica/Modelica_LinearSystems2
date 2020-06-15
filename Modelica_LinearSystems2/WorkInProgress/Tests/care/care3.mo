within Modelica_LinearSystems2.WorkInProgress.Tests.care;
function care3 "Example 3 from Benner benchmarks"
  extends Modelica.Icons.Function;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices;

  input String outputFile = "";

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
  Real condH=MatricesMSL.conditionNumber(H);
  Real normH=MatricesMSL.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;
  Real Qr1[4,4];
  Real Qr2[4,4];
  Real Qr3[4,4];
  Real deltaQ1;
  Real deltaQ2;
  Real deltaQ3;
  Real resX1;
  Real resX2;
public
  output Real X1[4,4]=Matrices.care(A, B, R, Q, false);
  output Real X2[4,4]=Matrices.care(A, B, R, Q, true);
  output Real X3[4,4]=[
  1.3238595718183990,               9.0153284952164081e-001,   5.4663403916715425e-001,  -1.7672385587639674;
  9.0153284952164081e-001,   9.6068122262991273e-001,   4.3342816873410484e-001,  -1.1989126854651069;
  5.4663403916715425e-001,   4.3342816873410484e-001,   4.6054882548934967e-001,  -1.3632873589876660;
 -1.7672385587639674,              -1.1989126854651069,     -1.3632873589876660,              4.4611816254580852];

  output Real ku1=Matrices.Internal.k_care_u(A, Q, G, X1);
  output Real ku2=Matrices.Internal.k_care_u(A, Q, G, X2);
  output Real ku3=Matrices.Internal.k_care_u(A, Q, G, X3);

algorithm
  Qr1 := X1*G*X1-transpose(A)*X1-X1*A;
  Qr2 := X2*G*X2-transpose(A)*X2-X2*A;
  Qr3 := X3*G*X3-transpose(A)*X3-X3*A;
  deltaQ1 := MatricesMSL.norm(Q-Qr1)/MatricesMSL.norm(Q);
  deltaQ2 := MatricesMSL.norm(Q-Qr2)/MatricesMSL.norm(Q);
  deltaQ3 := MatricesMSL.norm(Q-Qr3)/MatricesMSL.norm(Q);
  resX1 := MatricesMSL.norm(X1-X3)/MatricesMSL.norm(X3);
  resX2 := MatricesMSL.norm(X2-X3)/MatricesMSL.norm(X3);

  condX1 := MatricesMSL.conditionNumber(X1);
  normX1 := MatricesMSL.norm(X1, 2);
  condX2 := MatricesMSL.conditionNumber(X2);
  normX2 := MatricesMSL.norm(X2, 2);
  condX3 := MatricesMSL.conditionNumber(X3);
  normX3 := MatricesMSL.norm(X3, 2);

  Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement",outputFile);
  Modelica.Utilities.Streams.print(Matrices.printMatrix(X1, 16, "X1"),outputFile);
  Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement",outputFile);
  Modelica.Utilities.Streams.print(Matrices.printMatrix(X2, 16, "X2"),outputFile);
  Modelica.Utilities.Streams.print("MATLAB solution X3",outputFile);
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

//   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
//   Modelica.Utilities.Streams.print(Matrices.printMatrix(X1, 16, "X1"));
//   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
//   Modelica.Utilities.Streams.print(Matrices.printMatrix(X2, 16, "X2"));
//   Modelica.Utilities.Streams.print("MATLAB solution X3");
//   Modelica.Utilities.Streams.print(Matrices.printMatrix(X3, 16, "X3"));
//   Modelica.Utilities.Streams.print("\n normH = " + String(normH));
//   Modelica.Utilities.Streams.print("\n condH = " + String(condH));
//   Modelica.Utilities.Streams.print("\n normX1 = " + String(normX1));
//   Modelica.Utilities.Streams.print("\n condX1 = " + String(condX1));
//   Modelica.Utilities.Streams.print("\n ku1 = " + String(ku1));
//   Modelica.Utilities.Streams.print("\n normX2 = " + String(normX2));
//   Modelica.Utilities.Streams.print("\n condX2 = " + String(condX2));
//   Modelica.Utilities.Streams.print("\n ku2 = " + String(ku2));
//   Modelica.Utilities.Streams.print("\n normX3 = " + String(normX3));
//   Modelica.Utilities.Streams.print("\n condX3 = " + String(condX3));
//   Modelica.Utilities.Streams.print("\n ku3 = " + String(ku3));
//   Modelica.Utilities.Streams.print("\n deltaQ1 = " + String(deltaQ1));
//   Modelica.Utilities.Streams.print("\n deltaQ2 = " + String(deltaQ2));
//   Modelica.Utilities.Streams.print("\n deltaQ3 = " + String(deltaQ3));

end care3;

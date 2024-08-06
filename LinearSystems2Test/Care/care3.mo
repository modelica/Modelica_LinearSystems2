within LinearSystems2Test.Care;
function care3 "Example 3 from Benner benchmarks"
  extends Modelica.Icons.Function;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input String outputFile = "";
  output Boolean ok;

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
  ok := false;
  print("--  Test of " + getInstanceName() + " --");
  if not Modelica.Utilities.Strings.isEmpty(outputFile) then
    print("--  Test of " + getInstanceName() + " --", outputFile);
  end if;

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

  print("Solution X1 without subsequent Newton refinement", outputFile);
  print(MatricesMSL.toString(X1, "X1", 16), outputFile);
  print("Solution X2 with subsequent Newton refinement", outputFile);
  print(MatricesMSL.toString(X2, "X2", 16), outputFile);
  print("MATLAB solution X3", outputFile);
  print(MatricesMSL.toString(X3, "X3", 16), outputFile);
  print("Residum of solution X1: resX1 = "+String(resX1), outputFile);
  print("Residum of solution X2: resX2 = "+String(resX2), outputFile);

  print("\n normH = " + String(normH), outputFile);
  print("\n condH = " + String(condH), outputFile);
  print("\n normX1 = " + String(normX1), outputFile);
  print("\n condX1 = " + String(condX1), outputFile);
  print("\n ku1 = " + String(ku1), outputFile);
  print("\n normX2 = " + String(normX2), outputFile);
  print("\n condX2 = " + String(condX2), outputFile);
  print("\n ku2 = " + String(ku2), outputFile);
  print("\n normX3 = " + String(normX3), outputFile);
  print("\n condX3 = " + String(condX3), outputFile);
  print("\n ku3 = " + String(ku3), outputFile);
  print("\n deltaQ1 = " + String(deltaQ1), outputFile);
  print("\n deltaQ2 = " + String(deltaQ2), outputFile);
  print("\n deltaQ3 = " + String(deltaQ3), outputFile);

//   print("Solution X1 without subsequent Newton refinement");
//   print(MatricesMSL.toString(X1, "X1", 16));
//   print("Solution X2 with subsequent Newton refinement");
//   print(MatricesMSL.toString(X2, "X2", 16));
//   print("MATLAB solution X3");
//   print(MatricesMSL.toString(X3, "X3", 16));
//   print("\n normH = " + String(normH));
//   print("\n condH = " + String(condH));
//   print("\n normX1 = " + String(normX1));
//   print("\n condX1 = " + String(condX1));
//   print("\n ku1 = " + String(ku1));
//   print("\n normX2 = " + String(normX2));
//   print("\n condX2 = " + String(condX2));
//   print("\n ku2 = " + String(ku2));
//   print("\n normX3 = " + String(normX3));
//   print("\n condX3 = " + String(condX3));
//   print("\n ku3 = " + String(ku3));
//   print("\n deltaQ1 = " + String(deltaQ1));
//   print("\n deltaQ2 = " + String(deltaQ2));
//   print("\n deltaQ3 = " + String(deltaQ3));

  ok := true;
end care3;

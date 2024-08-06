within LinearSystems2Test.Care;
function care11 "Example 11  from Benner benchmarks"
  extends Modelica.Icons.Function;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input String outputFile = "";
  output Boolean ok;

protected
  Real eps=0;
  Real A[2,2]=[3-eps, 1.0; 4.0, 2-eps];
  Real B[2,1]=[1.0; 1.0];
  Real R[1,1]=[1];
  Real G[2,2]=B*transpose(B);
  Real Q[2,2]=[4*eps-11, 2*eps-5; 2*eps-5, 2*eps-2];
  Real Qr1[2,2];
  Real Qr2[2,2];
  Real Qr3[2,2];
  Real deltaQ1;
  Real deltaQ2;
  Real deltaQ3;
  Real H[4,4]=[A,-G; -Q,-transpose(A)];
  Real condH=MatricesMSL.conditionNumber(H);
  Real normH=MatricesMSL.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;
  Real resX1;
  Real resX2;

public
  output Real X1[2,2];
  output Real X2[2,2];
  output Real X3[2,2];
  output Real ku1;
  output Real ku2;
  output Real ku3;

algorithm
  ok := false;
  print("--  Test of " + getInstanceName() + " --");
  if not Modelica.Utilities.Strings.isEmpty(outputFile) then
    print("--  Test of " + getInstanceName() + " --", outputFile);
  end if;

   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);
   X3:=[2, 1; 1, 1];
//    ku1:=Matrices.Internal.k_care_u(A, Q, G, X1);
//    ku2:=Matrices.Internal.k_care_u(A, Q, G, X2);
//    ku3:=Matrices.Internal.k_care_u(A, Q, G, X3);

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
  print("Exact solution X3", outputFile);
  print(MatricesMSL.toString(X3, "X3", 16), outputFile);
  print("Residum of solution X1: resX1 = "+String(resX1), outputFile);
  print("Residum of solution X2: resX2 = "+String(resX2), outputFile);

  print("\n normH = " + String(normH), outputFile);
  print("\n condH = " + String(condH), outputFile);
  print("\n normX1 = " + String(normX1), outputFile);
  print("\n condX1 = " + String(condX1), outputFile);
//  print("\n ku1 = " + String(ku1), outputFile);
  print("\n normX2 = " + String(normX2), outputFile);
  print("\n condX2 = " + String(condX2), outputFile);
//  print("\n ku2 = " + String(ku2), outputFile);
  print("\n normX3 = " + String(normX3), outputFile);
  print("\n condX3 = " + String(condX3), outputFile);
//  print("\n ku3 = " + String(ku3), outputFile);
  print("\n deltaQ1 = " + String(deltaQ1), outputFile);
  print("\n deltaQ2 = " + String(deltaQ2), outputFile);
  print("\n deltaQ3 = " + String(deltaQ3), outputFile);

  ok := true;
end care11;

within Modelica_LinearSystems2.Tests.care;
function care4 "Example 4 from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real A[8,8]=[  -9.910e-01,   5.290e-01,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0;
                  5.220e-01,  -1.051,   5.960e-01,   0.0,   0.0, 0.0,   0.0,   0.0;
                  0.0,   5.220e-01,  -1.118,   5.960e-01,   0.0,  0.0,   0.0,   0.0;
                  0.0,   0.0,   5.220e-01,  -1.548,   7.180e-01,  0.0,   0.0,   0.0;
                  0.0,   0.0,   0.0,   9.220e-01,  -1.640, 7.990e-01,   0.0,   0.0;
                  0.0,   0.0,   0.0,   0.0,   9.220e-01, -1.721,   9.010e-01,   0.0;
                  0.0,   0.0,   0.0,   0.0,   0.0, 9.220e-01,  -1.823,   1.021;
                  0.0,   0.0,   0.0,   0.0,   0.0, 0.0,   9.220e-01,  -1.943];

  Real B[8,2]=1e-3*[3.840,  -2.880;
                    4.0,  -3.040;
                   37.60,  -2.800;
                   3.080,  -2.320;
                   2.360,  -3.320;
                   2.880,  -3.820;
                   3.080,  -4.120;
                   3.0,  -3.960];

  Real R[2,2]=[1, 0; 0, 1];
  Real Q[8,8]=[  1.0,   0.0,   0.0,   0.0,   5.0e-01, 0.0,   0.0,   1.0e-01;
                 0.0,   1.0,   0.0,   0.0,   1.0e-01, 0.0,   0.0,   0.0;
                 0.0,   0.0,   1.0,   0.0,   0.0, 5.0e-01,   0.0,   0.0;
                 0.0,   0.0,   0.0,   1.0,   0.0, 0.0,   0.0,   0.0;
                 5.0e-01,   1.0e-01,   0.0,   0.0,   1.0e-01, 0.0,   0.0,   0.0;
                 0.0,   0.0,   5.0e-01,   0.0,   0.0, 1.0e-01,   0.0,   0.0;
                 0.0,   0.0,   0.0,   0.0,   0.0, 0.0,   1.0e-01,   0.0;
                 1.0e-01,   0.0,   0.0,   0.0,   0.0, 0.0,   0.0,   1.0e-01];
  Real G[8,8]=B*transpose(B);
  Real H[16,16]=[A,-G; -Q,-transpose(A)];
  Real condH=Matrices.Internal.conditionNumber(H);
  Real normH=Matrices.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;
  Real Qr1[8,8];
  Real Qr2[8,8];
  Real deltaQ1;
  Real deltaQ2;
public
  output Real X1[8,8]=Matrices.care(A, B, R, Q, false);
  output Real X2[8,8]=Matrices.care(A, B, R, Q, true);

  output Real X3[8,8]=[8.9189179333315727e-001,   7.3663786199205217e-001,   6.0233803898550886e-001,   5.2119369438941698e-001,   5.9291413705093765e-001,   3.4884170426731775e-001,   2.1987251788724343e-001,   1.4147828732718845e-001;
                       7.3663786199205217e-001,   1.3795283147296966e+000,   1.0764852472252566e+000,   8.0385954540875648e-001,   7.0054708500170060e-001,   5.1907694010695571e-001,   3.3482970113161148e-001,   1.7436126307010602e-001;
                       6.0233803898550886e-001,   1.0764852472252566e+000,   1.4919719225454409e+000,   1.0138055847275820e+000,   8.0136227954219730e-001,   7.4348008693508494e-001,   4.1924715019049130e-001,   2.0305623101755202e-001;
                       5.2119369438941698e-001,   8.0385954540875648e-001,   1.0138055847275820e+000,   1.1487703436604533e+000,   7.3268221838065584e-001,   5.3128128249432527e-001,   3.4101748640719931e-001,   1.7321992448577761e-001;
                       5.9291413705093765e-001,   7.0054708500170060e-001,   8.0136227954219730e-001,   7.3268221838065584e-001,   5.9205254501700877e-001,   4.2933248637928978e-001,   2.8472012561858129e-001,   1.4763528881298796e-001;
                       3.4884170426731775e-001,   5.1907694010695571e-001,   7.4348008693508494e-001,   5.3128128249432527e-001,   4.2933248637928978e-001,   3.5530939867621736e-001,   2.3769595202736848e-001,   1.2407312100669236e-001;
                       2.1987251788724343e-001,   3.3482970113161148e-001,   4.1924715019049130e-001,   3.4101748640719931e-001,   2.8472012561858129e-001,   2.3769595202736848e-001,   1.9654065110986793e-001,   1.0236228817160951e-001;
                       1.4147828732718845e-001,   1.7436126307010602e-001,   2.0305623101755202e-001,   1.7321992448577761e-001,   1.4763528881298796e-001,   1.2407312100669236e-001,   1.0236228817160951e-001,   7.9489693942802808e-002];
  output Real ku1=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X1);
  output Real ku2=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X2);
  output Real ku3=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X3);

algorithm
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

end care4;

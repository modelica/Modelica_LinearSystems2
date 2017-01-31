within Modelica_LinearSystems2.WorkInProgress.Tests.care;
function care5 "Example 5 from Benner benchmarks"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;
  input String outputFile = "";

protected
  Real A[9,9]=[-4.019,    5.120,   0.0,   0.0,  -2.082,   0.0,  0.0,   0.0,   8.700e-01;
               -3.460e-01,    9.860e-01,   0.0,   0.0,  -2.340,   0.0,  0.0,   0.0,   9.700e-01;
               -7.909,   1.5407e+01,  -4.069,   0.0,  -6.450,   0.0,  0.0,   0.0,   2.680;
               -2.1816e+01,  3.5606e+01,  -3.390e-01,  -3.870,  -1.780e+01,   0.0,  0.0,   0.0,   7.390;
               -6.0196e+01,  9.8188e+01,  -7.907,   3.400e-01,  -5.3008e+01,  0.0,  0.0,   0.0,   2.040e+01;
                0.0,    0.0,   0.0,   0.0,   9.400e+01,  -1.472e+02,  0.0,   5.320e+01,   0.0;
                0.0,    0.0,   0.0,   0.0,   0.0,   9.400e+01, -1.472e+02,   0.0,   0.0;
                0.0,    0.0,   0.0,   0.0,   0.0,   1.280e+01,  0.0,  -3.160e+01,   0.0;
                0.0,    0.0,   0.0,   0.0,   1.280e+01,   0.0,  0.0,   1.880e+01,  -3.160e+01];

  Real B[9,3]=[  0.0100,   -0.0110,   -0.1510;
                 0.0030,   -0.0210,         0;
                 0.0090,   -0.0590,         0;
                 0.0240,   -0.1620,         0;
                 0.0680,   -0.4450,         0;
                 0,         0,         0;
                 0,         0,         0;
                 0,         0,         0;
                 0,         0,         0];

  Real R[3,3]=identity(3);
  Real Q[9,9]=identity(9);
  Real G[9,9]=B*transpose(B);
  Real H[18,18]=[A,-G; -Q,-transpose(A)];
  Real condH=Modelica_LinearSystems2.Math.Matrices.conditionNumber(
                                               H);
  Real normH=Matrices.norm(H, 2);
  Real condX1;
  Real normX1;
  Real condX2;
  Real normX2;
  Real condX3;
  Real normX3;
  Real Qr1[9,9];
  Real Qr2[9,9];
  Real Qr3[9,9];
  Real deltaQ1;
  Real deltaQ2;
  Real deltaQ3;

public
  output Real X1[9,9]=Matrices.care(A, B, R, Q, false);
  output Real X2[9,9]=Matrices.care(A, B, R, Q, true);

  output Real X3[9,9]=[1.8813417073612682e+000,   3.9350887811241636e-001,   3.0066821723622605e-001,  -5.6646472831415537e-002,  -1.3921476885988587e-001,  -4.0430720912527553e-003,  -2.0593909138650463e-004,  -3.2846155476217977e-002,  -2.5012667587779286e-002;
                       3.9350887811241636e-001,   2.4445142440433059e+000,   2.9761009646012448e-001,   1.1904010124814599e-001,  -1.4000314732846808e-001,   8.6371528201907403e-003,   3.3855145579537928e-004,   7.6337899432566164e-002,   6.9126140585845194e-002;
                       3.0066821723622605e-001,   2.9761009646012448e-001,   2.5625289349920533e-001,   3.5608619331945784e-002,  -7.0291882601897213e-002,  -2.9479293666103855e-004,  -2.7109320576564369e-005,  -1.5876949090013865e-003,   6.6183235965156960e-004;
                      -5.6646472831415537e-002,   1.1904010124814599e-001,   3.5608619331945784e-002,   1.2531244600262245e-001,  -4.4098271439652045e-002,   2.8906009190787904e-004,   1.1585504034027231e-006,   3.3312522262940710e-003,   5.3504545840022469e-003;
                      -1.3921476885988587e-001,  -1.4000314732846808e-001,  -7.0291882601897213e-002,  -4.4098271439652045e-002,   5.0544674045016849e-002,   2.7384271164510657e-003,   5.1877285822234109e-004,   6.6383915326013159e-003,   5.1811830100945500e-003;
                      -4.0430720912527553e-003,   8.6371528201907403e-003,  -2.9479293666103855e-004,   2.8906009190787904e-004,   2.7384271164510657e-003,   4.4145119317543710e-003,   1.0990456156319732e-003,   3.6333616685830985e-003,   1.1417260453168568e-003;
                      -2.0593909138650463e-004,   3.3855145579537928e-004,  -2.7109320576564369e-005,   1.1585504034027231e-006,   5.1877285822234109e-004,   1.0990456156319732e-003,   3.3967389367722083e-003,   3.3327489291592403e-004,   5.9658514141889494e-005;
                      -3.2846155476217977e-002,   7.6337899432566164e-002,  -1.5876949090013865e-003,   3.3312522262940710e-003,   6.6383915326013159e-003,   3.6333616685830985e-003,   3.3327489291592403e-004,   2.8282052511283625e-002,   1.0661764834431937e-002;
                      -2.5012667587779286e-002,   6.9126140585845194e-002,   6.6183235965156960e-004,   5.3504545840022469e-003,   5.1811830100945500e-003,   1.1417260453168568e-003,   5.9658514141889494e-005,   1.0661764834431937e-002,   2.1907727243964663e-002];

  output Real ku1=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X1);
  output Real ku2=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X2);
  output Real ku3=Modelica_LinearSystems2.Math.Matrices.Internal.k_care_u(A, Q, G, X3);
algorithm
  Qr1 := X1*G*X1-transpose(A)*X1-X1*A;
  Qr2 := X2*G*X2-transpose(A)*X2-X2*A;
  Qr3 := X3*G*X3-transpose(A)*X3-X3*A;
  deltaQ1 := Modelica.Math.Matrices.norm(Q-Qr1)/Modelica.Math.Matrices.norm(Q);
  deltaQ2 := Modelica.Math.Matrices.norm(Q-Qr2)/Modelica.Math.Matrices.norm(Q);
  deltaQ3 := Modelica.Math.Matrices.norm(Q-Qr3)/Modelica.Math.Matrices.norm(Q);
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
  Modelica.Utilities.Streams.print("MATLAB solution X3",outputFile);
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

end care5;

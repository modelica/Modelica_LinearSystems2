within Modelica_LinearSystems2.Tests.care;
function care14 "Example 14  from Benner benchmarks"
//needs another algorithm; ev of Hamiltonian close to imaginary axis
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices;

protected
  Real eps=1e-6;
  Real A[4,4]=[-eps,  1.0,   0.0,  0.0;
               -1.0,  -eps,  0.0,  0.0;
                0.0,   0.0,  eps,  1.0;
                0.0,  0.0,  -1.0,  eps];

  Real B[4,1]=[1; 1; 1; 1];
  Real R[1,1]=[1];
  Real C[1,4]=[1, 1, 1, 1];
  Real Q[4,4]=transpose(C)*C;

public
  output Real X1[4,4];
  output Real X2[4,4];
  output Real X3[4,4]=[1.0010401156071251e+000,  -3.0124557362509741e-016,  -1.0421176888569740e-003,   1.0010411167310994e-006;
 -3.0124557362509741e-016,   1.0010421176893556e+000,  -1.0010431186887811e-006,  -1.0421176888547718e-003;
 -1.0421176888569740e-003,  -1.0010431186887811e-006,   1.0010441197755955e+000,  -9.2861916073383455e-017;
  1.0010411167310994e-006,  -1.0421176888547718e-003,  -9.2861916073383455e-017,   1.0010421176893556e+000];
algorithm
   X1:=Matrices.care(A, B, R, Q, false);
   X2:=Matrices.care(A, B, R, Q, true);

   Modelica.Utilities.Streams.print("Solution X1 without subsequent Newton refinement");
   Matrices.printMatrix(X1, 16, "X1");
   Modelica.Utilities.Streams.print("Solution X2 with subsequent Newton refinement");
   Matrices.printMatrix(X2, 16, "X2");
   Modelica.Utilities.Streams.print("MATLAB solution X3");
   Matrices.printMatrix(X3, 16, "X3");

end care14;

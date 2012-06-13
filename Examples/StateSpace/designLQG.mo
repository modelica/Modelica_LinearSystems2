within Modelica_LinearSystems2.Examples.StateSpace;
function designLQG "Example for LQG controller design"

  import Modelica_LinearSystems2.StateSpace;

 extends Modelica.Icons.Function;

input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
      A=[-0.02, 0.005, 2.4,  -32; -0.14,  0.44,  -1.3,  -30; 0,  0.018,  -1.6,  1.2; 0, 0, 1, 0],
      B=[0.14,  -0.12; 0.36, -8.6; 0.35, 0.009; 0, 0],
      C=[0, 1, 0, 0; 0, 0, 0, 57.3],
      D=[0,0; 0,0]);

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
   annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
                                                                           annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                     caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"  annotation(Dialog(group="system data definition",enable = systemOnFile));

  output Boolean ok;

protected
 StateSpace ss=if systemOnFile then Modelica_LinearSystems2.StateSpace.Import.fromFile(
                                                         fileName) else ssi;
 Real Q[:,:] = transpose(ss.C)*ss.C " state weighting matrix";
 Real R[:,:] = identity(2) " input weighting matrix";
 Real V[:,:] = identity(2) " covariance output noise matrix";
 Real W[:,:] = ss.B*transpose(ss.B) " covariance input noise matrix";
 Real Kc[size(ss.B, 2),size(ss.A, 1)] "Controller feedback gain matrix";
 Real Kf[size(ss.A, 1),size(ss.C, 1)] "Kalman feedback gain matrix";
 StateSpace sslgq(
    redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
    redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
    redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
    redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

algorithm
 (Kc, Kf, sslgq) := StateSpace.Design.lqg(ss, Q, R, V, W);

 Modelica_LinearSystems2.Math.Matrices.printMatrix(Kc,6,"Kc");
 Modelica_LinearSystems2.Math.Matrices.printMatrix(Kf,6,"Kf");
 Modelica.Utilities.Streams.print(String(sslgq));
 ok := true;

  annotation (Documentation(info="<html>
This example demonstrates the computatrion of a LQG controllerr by calling function <b>StateSpace.Design.lqg()</b> 
</html>"));
end designLQG;

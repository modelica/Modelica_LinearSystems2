within Modelica_LinearSystems2.Examples.StateSpace;
function designLQR "Description"

  import Modelica_LinearSystems2.StateSpace;

 extends Modelica.Icons.Function;

input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
      A=[0, 1, 0, 0; 0, 0, 39.2, 0; 0, 0, 0, 1; 0, 0, 49, 0],
      B=[0; 1; 0; 1],
      C=[1, 0, 0, 0],
      D=[0]);

  input Boolean systemOnFile=false
    "true, if state space system is defined on file"
   annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
                                                                           annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                     caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"  annotation(Dialog(group="system data definition",enable = systemOnFile));

  output Boolean ok;

protected
  StateSpace ss=if systemOnFile then Modelica_LinearSystems2.StateSpace.Import.fromFile(
                                                          fileName) else ssi;
  Real Q[:,:]=identity(4) " state weighting matrix";
  Real R[:,:]=identity(1) " input weighting matrix";
  Real K[size(ss.B, 2),size(ss.A, 1)] "Controller feedback gain matrix";

algorithm
  K := StateSpace.Design.lqr(ss, Q, R);

  Modelica_LinearSystems2.Math.Matrices.printMatrix(K, 6, "K");

  ok := true;

end designLQR;

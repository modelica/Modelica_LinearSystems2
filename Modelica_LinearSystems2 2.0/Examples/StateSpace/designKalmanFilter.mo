within Modelica_LinearSystems2.Examples.StateSpace;
function designKalmanFilter "Description"
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
 extends Modelica.Icons.Function;

 input StateSpace ss=Modelica_LinearSystems2.StateSpace(TransferFunction(n={1}, d={1,2,3,4}));
protected
 Real Q[:,:]=identity(3);
 Real R[:,:]=[1];

public
 output Real L[:,:];
 output StateSpace kss(
   redeclare Real A[size(ss.A, 1),size(ss.A, 1)],
   redeclare Real B[size(ss.B, 1),size(ss.B, 2) + size(ss.C, 1)],
   redeclare Real C[size(ss.A, 1),size(ss.A, 2)],
   redeclare Real D[size(ss.A, 1),size(ss.B, 2) + size(ss.C,1)]);

algorithm
 (L,kss) := StateSpace.Design.kalmanFilter(
   ss,
   Q,
   R);

end designKalmanFilter;

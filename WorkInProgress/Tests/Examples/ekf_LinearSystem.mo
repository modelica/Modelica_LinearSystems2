within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function ekf_LinearSystem "Linear state space function"

//  import Modelica.Math.Matrices;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  extends Modelica_LinearSystems2.DiscreteStateSpace.Internal.ekfSystemBase;

  input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method";

//   parameter StateSpace ss=StateSpace(A=[-1.2822, 0, 0.98, 0; 0, 0, 1, 0; -5.4293, 0, -1.8366, 0; -128.2, 128.2, 0, 0],
//                                      B=[-0.3; 0; -17; 0], C=[0, 1, 0, 0; 0, 0, 0, 1; -128.2, 128.2, 0, 0], D=[0; 0; 0]);
//

  parameter Real A[:,:]=[-1.2822, 0, 0.98, 0; 0, 0, 1, 0; -5.4293, 0, -1.8366, 0; -128.2, 128.2, 0, 0];
  parameter Real B[:,:]=[-0.3; 0; -17; 0];
  parameter Real C[:,:]=[0, 1, 0, 0; 0, 0, 0, 1; -128.2, 128.2, 0, 0];
  parameter Real D[:,:]=[0; 0; 0];

protected
     DiscreteStateSpace dss(
       redeclare Real A[size(A, 1),size(A, 2)],
       redeclare Real B[size(B, 1),size(B, 2)],
       redeclare Real C[size(C, 1),size(C, 2)],
       redeclare Real D[size(D, 1),size(D, 2)],
       redeclare Real B2[size(B, 1),size(B, 2)]);
//  DiscreteStateSpace dss;//=DiscreteStateSpace(A, B, C, D, Ts, method);

algorithm
    dss := DiscreteStateSpace(A, B, C, D,  Ts,  method);
    Ak := dss.A;
    Ck := dss.C;
    x_new := Ak*x + dss.B*u;
    y := Ck*x_new + dss.D*u2;

end ekf_LinearSystem;

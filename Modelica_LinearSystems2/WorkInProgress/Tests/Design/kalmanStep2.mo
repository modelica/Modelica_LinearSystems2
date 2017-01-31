within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function kalmanStep2
  "Example for one recursion of the conventional Kalman filter equations"

  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

  output Real P_out[:,:];
  output Real K[:,:];
  output Real R_out[:,:];

protected
  Real P[:,:]=identity(4);
  Real A[:,:]= [0.239960151286054, 0, 0.178712872351482, 0;
  -0.372217567033019, 1, 0.270264106474929, 0;
  -0.990087548834595, 0, 0.138859726355787, 0;
  -48.935406546735, 64.1, 2.39923411171398, 1];
  Real B[:,:]= [-1.23464449680517;
  -1.43828223420858;
  -4.48282453996389;
  -1.79989042995267];
  Real C[:,:] = [0, 1, 0, 0;
  0, 0, 0, 1;
  -128.2, 128.2, 0, 0];
  Real Q[:,:]=[1];

  Real R[:,:]=identity(3);

   DiscreteStateSpace dss=DiscreteStateSpace(A=A, B=B, C=C, D=zeros(3,1));
    Real x_init[:]={0,0,0,0};
   Real x_est[4];
   Real u[:]={0};
   Real y[:]={0,0,0};

  Integer info;

algorithm
  (K, P_out, R_out) := Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Internal.kfStepMatrices(A, B, C, P, Q, R);

//  Matrices.printMatrix(K,6,"K");
//  Matrices.printMatrix(P_out,6,"P");
  Matrices.printMatrix(R_out,6,"R");

  x_est := dss.A*x_init - K*(dss.C*x_init - y);
//  Modelica_LinearSystems2.Math.Vectors.printVector(x_est,6,"x_est");

 (x_est,K,P_out) :=DiscreteStateSpace.Internal.kfStepState(
    dss,
    P,
    Q,
    R,
    x_init,
    u,
    y);
// Modelica_LinearSystems2.Math.Vectors.printVector(x_est,6,"x_est");

//  Matrices.printMatrix(K,6,"K");
//  Matrices.printMatrix(P_out,6,"P");
  Matrices.printMatrix(R_out,6,"R");
end kalmanStep2;

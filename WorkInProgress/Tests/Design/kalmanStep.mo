within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function kalmanStep
  "Example for one recursion of the conventional Kalman filter equations"

  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

  output Real P_out[:,:];
  output Real K[:,:];
  output Real R_out[:,:];

protected
  Real P[:,:]=[0.5015,  0.4368,  0.2693,  0.6325;
               0.4368, 0.4818,  0.2639,  0.4148;
               0.2693,  0.2639,  0.1121,  0.6856;
               0.6325,  0.4148,  0.6856,  0.8906];
  Real A[:,:]=[0.2113,  0.8497,  0.7263,  0.8833;
               0.7560,  0.6857,  0.1985,  0.6525;
               0.0002,  0.8782,  0.5442,  0.3076;
               0.3303,  0.0683,  0.2320,  0.9329];
  Real B[:,:]=[0.0437,  0.7783,  0.5618;
               0.4818,  0.2119,  0.5896;
               0.2639,  0.1121,  0.6853;
               0.4148,  0.6856,  0.8906];
  Real Q[:,:]=[0.9329,  0.2146,  0.3126;
               0.2146,  0.2922,  0.5664;
               0.3126,  0.5664,  0.5935];
  Real C[:,:]=[0.3873,  0.9488,  0.3760,  0.0881;
               0.9222,  0.3435,  0.7340,  0.4498];
  Real R[:,:]=[1.0,  0.0;
               0.0,  1.0];

   DiscreteStateSpace dss=DiscreteStateSpace(A=A, B=B, C=C, D=zeros(2,3));
   Real x_init[:]={1,2,3,4};
   Real x_est[4];
   Real u[:]={1,2,3};
   Real y[:]={2,2};

  Integer info;

algorithm
  (K, P_out, R_out) := Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Internal.kfStepMatrices(A, B, C, P, Q, R);

  Matrices.printMatrix(K,6,"K");
  Matrices.printMatrix(P_out,6,"P");
  Matrices.printMatrix(R_out,6,"R");

  x_est := dss.A*x_init - K*(dss.C*x_init - y);
  Modelica_LinearSystems2.Math.Vectors.printVector(x_est,6,"x_est");

 (x_est,K,P_out) :=DiscreteStateSpace.Internal.kfStepState(
    dss,
    P,
    Q,
    R,
    x_init,
    u,
    y);
 Modelica_LinearSystems2.Math.Vectors.printVector(x_est,6,"x_est");

  Matrices.printMatrix(K,6,"K");
  Matrices.printMatrix(P_out,6,"P");
  Matrices.printMatrix(R_out,6,"R");
end kalmanStep;

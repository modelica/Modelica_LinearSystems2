within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function ekfSystem_pendular "Pendular state function"

  import Modelica.Math.Matrices;

  extends Modelica_LinearSystems2.DiscreteStateSpace.Internal.ekfSystemBase;

protected
  Real F[nx] "System state function";
  Real F_x[nx,nx] "Jacobian matrix of system";
  Real dFdx_1[nx];
  Real dFdx_2[nx];
  Real dFdx_3[nx]={0, 0, 0, 0};
  Real dFdx_4[nx]={0, 0, 1, 0};

  Real delta_x[nx] "Solution of tustin approximation";
  Real LU[nx,nx] "LU decomposition";
  Integer pivots[nx] "Pivots of LU decomposition";

algorithm
  F := {x[2], -(4000*cos(x[1])*sin(x[1])*x[2]^2 + 4905*sin(x[1]) + (u[1]
    *cos(x[1]))/10)/(1000*(4*sin(x[1]) + 1)), x[4], (u[1]/1000 - 10*x[2]^2 +
    (981*sin(2*x[1]))/50)/(4*sin(x[1]) + 1) + 10*x[2]^2};

  dFdx_1 := {0, -((981*cos(x[1]))/200 - u[1]/2500 + 4*x[2]^2*cos(2*x[1]) + 4*
    x[2]^2*sin(3*x[1]) - sin(x[1])*(12*x[2]^2 + u[1]/10000))/(4*sin(x[1]) +
    1)^2, 0, (40*cos(x[1])*x[2]^2 + (981*sin(3*x[1]))/25 - (1962*sin(x[1]))/
    25 - (u[1]*cos(x[1]))/250 + 8829/200)/(4*sin(x[1]) + 1)^2 - 981/200};
  dFdx_2 := {1, -(4*x[2]*sin(2*x[1]))/(4*sin(x[1]) + 1), 0, 20*x[2] - (20*x[2])/(4*sin(x[1]) + 1)};

  F_x := [dFdx_1,dFdx_2,dFdx_3,dFdx_4];

  (LU,pivots) := Matrices.LU(identity(nx) - (Ts/2)*F_x);

  delta_x := Matrices.LU_solve(LU, pivots, Ts*F);
  x_new := x + delta_x;
  y := {x_new[3] + 10*sin(x_new[1]),x_new[1]};
  Ak := Matrices.LU_solve2(LU, pivots, identity(nx) + (Ts/2)*(F_x));// Discretized Jacobi of F
  Ck := [10*cos(x_new[1]),0,1,0; 1,0,0,0];// Discretized Jacobi of y

end ekfSystem_pendular;

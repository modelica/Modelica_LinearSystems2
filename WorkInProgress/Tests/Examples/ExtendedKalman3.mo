within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model ExtendedKalman3 "Extended Kalman"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Types.Method;
  import Modelica.Math.Matrices.LU_solve2;
  import Modelica.Math.Matrices.solve2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  parameter Real xm_start[:]={0,0,0,0} "Start value of estimated state";
  parameter Real Q_vec[:] = {1e-5,1e-5,1e-5,1e-5} "Overall states confidence";
  parameter Real R_vec[:] = {20,2.13e-1} "Measurement confidence";
  parameter Real P_start[:,:] = diagonal(fill(0.5,nx));
  parameter Real Ts = 0.005 "Sample time of algorithm";
  parameter Integer nu = 1 "Number of system inputs";
  final parameter Integer nx = size(xm_start,1) "Number of system states";
  final parameter Integer ny = size(R_vec,1) "Number of observed measurements";
  final parameter Real Q[nx,nx] = diagonal(Q_vec) "Confidence of states";
  final parameter Real R[ny,ny] = diagonal(R_vec)
    "Confidence of measurements - large values low confidence | acts like coefficient of PT1-Filter";

  Real xp[nx] "A posteriori state estimation [phi,omega,x_1,vx_1]";
  Real xm[nx] "A priori state estimation [phi,omega,x_1,vx_1]";
//  Real F = u[1];
  Real Pp[nx,nx] "Init Covar-Matrix - large values = low confidence";
  Real Phi[nx,nx] "Transition matrix";
  Real F_x[nx,nx] "Jacobian matrix of system";
  Real dFdx_1[nx,1];
  Real dFdx_2[nx,1];
  Real dFdx_3[nx,1] = [0;0;0;0];
  Real dFdx_4[nx,1] = [0;0;1;0];
  Real H[ny,nx] "Jacobian matrix of system output";
  Real F_x_h[nx,1] "FOS system function RHS";
  Real Sol[nx,1] "Solution of tustin approximation";
  Real K[nx,ny] "Kalman gain";
  Real y_h[ny];
  Real P_det;

// Matrix operation stuff
  Real LU[nx,nx] "LU decomposition";
  Integer pivots[nx] "Pivots of LU decomposition";
  Modelica.Blocks.Interfaces.RealInput u[nu] "system inputs" 
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealOutput x[nx] = xm "corrected system states" 
    annotation (Placement(transformation(extent={{100,30},{120,50}})));
  Modelica.Blocks.Interfaces.RealInput y[ny] "system measurements" 
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));

  Modelica.Blocks.Interfaces.RealOutput y_out[:]=y_h "corrected system outputs"
    annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
initial equation
  Pp = P_start;
  xp = xm_start;
  dFdx_1 = [0; -((981*cos(xp[1]))/200 - u[1]/2500 + 4*xp[2]^2*cos(2*xp[1]) + 4*
    xp[2]^2*sin(3*xp[1]) - sin(xp[1])*(12*xp[2]^2 + u[1]/10000))/(4*sin(xp[1]) +
    1)^2; 0; (40*cos(xp[1])*xp[2]^2 + (981*sin(3*xp[1]))/25 - (1962*sin(xp[1]))/
    25 - (u[1]*cos(xp[1]))/250 + 8829/200)/(4*sin(xp[1]) + 1)^2 - 981/200];
  dFdx_2 = [1; -(4*xp[2]*sin(2*xp[1]))/(4*sin(xp[1]) + 1); 0; 20*xp[2] - (20*
    xp[2])/(4*sin(xp[1]) + 1)];
  F_x = [dFdx_1,dFdx_2,dFdx_3,dFdx_4];

  Phi = solve2(identity(nx) - (Ts/2)*F_x,identity(nx) + (Ts/2)*F_x);

equation
  when sample(0, Ts) then

    F_x_h = [xp[2]; -(4000*cos(xp[1])*sin(xp[1])*xp[2]^2 + 4905*sin(xp[1]) + (u[
      1]*cos(xp[1]))/10)/(1000*(4*sin(xp[1]) + 1)); xp[4]; (u[1]/1000 - 10*xp[2]
      ^2 + (981*sin(2*xp[1]))/50)/(4*sin(xp[1]) + 1) + 10*xp[2]^2];
    dFdx_1 = [0; -((981*cos(xp[1]))/200 - u[1]/2500 + 4*xp[2]^2*cos(2*xp[1]) + 4
      *xp[2]^2*sin(3*xp[1]) - sin(xp[1])*(12*xp[2]^2 + u[1]/10000))/(4*sin(xp[1])
       + 1)^2; 0; (40*cos(xp[1])*xp[2]^2 + (981*sin(3*xp[1]))/25 - (1962*sin(xp[
      1]))/25 - (u[1]*cos(xp[1]))/250 + 8829/200)/(4*sin(xp[1]) + 1)^2 - 981/200];
    dFdx_2 = [1; -(4*xp[2]*sin(2*xp[1]))/(4*sin(xp[1]) + 1); 0; 20*xp[2] - (20
      *xp[2])/(4*sin(xp[1]) + 1)];
    F_x = [dFdx_1,dFdx_2,dFdx_3,dFdx_4];
    H = [10*cos(xp[1]),0,1,0; 1,0,0,0];

  //Calculate transition matrix and a priori system state
  (LU,pivots) = Modelica.Math.Matrices.LU(identity(nx) - (Ts/2)*F_x);
  Phi = LU_solve2(LU, pivots, identity(nx) + (Ts/2)*(F_x));
  Sol = LU_solve2(LU, pivots, Ts*F_x_h);
  xm = pre(xp) + pre(Sol[:,1]);

  (K, Pp) = DiscreteStateSpace.Design.kfStepMatrices2(A=Phi, C=H, P=pre(Pp), Q=Q, R=R);

  y_h = {xm[3] + 10*sin(xm[1]), xm[1]};
  xp = xm + K*(y - y_h);

  P_det = Modelica.Math.Matrices.det(Pp);
end when;

  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),                      Text(
          extent={{-88,38},{90,-46}},
          lineColor={0,0,0},
          textString="EKF",
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}));
end ExtendedKalman3;

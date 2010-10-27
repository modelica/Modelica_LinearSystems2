within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model SR_EKF "Extended square root Kalman"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica.Math.Matrices.LU_solve2;
  import Modelica.Math.Matrices.solve2;
  import Modelica_LinearSystems2.Math.Matrices;

 parameter Real xm_start[:]={0,0,0,0} "Start value of estimated state";
  parameter Real P_start[:,:] = 0.5*identity(nx);
  final parameter Real S_start[:,:]=Matrices.cholesky(P_start,false);
  parameter Real Ts = 0.005 "Sample time of algorithm";
  parameter Integer nu = 1 "Number of system inputs";
  final parameter Integer nx = size(xm_start,1) "Number of system states";
  final parameter Integer ny = size(R,1) "Number of observed measurements";
  parameter Real Q[nx,nx] "Confidence of states";
  parameter Real R[:,:]
    "Confidence of measurements - large values low confidence | acts like coefficient of PT1-Filter";
  final parameter Real Cq[nx,nx] = Matrices.cholesky(Q,false)
    "Lower Cholesky factor of Q";
  final parameter Real Cr[ny,ny] = Matrices.cholesky(R,false)
    "Lower Cholesky factor of R";

  Real xp[nx] "A posteriori state estimation [phi,omega,x_1,vx_1]";
  Real xm[nx] "A priori state estimation [phi,omega,x_1,vx_1]";
//  Real F = u[1];
//  Real Pm[nx,nx] "A prori covaiance matrix";
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
  Real S[nx,nx] "Lower Cholesky factor of covariance matrix";
//  Real Crm[ny,ny] "modified Cr";
Real B[nx,1];

 Real Cr_new[ny,ny]
    "Modified Cholesky output or measurement noise covariance matrix";
//    output Real M[size(C, 1)+size(A, 1), size(C, 1)+2*size(A, 1)];
    Real M[ny+nx, nx+ny];

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
//  S = Matrices.cholesky(P_start,false);
  xp = xm_start;
   F_x_h = [xp[2]; -(4000*cos(xp[1])*sin(xp[1])*xp[2]^2 + 4905*sin(xp[1]) + (u[
      1]*cos(xp[1]))/10)/(1000*(4*sin(xp[1]) + 1)); xp[4]; (u[1]/1000 - 10*xp[2]
      ^2 + (981*sin(2*xp[1]))/50)/(4*sin(xp[1]) + 1) + 10*xp[2]^2];
  dFdx_1 = [0; -((981*cos(xp[1]))/200 - u[1]/2500 + 4*xp[2]^2*cos(2*xp[1]) + 4*
    xp[2]^2*sin(3*xp[1]) - sin(xp[1])*(12*xp[2]^2 + u[1]/10000))/(4*sin(xp[1]) +
    1)^2; 0; (40*cos(xp[1])*xp[2]^2 + (981*sin(3*xp[1]))/25 - (1962*sin(xp[1]))/
    25 - (u[1]*cos(xp[1]))/250 + 8829/200)/(4*sin(xp[1]) + 1)^2 - 981/200];
  dFdx_2 = [1; -(4*xp[2]*sin(2*xp[1]))/(4*sin(xp[1]) + 1); 0; 20*xp[2] - (20*
    xp[2])/(4*sin(xp[1]) + 1)];
  F_x = [dFdx_1,dFdx_2,dFdx_3,dFdx_4];

  Phi = solve2(identity(nx) - (Ts/2)*F_x,identity(nx) + (Ts/2)*F_x);
   H = [10*cos(xp[1]),0,1,0; 1,0,0,0];
  B= [0;
      -(cos(xp[1]))/10/(1000*(4*sin(xp[1]) + 1));
        0;
        (1/1000)/(4*sin(xp[1]) + 1)];
  // Sol = solve2(identity(nx) - (Ts/2)*F_x, Ts*F_x_h);
  // xm = xp + Sol[:,1];
//M2 = [Cr, (H)*(S), zeros(ny,nu); zeros(nx,ny), pre(Phi)*pre(S), 10*B];
  M = Matrices.LQ([Cr, H*S_start, zeros(ny,nx); zeros(nx,ny),Phi*S_start, Cq]);
  S = Matrices.triangle(M[ny+1:ny+nx,ny+1:nx+ny],false);
  Cr_new =  Matrices.triangle(M[1:ny,1:ny],false);
   K =  Matrices.Internal.solve2rSym(Cr_new,M[ny+1:ny+nx,1:ny],true,false);
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
    B = [0;
        -(cos(xp[1]))/10/(1000*(4*sin(xp[1]) + 1));
         0;
         (1/1000)/(4*sin(xp[1]) + 1)];
  //Calculate transition matrix and a priori system state
  (LU,pivots) = Modelica.Math.Matrices.LU(identity(nx) - (Ts/2)*F_x);
  Phi = LU_solve2(LU, pivots, identity(nx) + (Ts/2)*F_x);
  Sol = LU_solve2(LU, pivots, Ts*F_x_h);
  xm = pre(xp) + pre(Sol[:,1]);
  y_h = {xm[3] + 10*sin(xm[1]), xm[1]};

//  (K, S, Cr_new, M) =  Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design.sr_kfStepMatrices2(A=pre(Phi), C=pre(H), S=pre(S), Cq=B, Cr=Cr);

// M = Matrices.LAPACK.dgelqf([Cr, pre(H)*pre(S), zeros(ny,nx); zeros(nx,ny), pre(Phi)*pre(S), Cq]);
  M = Matrices.LQ([Cr, pre(H)*pre(S), zeros(ny,nx); zeros(nx,ny), pre(Phi)*pre(S), Cq]);
  S = Matrices.triangle(M[ny+1:ny+nx,ny+1:nx+ny],false);
  Cr_new =  Matrices.triangle(M[1:ny,1:ny],false);
   K =  Matrices.Internal.solve2rSym(Cr_new,M[ny+1:ny+nx,1:ny],true,false);

  xp = xm + K*(y - y_h);

end when;

  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),                      Text(
          extent={{-84,38},{94,-46}},
          lineColor={0,0,0},
          textString="EKF",
          fillPattern=FillPattern.Solid,
          fillColor={0,0,0}),
        Line(
          points={{-96,68},{-88,68},{-82,-70},{-70,68},{86,68},{86,54}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5)}),
    experiment(StopTime=10),
    experimentSetupOutput);
end SR_EKF;

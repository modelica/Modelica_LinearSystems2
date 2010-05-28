within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model SR_UnscentedKalman2 "Unscented Kalman filter - SR computation"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Types.Method;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Math.Matrices.LU_solve2;
  import Modelica.Math.Matrices.solve2;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

 parameter Real xm_start[:]={0,0,0,0} "Start value of estimated state";
  parameter Real Q_vec[:] = {1e-5,1e-5,1e-5,1e-5} "Overall states confidence";
  parameter Real R_vec[:] = {0.1,2.13e-1} "Measurement confidence";
  parameter Real P_start[:,:] = diagonal(fill(0.5,nx));
  parameter Real Ts = 0.005 "Sample time of algorithm";
  parameter Integer nu = 1 "Number of system inputs";
  final parameter Integer nx = size(xm_start,1) "Number of system states";
  final parameter Integer ny = size(R_vec,1) "Number of observed measurements";
  final parameter Real Q[nx,nx] = diagonal(Q_vec) "Confidence of states";
  final parameter Real R[ny,ny] = diagonal(R_vec)
    "Confidence of measurements - large values low confidence | acts like coefficient of PT1-Filter";

  final parameter Real CfQ[nx,nx] =  Matrices.cholesky(Q,false)
    "Confidence of states";
  final parameter Real CfR[ny,ny] =  Matrices.cholesky(R,false)
    "Confidence of measurements - large values low confidence | acts like coefficient of PT1-Filter";

  parameter Real alpha=0.1;
  parameter Real beta=2;
  parameter Real kappa=0 "Scaling parameter";

  Real xp[nx] "A posteriori state estimation [phi,omega,x_1,vx_1]";
  Real Sp[nx,nx] "Init Cholesky Covar-Matrix - large values = low confidence";

  Real F_x[nx,nx] "Jacobian matrix of system";
  Real dFdx_1[nx,1];
  Real dFdx_2[nx,1];
  Real dFdx_3[nx,1]=[0; 0; 0; 0];
  Real dFdx_4[nx,1]=[0; 0; 1; 0];
//  Real H[ny,nx] "jacobian matrix of system output";
  Real F_x_h[nx,1] "FOS system function RHS";
  Real Sol[nx] "Solution of tustin approximation";
  Real K[nx,ny] "Kalman gain";
//  Real y_h[ny];
  Real P_det;

  Real mux[nx] "Predicted mean";
  Real muy[ny] "Predicted mean";
  Real Skx[nx,nx] "Cholesky factor of the transformed covariance matrix";
  Real Syy[ny,size(Syy, 1)] "Transformed covariance matrix";
  Real Rxy[nx,ny] "Transformed cross covariance matrix";

  Modelica.Blocks.Interfaces.RealInput u[nu] "system inputs" 
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealOutput x[nx]=vector(mux)
    "corrected system states" 
    annotation (Placement(transformation(extent={{100,30},{120,50}})));
  Modelica.Blocks.Interfaces.RealInput y[ny,1] "system measurements" 
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));

  Modelica.Blocks.Interfaces.RealOutput y_out[ny] "corrected system outputs" 
    annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
initial equation
  xp = xm_start;
  Sp = Matrices.cholesky(P_start);

    dFdx_1 = [0; -((981*cos(xp[1]))/200 - u[1]/2500 + 4*xp[2]^2*cos(2*xp[1]) + 4
      *xp[2]^2*sin(3*xp[1]) - sin(xp[1])*(12*xp[2]^2 + u[1]/10000))/(4*sin(xp[1])
       + 1)^2; 0; (40*cos(xp[1])*xp[2]^2 + (981*sin(3*xp[1]))/25 - (1962*sin(xp[
      1]))/25 - (u[1]*cos(xp[1]))/250 + 8829/200)/(4*sin(xp[1]) + 1)^2 - 981/200];
    dFdx_2 = [1; -(4*xp[2]*sin(2*xp[1]))/(4*sin(xp[1]) + 1); 0; 20*xp[2] - (20*
      xp[2])/(4*sin(xp[1]) + 1)];
    F_x = [dFdx_1,dFdx_2,dFdx_3,dFdx_4];
equation
  when sample(0, Ts) then
    F_x_h = [xp[2]; -(4000*cos(xp[1])*sin(xp[1])*xp[2]^2 + 4905*sin(xp[1]) + (u[
      1]*cos(xp[1]))/10)/(1000*(4*sin(xp[1]) + 1)); xp[4]; (u[1]/1000 - 10*xp[2]
      ^2 + (981*sin(2*xp[1]))/50)/(4*sin(xp[1]) + 1) + 10*xp[2]^2];
    dFdx_1 = [0; -((981*cos(xp[1]))/200 - u[1]/2500 + 4*xp[2]^2*cos(2*xp[1]) + 4
      *xp[2]^2*sin(3*xp[1]) - sin(xp[1])*(12*xp[2]^2 + u[1]/10000))/(4*sin(xp[1])
       + 1)^2; 0; (40*cos(xp[1])*xp[2]^2 + (981*sin(3*xp[1]))/25 - (1962*sin(xp[
      1]))/25 - (u[1]*cos(xp[1]))/250 + 8829/200)/(4*sin(xp[1]) + 1)^2 - 981/200];
    dFdx_2 = [1; -(4*xp[2]*sin(2*xp[1]))/(4*sin(xp[1]) + 1); 0; 20*xp[2] - (20*
      xp[2])/(4*sin(xp[1]) + 1)];
    F_x = [dFdx_1,dFdx_2,dFdx_3,dFdx_4];
//    H = [10*sin(xp[1]),0,1,0; 1,0,0,0];

    Sol = Modelica.Math.Matrices.solve(identity(nx) - (Ts/2)*F_x,Ts*vector(F_x_h));

//    (mux,Skx) = DiscreteStateSpace.Design.sr_ukfPredict(pre(xp), transpose(pre(Sp)), CfQ, u, pre(Sol), alpha, beta, kappa);
//    (mux,Skx) = DiscreteStateSpace.Design.sr_ukfPredict2(pre(xp), transpose(pre(Sp)), CfQ, u, pre(F_x), alpha, beta, kappa,Ts);
    (mux,Skx) = DiscreteStateSpace.Design.sr_ukfPredict3(pre(xp), transpose(pre(Sp)), CfQ, u, alpha, beta, kappa, Ts);

    (muy,Syy,Rxy) = DiscreteStateSpace.Design.sr_ukfUpdate(mux, transpose(Skx), CfR, u, alpha, beta, kappa);
    (K,Sp,xp) = DiscreteStateSpace.Design.sr_ukfEstimate(vector(y), mux, muy, Skx, Syy, Rxy);

    y_out = {mux[3] + 10*sin(mux[1]),mux[1]};

    P_det = Modelica_LinearSystems2.Math.Matrices.det(transpose(Sp)*Sp);
  end when;

  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-88,24},{110,-62}},
          lineColor={0,0,0},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="UKF"),
        Line(
          points={{-94,58},{-86,58},{-80,-80},{-68,58},{88,58},{88,44}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5)}));
end SR_UnscentedKalman2;

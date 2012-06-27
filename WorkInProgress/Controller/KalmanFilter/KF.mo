within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model KF "Discrete State Space block"

  import Modelica;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;
  import Modelica_LinearSystems2.WorkInProgress.MPC;
  import Modelica_LinearSystems2.WorkInProgress.Optimizer;
  import Modelica_LinearSystems2.StateSpace;

  parameter StateSpace ss "Continuous linear system model";
  parameter Modelica.SIunits.Time sampleTime=0.5
    "Base sample time for discrete blocks";
  parameter Real Q[size(dss.B, 2),size(dss.B, 2)]=identity(size(dss.B, 2))
    "Input or process noise covariance matrix of the previous instant";
  parameter Real R[size(dss.C, 1),size(dss.C, 1)]=identity(size(dss.C, 1))
    "Output or measurement noise covariance matrix of the previous instant";
  parameter Real P0[size(dss.A, 1),size(dss.A, 1)]=10*identity(size(dss.A, 1))
    "Initial state covariance matrix of the previous instant";
  parameter Real x_init[size(dss.A,1)]=fill(0,size(dss.A,1));

  parameter DiscreteStateSpace dss=DiscreteStateSpace(ss,sampleTime,method=Modelica_LinearSystems2.Types.Method.StepExact)
    "Discrete linear system model";

public
  Real P[size(dss.A, 1),size(dss.A, 1)](start=P0)
    "State covariance matrix of the previous instant";
  Real P_pre[size(dss.A, 1),size(dss.A, 1)](start=P0)
    "State covariance matrix of the previous instant";
  Real K[size(dss.A, 1),size(dss.C, 1)] "Kalman filter gain matrix";

  discrete Modelica.Blocks.Interfaces.RealOutput x_est[size(dss.A, 1)]
    "Current estmated state vector"
    annotation (extent=[100, -10; 120, 10], Placement(transformation(extent={{100,
            -10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput x_est_pre[size(dss.A, 1)](start=x_init)
    "Preceding estimated state vector"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,40},
            {-100,80}})));
  Modelica.Blocks.Interfaces.RealInput y[size(dss.C, 1)] "measured output"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,
            -80},{-100,-40}}), iconTransformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealInput u[size(dss.B, 2)] "Current input"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,
            -20},{-100,20}})));

  Real z[size(dss.A,1)];
// Real M[size(dss.C, 1),size(dss.C, 1)]
//     "State covariance matrix of the previous instant";
//   Real PCT[size(dss.A, 1),size(dss.C, 1)] "Kalman filter gain matrix";
initial equation
//  (x_est,K,P) = DiscreteStateSpace.Internal.kfStepState(dss, P0, Q, R, x_init, u, y);

equation
//   when initial() then
//     (x_est, K, P) = DiscreteStateSpace.Internal.kfStepState(dss, P0,Q,R, x_init,u,y);
//   end when;

  when sample(sampleTime, sampleTime) then
   P_pre=pre(P);
   z =  dss.A*x_est_pre +dss.B*u;
//        (x_est, K, P) = DiscreteStateSpace.Internal.kfStepState(dss, P_pre,Q,R, x_est_pre,u,y);
        (K,P) =  DiscreteStateSpace.Internal.kfStepMatrices(dss.A, dss.B, dss.C, P_pre, Q, R);
//          x_est =  dss.A*x_est_pre - dss.A*K*(dss.C*x_est_pre - y);
//          x_est =  dss.A*x_est_pre +dss.B*u - K*(dss.C*x_est_pre - y);
          x_est =  dss.A*x_est_pre +dss.B*u - K*(dss.C*z - y);
  end when;
  annotation (Documentation(info="<html>

</HTML>
"), Diagram(graphics),
    Icon(graphics={Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
          Text(
          extent={{-60,30},{60,-30}},
          lineColor={0,0,0},
          textString="K F"),
        Text(
          extent={{-140,100},{-100,80}},
          lineColor={0,0,255},
          textString="x_est_pre"),
        Text(
          extent={{-140,-20},{-100,-40}},
          lineColor={0,0,255},
          textString="y"),
        Text(
          extent={{-140,40},{-100,20}},
          lineColor={0,0,255},
          textString="u")}));
end KF;

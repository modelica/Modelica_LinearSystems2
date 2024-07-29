within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model KF_inner "Discrete State Space block"

  import Modelica;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.WorkInProgress.Optimizer;

  parameter DiscreteStateSpace dss=DiscreteStateSpace(A=[1],B=[1],C=[1],D=[1])
    "Discrete linear system model";
  parameter Modelica.Units.SI.Time sampleTime=0.5
    "Base sample time for discrete blocks";
  parameter Real wB[size(dss.B, 1),size(dss.B, 2)]=ones(size(dss.B, 1),size(dss.B, 2))
    "Wheighting matrix for input noise covariance matrix of the previous instant";
  parameter Real Q[size(dss.B, 2),size(dss.B, 2)]=identity(size(dss.B, 2))
    "Input or process noise covariance matrix of the previous instant";
  parameter Real R[size(dss.C, 1),size(dss.C, 1)]=identity(size(dss.C, 1))
    "Output or measurement noise covariance matrix of the previous instant";
  parameter Real P0[size(dss.A, 1),size(dss.A, 1)]=10*identity(size(dss.A, 1))
    "Initial state covariance matrix of the previous instant";
  parameter Real x_init[size(dss.A,1)]=fill(0,size(dss.A,1));

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
  Modelica.Blocks.Interfaces.RealInput y[size(dss.C, 1)] "measured output"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,
            -80},{-100,-40}}), iconTransformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealInput u[size(dss.B, 2)] "Current input"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,40},
            {-100,80}})));

  Real z[size(dss.A,1)];

equation
  when sample(0, sampleTime) then
   P_pre=pre(P);
   z =  dss.A*pre(x_est) +dss.B*pre(u);
//        (x_est, K, P) = DiscreteStateSpace.Internal.kfStepState(dss, P_pre,Q,R, x_est_pre,u,y);
        (K,P) =  Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design.kfStepMatrices(
                                                          dss.A, wB, dss.C, P_pre, Q, R);
//          x_est =  dss.A*x_est_pre - dss.A*K*(dss.C*x_est_pre - y);
//          x_est =  dss.A*x_est_pre +dss.B*u - K*(dss.C*x_est_pre - y);
//          x_est =  dss.A*pre(x_est) +dss.B*pre(u) - K*(dss.C*z + dss.D*pre(u) - y);
          x_est =  z - K*(dss.C*z + dss.D*pre(u) - y);
  end when;
  annotation (Documentation(info="<html>
</html>"),    Icon(graphics={Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-60,30},{60,-30}},
          textColor={0,0,0},
          textString="K F"),
        Text(
          extent={{-140,-20},{-100,-40}},
          textColor={0,0,255},
          textString="y"),
        Text(
          extent={{-140,40},{-100,20}},
          textColor={0,0,255},
          textString="u")}));
end KF_inner;

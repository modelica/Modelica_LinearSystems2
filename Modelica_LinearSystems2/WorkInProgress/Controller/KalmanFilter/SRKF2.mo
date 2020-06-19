within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model SRKF2 "Discrete State Space block"
  import Modelica_LinearSystems2;

  import Modelica;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.WorkInProgress.MPC;
  import Modelica_LinearSystems2.WorkInProgress.Optimizer;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.WorkInProgress.Controller;

  parameter StateSpace ss "Continuous linear system model";
  parameter Modelica.Units.SI.Time sampleTime=0.5
    "Base sample time for discrete blocks";
  parameter DiscreteStateSpace dss=DiscreteStateSpace(ss,sampleTime,method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact)
    "Discrete linear system model";
  parameter Real Q[size(dss.A, 1),size(dss.A, 1)]=identity(size(dss.A, 1))
    "Input or process noise covariance matrix of the previous instant";

  parameter Real R[size(dss.C, 1),size(dss.C, 1)]=identity(size(dss.C, 1))
    "Output or measurement noise covariance matrix of the previous instant";
  parameter Real P0[size(dss.A, 1),size(dss.A, 1)]=10*identity(size(dss.A, 1))
    "Initial state covariance matrix of the previous instant";
  parameter Real x_init[size(dss.A, 1)]=fill(0, size(dss.A, 1));

public
  Modelica.Blocks.Interfaces.RealOutput x_est[size(dss.A, 1)]
    "Current estmated state vector"
    annotation (extent=[100, -10; 120, 10], Placement(transformation(extent={{100,
            -10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput y[size(dss.C, 1)] "measured output"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,
            -80},{-100,-40}}), iconTransformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealInput u[size(dss.B, 2)] "Current input"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(extent={{-140,40},
            {-100,80}})));

  Controller.KalmanFilter.SRKF_inner2 kf_int(
    sampleTime=sampleTime,
    dss=dss,
    Q=Q,
    R=R,
    P0=P0,
    x_init=x_init) annotation (Placement(transformation(extent={{-40,-30},{20,30}})));
equation

  connect(kf_int.x_est, x_est) annotation (Line(
      points={{23,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(kf_int.y, y) annotation (Line(
      points={{-46,-18},{-80,-18},{-80,-60},{-120,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, kf_int.u) annotation (Line(
      points={{-120,60},{-100,60},{-100,62},{-80,62},{-80,18},{-46,18}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Documentation(info="<html>

</html>"),    Icon(graphics={Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
          Text(
          extent={{-60,30},{60,-30}},
          lineColor={0,0,0},
          textString="K F"),
        Text(
          extent={{-140,-20},{-100,-40}},
          lineColor={0,0,255},
          textString="y"),
        Text(
          extent={{-140,92},{-100,72}},
          lineColor={0,0,255},
          textString="u")}));
end SRKF2;

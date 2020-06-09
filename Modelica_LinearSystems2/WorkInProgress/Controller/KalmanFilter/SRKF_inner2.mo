within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model SRKF_inner2 "Discrete State Space block"

  import Modelica;
  import MSLcholesky = Modelica.Math.Matrices.cholesky;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Math.Matrices.triangle;

  parameter DiscreteStateSpace dss=DiscreteStateSpace(
      A=[1],
      B=[1],
      C=[1],
      D=[1]) "Discrete linear system model";
  parameter Modelica.Units.SI.Time sampleTime=0.5
    "Base sample time for discrete blocks";
  parameter Real Q[size(dss.A, 1),size(dss.A, 1)]=identity(size(dss.A, 1))
    "Input or process noise covariance matrix of the previous instant";
  final parameter Real Cq[size(dss.A, 1),size(dss.A, 1)]=MSLcholesky(Q,false)
    "Input or process noise covariance matrix of the previous instant";
  parameter Real R[size(dss.C, 1),size(dss.C, 1)]=identity(size(dss.C, 1))
    "Output or measurement noise covariance matrix of the previous instant";
  final parameter Real Cr[size(dss.C, 1),size(dss.C, 1)]=MSLcholesky(R,false)
    "Output or measurement noise covariance matrix of the previous instant";
  parameter Real P0[size(dss.A, 1),size(dss.A, 1)]=10*identity(size(dss.A, 1))
    "Initial state covariance matrix of the previous instant";
  parameter Real x_init[size(dss.A, 1)]=fill(0, size(dss.A, 1));

public
  Real S[size(dss.A, 1),size(dss.A, 1)](start=MSLcholesky(P0,false))
    "State covariance matrix of the previous instant";
  Real K[size(dss.A, 1),size(dss.C, 1)] "Kalman filter gain matrix";
//  Real K1[size(dss.A, 1),size(dss.C, 1)] "Kalman filter gain matrix";

  Real Cr_out[size(dss.C, 1),size(dss.C, 1)];
  Real M[size(dss.C, 1)+size(dss.A, 1), size(dss.C, 1)+2*size(dss.A, 1)];

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

  Real z[size(dss.A, 1)];

initial equation
  S=MSLcholesky(P0,false);
  x_est=x_init;

equation
 when sample(0, sampleTime) then
    z = dss.A*pre(x_est) + dss.B*pre(u);
    (K,S, Cr_out,M) = Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design.sr_kfStepMatrices2(
                                                                   dss.A, dss.C, pre(S), Cq, Cr);
//    K=Matrices.Internal.solve2rSym(Cr_out,K1,true,false);
    x_est = z - K*(dss.C*pre(x_est) + dss.D*pre(u) - y);
  end when;
  annotation (
    Documentation(info="<html>

</html>"),    Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
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
          extent={{-140,40},{-100,20}},
          lineColor={0,0,255},
          textString="u")}));
end SRKF_inner2;

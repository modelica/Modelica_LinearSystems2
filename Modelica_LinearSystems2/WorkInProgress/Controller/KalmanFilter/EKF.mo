within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model EKF "Extended Kalman filter"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

 extends Modelica_LinearSystems2.Controller.Interfaces.PartialDiscreteBlock(final
      initType = Modelica_LinearSystems2.Controller.Types.Init.InitialState);

  replaceable function ekfFunction =
            Modelica_LinearSystems2.DiscreteStateSpace.Internal.ekfSystemDummy
            constrainedby Modelica_LinearSystems2.DiscreteStateSpace.Internal.ekfSystemBase(
                                                                      redeclare input Real x[4])
    "Function to calculate xk, yk, Ak, Ck"
    annotation (
      choicesAllMatching,
            Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-10-25</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));

  parameter Real x_est_init[:] "Initial value for state estimation";
  parameter Real Q[:,size(Q,1)]=identity(size(x_est_init,1))
    "Covariance matrix of the associated process noise";
  parameter Real G[size(x_est_init,1),size(Q,1)]
    "Process noise weight matrix (usually identity or input matrix)";
  final parameter Real Q2[size(x_est_init,1),size(x_est_init,1)]=G*Q*transpose(G)
    "Appropriately weighted process noise covariance matrix";
  parameter Real R[:,size(R,1)] "Covariance matrix of the measurement noise";
  parameter Real M_init[:,:] = identity(size(x_est_init,1))
    "Initial value of the Riccati matrix M";
  parameter Integer nu = 1 "Number of system inputs";
  final parameter Integer nx = size(x_est_init,1) "Number of system states";
  final parameter Integer ny = size(R,1) "Number of observed measurements";

  Real M[nx,nx] "Solution of the discrete Riccati equation";
  Real K[nx,ny] "Kalman filter gain matrix";

  Modelica.Blocks.Interfaces.RealInput u[nu] "System input vector"
    annotation (Placement(transformation(extent={{-140,30},{-100,70}})));
  Modelica.Blocks.Interfaces.RealOutput x_est[nx] "Estimated state vector"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput y_measure[ny] "Measured outputs"
    annotation (Placement(transformation(extent={{-140,-70},{-100,-30}})));

  Modelica.Blocks.Interfaces.RealOutput y_est[ny] "Estimated system outputs"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-110})));

protected
  outer Modelica_LinearSystems2.Controller.SampleClock sampleClock
    "Global options";
initial equation
  x_est = x_est_init;
  M = M_init;

equation
  when {sampleTrigger} then
    (x_est, y_est, M, K) = Design.EKF(function ekfFunction(), pre(x_est), pre(u), u,
    y_measure, pre(M), Q2, R, Ts);
  end when;

  annotation ( Icon(graphics={
        Text(
          extent={{-84,34},{76,-40}},
          lineColor={0,0,0},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="EKF"),
        Text(
          extent={{52,30},{114,-26}},
          lineColor={0,0,255},
          textString="x"),
        Line(
          points={{76,18},{84,28},{92,18}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-2,-56},{6,-46},{14,-56}},
          color={0,0,255},
          smooth=Smooth.None),
        Text(
          extent={{-24,-42},{38,-98}},
          lineColor={0,0,255},
          textString="y"),
        Text(
          extent={{-118,82},{-56,26}},
          lineColor={0,0,255},
          textString="u"),
        Text(
          extent={{-118,-18},{-56,-74}},
          lineColor={0,0,255},
          textString="y"),
        Text(
          extent={{-106,-54},{-20,-84}},
          lineColor={0,0,255},
          textString="m")}),
    Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-10-25</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end EKF;

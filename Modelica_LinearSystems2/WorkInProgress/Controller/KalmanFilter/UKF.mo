within Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter;
model UKF "Unscented Kalman filter"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design;

  extends Modelica_LinearSystems2.Controllers.Interfaces.PartialDiscreteBlock(
    final initType = Modelica_LinearSystems2.Controllers.Types.Init.InitialState);

  replaceable function F_function =
    Modelica_LinearSystems2.DiscreteStateSpace.Internal.fSigmaDummy
    constrainedby Modelica_LinearSystems2.DiscreteStateSpace.Internal.fBase(
      redeclare input Real x[4])
    "Function F() in x_k+1 = F(x_k, u_k, Ts)" annotation(choicesAllMatching);
  replaceable function H_function =
    Modelica_LinearSystems2.DiscreteStateSpace.Internal.hSigmaDummy
    constrainedby Modelica_LinearSystems2.DiscreteStateSpace.Internal.hBase(
      redeclare input Real x[4])
    "Function H() in y_k = H(x_k, u_k, Ts)" annotation(choicesAllMatching);

  parameter Real x_est_init[:] "Initial value for state estimation";
  parameter Real Q[:,size(Q,1)]=identity(size(x_est_init,1))
    "Covariance matrix of the associated process noise";
  parameter Real G[size(x_est_init,1),size(Q,1)]
    "Process noise weight matrix (usually identity or input matrix)";
  final parameter Real Q2[size(x_est_init,1),size(x_est_init,1)]=G*Q*transpose(G)
    "Appropriately weighted process noise covariance matrix";
  parameter Real R[:,size(R,1)] "Covariance matrix of the measurement noise";
  parameter Real P_init[:,:] = identity(size(x_est_init,1))
    "Initial value of the error covariance matrix";
  parameter Integer nu = 1 "Number of system inputs";
  parameter Real alpha=0.1 "Spread of sigma points";
  parameter Real beta=2 "Characteristic of the distribution of x";
  parameter Real kappa=0 "Scaling kurtosis of sigma point distribution";
  final parameter Integer nx = size(x_est_init,1) "Number of system states";
  final parameter Integer ny = size(R,1) "Number of observed measurements";

  Real P[nx,nx] "Error covariance matrix";
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
  outer Modelica_LinearSystems2.Controllers.SampleClock sampleClock
    "Global options";
initial equation
  x_est = x_est_init;
  P = P_init;
equation
  when {sampleTrigger} then
    (x_est, y_est, P, K) = Design.UKF(function F_function(), function H_function(), pre(x_est), pre(u), y_measure, pre(P), Q2, R, alpha, beta, kappa, Ts);
  end when;

  annotation ( Icon(graphics={
        Text(
          extent={{-84,34},{76,-40}},
          textColor={0,0,0},
          textString="UKF"),
        Text(
          extent={{52,30},{114,-26}},
          textColor={0,0,255},
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
          textColor={0,0,255},
          textString="y"),
        Text(
          extent={{-118,82},{-56,26}},
          textColor={0,0,255},
          textString="u"),
        Text(
          extent={{-118,-18},{-56,-74}},
          textColor={0,0,255},
          textString="y"),
        Text(
          extent={{-106,-54},{-20,-84}},
          textColor={0,0,255},
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
end UKF;

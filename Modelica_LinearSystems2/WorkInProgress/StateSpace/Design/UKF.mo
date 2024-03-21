within Modelica_LinearSystems2.WorkInProgress.StateSpace.Design;
function UKF "Unscented Kalman filter design function"

  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input Real xpre[:] "State at instant k-1";
  input Real upre[:] "Input at instant k-1";
  input Real y[:] "Output at instant k";
  input Real Ppre[size(xpre,1),size(xpre,1)]
    "Error covariance matrix at instant k-1";
  input Real Q[size(xpre,1),size(xpre,1)] = identity(size(xpre,1))
    "Weighted covariance matrix of the associated process noise (F*Q*F')";
  input Real R[size(y,1),size(y,1)] = identity(size(y,1))
    "Covariance matrix of the measurement noise";
  input Real alpha=0.1 "Spread of sigma points";
  input Real beta=2 "Characteristic of the distribution of x";
  input Real kappa=0 "Kurtosis scaling of sigma point distribution";
  input Modelica.Units.SI.Time Ts "Sample time";

  output Real x_est[size(xpre,1)] "Estimated state vector";
  output Real y_est[size(y,1)] "Estimated output";
  output Real P[size(Ppre,1),size(Ppre,1)] "Error covariance matrix";
  output Real K[size(xpre,1),size(y,1)] "Kalman filter gain matrix";

  replaceable function predict=DiscreteStateSpace.Internal.ukfPredict(redeclare
        Modelica_LinearSystems2.DiscreteStateSpace.Internal.fSigmaDummy
          fSigma) annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-06-11</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));

  replaceable function update=DiscreteStateSpace.Internal.ukfUpdate(redeclare
        Modelica_LinearSystems2.DiscreteStateSpace.Internal.hSigmaDummy
          hSigma) annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-06-11</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));

  replaceable function estimate=DiscreteStateSpace.Internal.ukfEstimate(redeclare
        Modelica_LinearSystems2.DiscreteStateSpace.Internal.hSigmaDummy
          yOut) annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-06-11</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));

protected
  Integer nx = size(xpre,1) "Number of system states";
  Integer ny = size(y,1) "Number of observed measurements";

  Real mux[nx] "Predicted state mean";
  Real muy[ny] "Predicted output mean";
  Real Rxx[nx,nx] "Transformed covariance matrix";
  Real Ryy[ny,ny] "Transformed covariance matrix";
  Real Rxy[nx,ny] "Transformed cross covariance matrix";

algorithm
    (mux,Rxx) := predict(xpre, upre, Ppre, Q, alpha, beta, kappa, Ts);
    (muy,Ryy,Rxy) := update(mux, upre, Rxx, R,  alpha, beta, kappa, Ts);
    (K,P,x_est, y_est) := estimate(y, mux, muy, upre, Rxx, Ryy, Rxy, Ts);
    annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-06-11</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>", info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(x_est, y_est, P, K) = DiscreteStateSpace.Design.<strong>UKF</strong>(x_pre, u_pre, y, P_pre, Q, R, alpha, beta, kappa, Ts)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <strong>UKF</strong> computes one recursion of the Unscented Kalman filter. Unscented Kalman filters are similar to Extended Kalman filters
but using statistical linearization where extended Kalman filter apply the user-provided derivation of the system equation. Instead of explicit derivation
linear regression between spcifically chosen sample points (sigma points). See [1] for more information.
</p>
<p>
See also <a href=\"modelica://Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design.UKF_SR\">UKF_SR</a>, where the square root method to deal with positive definte matrices is applied to solve the mathematically identical problem.
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1]</dt>
<dd> <strong>http://en.wikipedia.org/wiki/Kalman_filter#Unscented_Kalman_filter</strong>.
    <br>&nbsp;</dd>
</dl>

</html>"));
end UKF;

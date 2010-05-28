within Modelica_LinearSystems2.WorkInProgress.Tests.DiscreteSystems;
function testUKF
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

protected
  Real x[:]={1,2,3} "Estimated vector of previous instant";
  Real P[size(x, 1),size(x, 1)]= [1,3,0;5,-1,2;3,1,0]*transpose([1,3,0;5,-1,2;3,1,0])
    "State covariance matrix of the previous instant";
  Real S[size(P, 1),size(P, 1)] = Matrices.cholesky(P);
  Real Q[size(P, 1),:]= P/3 "Covariance matrix of the process noise";
  Real CfQ[size(Q, 1),size(Q, 1)] = Matrices.cholesky(Q);
  Real R[:,:]= [2,3;4,3]*transpose([2,3;4,3])
    "Covariance matrix of the process noise";
  Real CfR[size(R, 1),size(R, 1)] = Matrices.cholesky(R);
  Real uk[:]={1} "Input at instant k";
  Real alpha=0.5;
  Real beta=2;
  Real kappa=0.1 "Scaling parameter";
  Real y[2]={1,2};

public
  output Real mu[size(x, 1)] "Predicted mean";
  output Real Pp[size(x, 1),size(x, 1)] "Transformed covariance matrix";
  output Real ym[size(R, 1)] "Predicted mean";
  output Real Ryy[size(R, 1),size(R, 2)] "Transformed covariance matrix";
  output Real Rxy[size(x, 1),size(R, 2)] "Transformed covariance matrix";
  output Real K[size(x,1),size(y,1)];
  output Real Pu[size(x,1),size(x,1)];
  output Real xmu[size(x,1)];

algorithm
  (mu, Pp) := DiscreteStateSpace.Design.ukfPredict(x, P,Q, uk, alpha, beta, kappa);
  (ym, Ryy, Rxy) := DiscreteStateSpace.Design.ukfUpdate(mu, Pp,R, uk, alpha, beta, kappa);
  (K, Pu, xmu) := Modelica_LinearSystems2.DiscreteStateSpace.Design.ukfEstimate(
                                                         y, mu,ym, Pp, Ryy, Rxy);

end testUKF;

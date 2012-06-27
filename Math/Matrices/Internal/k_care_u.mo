within Modelica_LinearSystems2.Math.Matrices.Internal;
function k_care_u
  "Calculate the upper bound of the CARE, i.e. Q + A'*X + X*A - X*G*X = 0  condition number using Lyapunov equations"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,size(A, 1)] "care-matrix A";
  input Real Q[:,size(Q, 1)] "care-matrix Q";
  input Real G[size(A, 1),size(A, 2)] "care-matrix G";
  input Real X[size(A, 1),size(A, 2)] "solution of care";

  output Real ku "upper bound of the care condition number";

protected
  Real Z0[size(A, 1),size(A, 2)]=Matrices.lyapunov(A - G*X, -identity(size(A, 1)))
    "solution of lyapunov equation H'*Z0+Z0*H=-I, H=A-G*X";
  Real Z2[size(A, 1),size(A, 2)]=Matrices.lyapunov(A - G*X, -X*X)
    "solution of lyapunov equation H'*Z2+Z2*H=-X*X, H=A-G*X";
  Real normA=Matrices.norm(A, 2) "spectral norm of matrix A (largest sv)";
  Real normQ=Matrices.norm(Q, 2) "spectral norm of matrix Q (largest sv)";
  Real normG=Matrices.norm(G, 2) "spectral norm of matrix G (largest sv)";
  Real normX=Matrices.norm(X, 2) "spectral norm of matrix X (largest sv)";
  Real normZ0=Matrices.norm(Z0, 2) "spectral norm of matrix Z0 (largest sv)";
  Real normZ2=Matrices.norm(Z2, 2) "spectral norm of matrix Z2 (largest sv)";
algorithm
  ku := (normZ0*normQ + 2*sqrt(normZ0*normZ2)*normA + normZ2*normG)/normX;
end k_care_u;

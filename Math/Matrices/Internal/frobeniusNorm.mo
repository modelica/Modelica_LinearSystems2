within Modelica_LinearSystems2.Math.Matrices.Internal;
function frobeniusNorm "Return the Frobenius norm of a matrix"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Input matrix";
  output Real result=if min(size(A))>0 then sqrt(sum(A.*A)) else -1e100
    "frobenius norm of matrix A";

algorithm
end frobeniusNorm;

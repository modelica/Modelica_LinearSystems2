within Modelica_LinearSystems2.Math.Matrices;
function frobeniusNorm "Return the Frobenius norm of a matrix"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Input matrix";
  output Real result=0.0 "frobenius norm of matrix A";

algorithm
  for i1 in 1:size(A, 1) loop
    for i2 in 1:size(A, 2) loop
      result := result + abs(A[i1, i2]);
    end for;
  end for;
  result := sqrt(result);
end frobeniusNorm;

within Modelica_LinearSystems2.ComplexMathAdds.Internal;
function frobeniusNorm "Return the Frobenius norm of a matrix"
  extends Modelica.Icons.Function;
  import Modelica.ComplexMath;

  input Complex A[:,:] "Input matrix";
  output Real result=0.0 "frobenius norm of matrix A";
algorithm
  for i1 in 1:size(A, 1) loop
    for i2 in 1:size(A, 2) loop
      result := result + ComplexMath.real(A[i1, i2]*ComplexMath.conj(A[i1, i2]));
    end for;
  end for;
  result := sqrt(result);
end frobeniusNorm;

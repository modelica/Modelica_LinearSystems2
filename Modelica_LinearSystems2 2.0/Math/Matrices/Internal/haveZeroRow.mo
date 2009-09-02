within Modelica_LinearSystems2.Math.Matrices.Internal;
function haveZeroRow
  "Boolean output is true if at least one matrix row is zero vector"

  input Real A[:,:]=fill(
        0,
        0,
        0) "input matrix";
  output Boolean result "is true if A has at least one zero row";

protected
  Boolean h;
  Integer i;
algorithm

  i := 1;
  h := false;
  while i < size(A, 1) + 1 and not h loop
    h := Modelica.Math.Vectors.isEqual(
        zeros(size(A, 2)),
        A[i, :],
        100*Modelica.Constants.eps);
    i := i + 1;
  end while;
  result := h;

end haveZeroRow;

within Modelica_LinearSystems2.Math.Matrices;
function det "Determinant of a matrix (computed by LU decomposition)"

  extends Modelica.Icons.Function;
  input Real A[:,size(A, 1)];
  output Real result "Determinant of matrix A";
protected
  Real LU[size(A, 1),size(A, 1)];
  Integer pivots[size(A, 1)];

algorithm
  if size(LU, 1) > 0 then
    (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(A);
    result := product(LU[i, i] for i in 1:size(A, 1))*product(if pivots[i] == i then
            1 else -1 for i in 1:size(pivots, 1));
  else
    result := -1e100;
  end if;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<b>det</b>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the determinant of matrix A
computed by a LU decomposition.
Usually, this function should never be used, because
there are nearly always better numerical algorithms
as by computing the determinant. E.g., use function
<a href=\"modelica://Modelica.Math.Matrices.rank\">Matrices.rank</a>
to compute the rank of a matrix.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.rank\">Matrices.rank</a>,
<a href=\"modelica://Modelica.Math.Matrices.solve\">Matrices.solve</a>
</p>
</html>"));
end det;

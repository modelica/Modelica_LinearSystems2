within Modelica_LinearSystems2.Math.Matrices;
function equalityLeastSquares
  "Solve a linear equality constrained least squares problem"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Minimize |A*x - a|^2";
  input Real a[size(A, 1)];
  input Real B[:,size(A, 2)] "Subject to B*x=b";
  input Real b[size(B, 1)];
  output Real x[size(A, 2)] "Solution vector";

protected
  Integer info;
algorithm
  if Modelica.Math.Matrices.rank(cat(
      1,
      A,
      B)) == size(A, 2) and min(size(A)) > 0 then

    (x,info) := Modelica.Math.Matrices.LAPACK.dgglse_vec(
      A,
      a,
      B,
      b);

    assert(info == 0, "Solving a linear equality-constrained least squares problem
with function \"Matrices.equalityLeastSquares\" failed.");
  else
    if size(A, 2) == 0 then
      x := fill(0, 0);
    else
      x := Modelica_LinearSystems2.Math.Matrices.leastSquares(cat(
        1,
        A,
        B), cat(
        1,
        a,
        b));
    end if;
  end if;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
x = Matrices.<b>equalityLeastSquares</b>(A,a,B,b);
</pre></blockquote>

<h4>Description</h4>
<p>
This function returns the
solution <b>x</b> of the linear equality-constrained least squares problem:
</p>
<blockquote>
<p>
min|<b>A</b>*<b>x</b> - <b>a</b>|^2 over <b>x</b>, subject to <b>B</b>*<b>x</b> = <b>b</b>
</p>
</blockquote>

<p>
It is required that the dimensions of A and B fulfill the following
relationship:
</p>

<blockquote>
size(B,1) &le; size(A,2) &le; size(A,1) + size(B,1)
</blockquote>

<h4>Note</h4>
<p>
The solution is computed with the LAPACK function &quot;dgglse&quot;
using the generalized RQ factorization under the assumptions that
B has full row rank (= size(B,1)) and the matrix [A;B] has
full column rank (= size(A,2)). In this case, the problem
has a unique solution.
</p>
</html>"));
end equalityLeastSquares;

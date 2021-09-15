within Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices;
function norm "Returns the norm of a matrix"
  extends Modelica.Icons.Function;
  import Modelica.ComplexMath;
  import Complex;
  import Modelica_LinearSystems2;

  input Complex A[:, :] "Input matrix";
  input Real p(min=1) = 2
    "Type of p-norm (only allowed: 1, 2 or Modelica.Constants.inf)";
  output Real result=0.0 "p-norm of matrix A";

algorithm
  if p == 1 then
    // column sum norm
    for i in 1:size(A, 2) loop
      result := max(result, sum(ComplexMath.abs(A[:, i])));

    end for;
  elseif p == 2 then
    // largest singular value
    result := max(Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_singularValues(A));
  elseif p==3 then
    // Frobenius norm
    result := Modelica_LinearSystems2.ComplexMathAdds.Internal.frobeniusNorm(A);
  elseif p == Modelica.Constants.inf then
    // row sum norm
    for i in 1:size(A, 1) loop
      result := max(result, sum(ComplexMath.abs(A[i, :])));
    end for;
  else
    assert(false, "Optional argument \"p\" of function \"norm\" must be
1, 2 or Modelica.Constants.inf");
  end if;
  annotation ( Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<strong>norm</strong>(A);
Matrices.<strong>norm</strong>(A, p=2);
</pre></blockquote>

<h4>Description</h4>
<p>
The function call &quot;<code>Matrices.norm(A)</code>&quot; returns the
2-norm of matrix A, i.e., the largest singular value of A.<br>
The function call &quot;<code>Matrices.norm(A, p)</code>&quot; returns the
p-norm of matrix A. The only allowed values for p are</p>
<ul>
<li> &quot;p=1&quot;: the largest column sum of A</li>
<li> &quot;p=2&quot;: the largest singular value of A</li>
<li> &quot;p=Modelica.Constants.inf&quot;: the largest row sum of A</li>
</ul>
<p>
Note, for any matrices A1, A2 the following inequality holds:
</p>
<blockquote><pre>
Matrices.<strong>norm</strong>(A1+A2,p) &le; Matrices.<strong>norm</strong>(A1,p) + Matrices.<strong>norm</strong>(A2,p)
</pre></blockquote>
<p>
Note, for any matrix A and vector v the following inequality holds:
</p>
<blockquote><pre>
Vectors.<strong>norm</strong>(A*v,p) &le; Matrices.<strong>norm</strong>(A,p)*Vectors.<strong>norm</strong>(A,p)
</pre></blockquote>
</html>"));
end norm;

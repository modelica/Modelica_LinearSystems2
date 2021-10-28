within Modelica_LinearSystems2.Math.Matrices;
function conditionNumber "Calculate the condition number norm(A)*norm(inv(A))"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Input matrix";
  input Real p(min=1) = 2
    "Type of p-norm (only allowed: 1, 2 or Modelica.Constants.inf)";
  output Real result=0.0 "p-norm of matrix A";

protected
  Real eps=1e-25;
  Real s[size(A, 1)] "singular values";

algorithm
  if min(size(A)) > 0 then
  if p == 2 then
    s := Modelica.Math.Matrices.singularValues(A);
    if min(s) < eps then
result := Modelica.Constants.inf;
    else
result := max(s)/min(s);
    end if;
  else
    result := Modelica.Math.Matrices.norm(A, p)*Modelica.Math.Matrices.norm(
      Modelica.Math.Matrices.inv(A), p);
    end if;
    end if;


  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.conditionNumber instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
r = Matrices.<strong>conditionNumber</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the the condition number
(norm(<strong>A</strong>) * norm(inv(<strong>A</strong>))) of a general real matrix <strong>A</strong>,
in either the 1-norm, 2-norm or the infinity-norm. In the case of 2-norm
the result is the ratio of the largest to the smallest singular value to <strong>A</strong>.
</p>

<h4>Example</h4>
<blockquote><pre>
  A = [1, 2
       2, 1];
  r = conditionNumber(A);

results in:

  r = 3.0
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.rcond\">Matrices.rcond</a>
</p>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end conditionNumber;

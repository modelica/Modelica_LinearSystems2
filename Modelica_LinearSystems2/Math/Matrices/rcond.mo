within Modelica_LinearSystems2.Math.Matrices;
function rcond "Reciprocal condition number"
  extends Modelica.Icons.Function;
  input Real A[:,size(A,1)] "Square real matrix";
  input Boolean inf = false
    "Is true if infinity norm is used and false for 1-norm";
  output Real rcond "Reciprocal condition number of A";
  output Integer info "Information";
protected
  Real LU[:,:] "LU factorization of matrix A, returned by dgetrf";
  Real anorm "Norm of matrix A";
  String normspec= if inf then "I" else "1" "Specifies the norm 1 or inf";

algorithm
  if min(size(A)) > 0 then
    (LU,,info) := Modelica.Math.Matrices.LAPACK.dgetrf(A);
    anorm := Modelica_LinearSystems2.Math.Matrices.LAPACK.dlange(A,normspec);
    (rcond,info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgecon(LU,inf,anorm);
  else
    rcond := Modelica.Constants.inf;
    info := 0;
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
        r = Matrices.<strong>rcond</strong>(A);
(r, info) = Matrices.<strong>rcond</strong>(A, false);
</pre></blockquote>

<h4>Description</h4>
<p>
This function estimates the reciprocal of the condition number
(norm(<strong>A</strong>) * norm(inv(<strong>A</strong>))) of a general real matrix <strong>A</strong>,
in either the 1-norm or the infinity-norm, using the LAPACK function DGECON.
</p>

<h4>Example</h4>
<blockquote><pre>
  A = [1, 2
       2, 1];
  r = rcond(A);

  results in:

  r = 0.3333
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.conditionNumber\">Matrices.conditionNumber</a>
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
end rcond;

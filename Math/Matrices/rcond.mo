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
  String normspec= if inf then "I" else "1" "specifies the norm 1 or inf";

algorithm
  if min(size(A)) > 0 then
    (LU,,info) := Modelica.Math.Matrices.LAPACK.dgetrf(A);
    anorm := Modelica_LinearSystems2.Math.Matrices.LAPACK.dlange(A,normspec);
    (rcond,info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgecon(LU,inf,anorm);
  else
    rcond := Modelica.Constants.inf;
    info := 0;
  end if;

  annotation (Documentation(info="<HTML>
  <h4>Syntax</h4>
<blockquote><pre>
r = Matrices.<b>rcond</b>(A);
(r, info) = Matrices.<b>rcond</b>(A, false);
</pre></blockquote>
<h4>Description</h4>
<p>
This function estimates the reciprocal of the condition number (norm(A) * norm(inv(A))) of a general real matrix A, in either the 1-norm or the infinity-norm, using the LAPACK function DGECON.   
</p>
<p>
<h4>Example</h4>
<blockquote><pre>
  A = [1, 2
       2, 1];
  r = rcond(A);
  
  results in:
  
  r = 0.3333
</pre></blockquote>
</p>
<h4>See also</h4>
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.conditionNumber\">Matrices.conditionNumber</a>
</HTML>", revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end rcond;

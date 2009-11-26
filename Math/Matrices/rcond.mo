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
(LU, pivots)       = Matrices.<b>LU</b>(A);
(LU, pivots, info) = Matrices.<b>LU</b>(A);
</pre></blockquote>
<h4>Description</h4>
<p>
This function call returns the
LU decomposition of a \"Real[m,n]\" matrix A, i.e.,
</p>
<blockquote>
<p>
<b>P</b>*<b>L</b>*<b>U</b> = <b>A</b>
</p>
</blockquote>
<p>
where <b>P</b> is a permutation matrix (implicitely
defined by vector <code>pivots</code>),
<b>L</b> is a lower triangular matrix with unit
diagonal elements (lower trapezoidal if m &gt; n), and
<b>U</b> is an upper triangular matrix (upper trapezoidal if m &lt; n).
Matrices <b>L</b> and <b>U</b> are stored in the returned
matrix <code>LU</code> (the diagonal of <b>L</b> is not stored).
With the companion function 
<a href=\"Modelica://Modelica.Math.Matrices.LU_solve\">Matrices.LU_solve</a>,
this decomposition can be used to solve
linear systems (<b>P</b>*<b>L</b>*<b>U</b>)*<b>x</b> = <b>b</b> with different right
hand side vectors <b>b</b>. If a linear system of equations with
just one right hand side vector <b>b</b> shall be solved, it is
more convenient to just use the function
<a href=\"Modelica://Modelica.Math.Matrices.solve\">Matrices.solve</a>.
</p>
<p>
The optional third (Integer) output argument has the following meaning:
<table border=0 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\">info = 0:</td
      <td valign=\"top\">successful exit</td></tr>
  <tr><td valign=\"top\">info &gt; 0:</td>
      <td valign=\"top\">if info = i, U[i,i] is exactly zero. The factorization
          has been completed, <br> 
          but the factor U is exactly
          singular, and division by zero will occur<br> if it is used
          to solve a system of equations.</td></tr>
</table>
</p>
<p>
The LU factorization is computed
with the LAPACK function \"dgetrf\",
i.e., by Gaussian elemination using partial pivoting
with row interchanges. Vector \"pivots\" are the
pivot indices, i.e., for 1 &le; i &le; min(m,n), row i of 
matrix A was interchanged with row pivots[i].
</p>
<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3; 
                 3,4,5;
                 2,1,4];
  Real b1[3] = {10,22,12};
  Real b2[3] = { 7,13,10};
  Real    LU[3,3];
  Integer pivots[3];
  Real    x1[3];
  Real    x2[3];
<b>algorithm</b>
  (LU, pivots) := Matrices.LU(A);
  x1 := Matrices.LU_solve(LU, pivots, b1);  // x1 = {3,2,1}
  x2 := Matrices.LU_solve(LU, pivots, b2);  // x2 = {1,0,2}
</pre></blockquote>
<h4>See also</h4>
<a href=\"Modelica://Modelica.Math.Matrices.LU_solve\">Matrices.LU_solve</a>, 
<a href=\"Modelica://Modelica.Math.Matrices.solve\">Matrices.solve</a>,
</HTML>"));
end rcond;

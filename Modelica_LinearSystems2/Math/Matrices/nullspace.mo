within Modelica_LinearSystems2.Math.Matrices;
function nullspace "Orthonormal nullspace of a matrix"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Input matrix";
  output Real Z[size(A, 2),:] "Orthonormal nullspace of matrix A";
  output Integer nullity "Nullity, i.e. the dimension of the nullspace";

protected
  Real V[size(A, 2),size(A, 2)] "Right orthogonal matrix ";
  Real sigma[min(size(A, 1), size(A, 2))] "Singular values";
  Integer rank "Rank of matrix A";
  Real eps "Tolerance for rank determination";
  Integer n=min(size(A, 1), size(A, 2));
  Integer i=n;

algorithm
  (sigma,,V) := Modelica.Math.Matrices.singularValues(A);
  V := transpose(V);
  // rank computation
  eps := max(size(A, 1), size(A, 2))*max(sigma)*Modelica.Constants.eps;
  rank := 0;
  if n > 0 then
    while i > 0 loop
      if sigma[i] > eps then
        rank := i;
        i := 0;
      end if;
      i := i - 1;
    end while;
  end if;
  Z := V[:,rank + 1:size(A,2)];// nullspace computation
  nullity := size(A,2) - rank;// nullity

  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.nullSpace instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
           Z = Matrices.<strong>nullspace</strong>(A);
(Z, nullity) = Matrices.<strong>nullspace</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates an orthonormal basis <strong>Z</strong>=[<strong>z</strong>_1, <strong>z</strong>_2, ...]
of the nullspace of a matrix <strong>A</strong>, i.e. <strong>A</strong>*<strong>z</strong>_i=0.
The nullspace is obtained by svd method. That is, matrix <strong>A</strong> is decomposed
into the matrices <strong>S</strong>, <strong>U</strong>, <strong>V</strong>:
</p>
<blockquote>
  <strong>A</strong> = <strong>U</strong><strong>S</strong><strong>V</strong><sup>T</sup>
</blockquote>
<p>
with the orthonormal matrices <strong>U</strong> and <strong>V</strong> and the matrix <strong>S</strong> with
</p>
<blockquote>
  <strong>S</strong> = [<strong>S</strong>1, <strong>0</strong>]
  <strong>S</strong>1 = [diag(s); <strong>0</strong>]
</blockquote>
<p>
with the singular values <strong>s</strong>={s1, s2, ..., sr} of <strong>A</strong> and r=rank(<strong>A</strong>).
Note, that <strong>S</strong> has the same size as <strong>A</strong>. Since <strong>U</strong> and <strong>V</strong> are
orthonormal, we may write
</p>
<blockquote>
  <strong>U</strong><sup>T</sup>*<strong>A</strong>*<strong>V</strong> = [<strong>S</strong>1, <strong>0</strong>].
</blockquote>
<p>
Matrix <strong>S</strong>1 obviously has full column rank and therefore, the left n-r rows
(n is the number of columns of <strong>A</strong> or <strong>S</strong>) of matrix <strong>V</strong> span
a nullspace of <strong>A</strong>.
</p>
<p>
The nullity of matrix <strong>A</strong> is the dimension of the nullspace of <strong>A</strong>.
In view of the above, it becomes clear that nullity holds
</p>
<blockquote><pre>
nullity = n - r
</pre></blockquote>
<p>
with
</p>
<blockquote>
n = number of columns of matrix <strong>A</strong> and <br>
r = rank(<strong>A</strong>).
</blockquote>

<h4>Example</h4>
<blockquote><pre>
  A = [1, 2,  3, 1;
       3, 4,  5, 2;
      -1, 2, -3, 3];
  (Z, nullity) = nullspace(A);

  results in:

  Z=[0.1715;
    -0.686;
     0.1715;
     0.686]

  nullity = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.singularValues\">Modelica.Math.Matrices.singularValues</a>
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
end nullspace;

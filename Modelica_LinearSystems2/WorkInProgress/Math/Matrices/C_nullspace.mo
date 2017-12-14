within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_nullspace "Orthonormal nullspace of a complex matrix"
  import Modelica_LinearSystems2;
  import Complex;
  extends Modelica.Icons.Function;

  input Complex A[:,:] "Input matrix";
  output Complex Z[size(A, 2),:]
    "Orthonormal nullspace of matrix A";
  output Integer nullity "Nullity, i.e. the dimension of the nullspace";

protected
  Complex V[size(A, 2),size(A, 2)]
    "Right orthogonal matrix ";
  Real sigma[min(size(A, 1), size(A, 2))] "Singular values";
  Integer rank "Rank of matrix A";
  Real eps "Tolerance for rank determination";
  Integer n=min(size(A, 1), size(A, 2));
  Integer i=n;

algorithm
  (sigma,,V) := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_singularValues(A);
  V :=Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.C_transpose(V);
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

annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
           Z = Matrices.<b>nullspace</b>(A);
                 or
(Z, nullity) = Matrices.<b>nullspace</b>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates a orthonormal basis Z=[z_1, z_2, ...] of the nullspace of a matrix A, i.e. A*z_i=0.
</p>
<p>
The nullspace is obtained by svd method. That is, matrix A is dcomposed into the matrices S, U, V:
</p>
<blockquote><pre>
<b>A</b> = <b>U</b><b>S</b><b>V</b><sup><big>H</big></sup>
</pre></blockquote>
<p>
with the orthonormal matrices U and V and the matrix S with
</p>
<blockquote><pre>
<b>S</b> = [<b>S</b>1, <b>0</b>]
<b>S</b>1 = [diag(s); <b>0</b>]
</pre></blockquote>
<p>
with the singular values s={s1, s2, ..., sr} of A and r=rank(A). Note, that S has the same size as A. Since, U and V are orthonormal, we may write
</p>
<blockquote><pre>
 <b>U</b><sup><big>H</big></sup>*A*V = [<b>S</b>1, <b>0</b>].
</pre></blockquote>
<p>
Matrix S1 obviously has full column rank and therefore, the left n-r rows (n is the number of columns of A or S) of matrix V span a nullspace of A.
</p>
<p>
The nullity of matrix A is the dimension of the nullspace of A. In view of the above, it becomes clear that nullity holds
</p>
<blockquote><pre>
nullity = n - r
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
n = number of columns of matrix A
r = rank(A)
</pre></blockquote>

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
<a href=\"modelica://Modelica.Math.Matrices.singularValues\">Matrices.singularValues</a>
</p>
</html>"));
end C_nullspace;

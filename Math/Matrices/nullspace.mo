within Modelica_LinearSystems2.Math.Matrices;
function nullspace "Orthonormal nullspace of a matrix"
  extends Modelica.Icons.Function;

  input Real A[:,:] "input matrix";
  output Real Z[size(A, 2),:] "orthonormal nullspace of matrix A";
  output Integer nullity "nullity, i.e. the dimension of the nullspace";

protected
  Real V[size(A, 2),size(A, 2)] "Right orthogonal matrix ";
  Real sigma[min(size(A, 1), size(A, 2))] "singular values";
  Integer rank "rank of matrix A";
  Real eps "tolerance for rank determination";
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

annotation (Documentation(info="<html>
  <h4>Syntax</h4>
<blockquote><pre>
Z = Matrices.<b>nullspace</b>(A);<br>
(Z, nullity) = Matrices.<b>nullspace</b>(A);
</pre></blockquote>
<h4>Description</h4>
<p>
This function calculates a orthonormal basis <b>Z</b>=[<b>z</b>_1, <b>z</b>_2, ...] of the nullspace of a matrix <b>A</b>, i.e. <b>A</b>*<b>z</b>_i=0.
</p>
The nullspace is obtained by svd method. That is, matrix <b>A</b> is decomposed into the matrices <b>S</b>, <b>U</b>, <b>V</b>: 
<blockquote><pre>
        T
 <b>A</b> = <b>U</b><b>S</b><b>V</b>
</pre></blockquote>
with the orthonormal matrices <b>U</b> and <b>V</b> and the matrix <b>S</b> with
<blockquote><pre>
 <b>S</b> = [<b>S</b>1, <b>0</b>]
 <b>S</b>1 = [diag(s); <b>0</b>]
</pre></blockquote>
with the singular values <b>s</b>={s1, s2, ..., sr} of <b>A</b> and r=rank(<b>A</b>). Note, that <b>S</b> has the same size as <b>A</b>. Since, <b>U</b> and <b>V</b> are orthonormal, we may write
<blockquote><pre>
  T
 <b>U</b>*<b>A</b>*<b>V</b> = [<b>S</b>1, <b>0</b>].
</pre></blockquote>
Matrix <b>S</b>1 obviously has full column rank and therefore, the left n-r rows (n is the number of columns of <b>A</b> or <b>S</b>) of matrix <b>V</b> span a nullspace of <b>A</b>.
</p>
<p>
The nullity of matrix <b>A</b> is the dimension of the nullspace of <b>A</b>. In view of the above, it becomes clear that nullity holds
<blockquote><pre>
 nullity = n - r
</pre></blockquote>
with
<blockquote><pre>
 n = number of columns of matrix <b>A</b>
 r = rank(<b>A</b>)
</pre></blockquote>

</p>
<p>
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
</p>
<h4>See also</h4>
<a href=\"modelica://Modelica.Math.Matrices.singularValues\">Matrices.singularValues</a>
</html>", revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end nullspace;

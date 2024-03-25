within Modelica_LinearSystems2.Math.Matrices.Internal;
function solveSymRight
  "Solve real system of linear equations X*A=B in X where A is symmetrix positive definite"

  extends Modelica.Icons.Function;
  import MatricesMSL = Modelica.Math.Matrices;

  input Real A[:,size(A, 1)] "Matrix A of X*A = B";
  input Real B[:,:] "Matrix B of X*op(A) = B";
  input Boolean isCholesky=false
    "True if the A is already lower Cholesky factor";
  input Boolean upper=false "True if A is upper triAngular";
  output Real X[size(B, 1),size(B, 2)]=B "Matrix X such that X*A = B";

protected
  Real H[size(A, 1),size(A, 2)];
//  Real AA[size(A,1),size(A,1)]=if upper then A else transpose(A);
  Real AA[size(A,1),size(A,1)]=A;

algorithm
  if not isCholesky then
    H := MatricesMSL.cholesky(AA, upper);
  else
    H := AA;
  end if;
  H := symmetric(H);
  X := MatricesMSL.LAPACK.dtrsm(H, X, 1, true, true, false, false);
  X := MatricesMSL.LAPACK.dtrsm(H, X, 1, true, false, false, false);

  annotation (Documentation(info="<html>
<p>
This function solves the equation
</p>

<blockquote><pre>
<strong>X</strong>*<strong>A</strong> = <strong>B</strong>
</pre></blockquote>

<p>
where matrix&nbsp;<strong>A</strong> with symmetric positive definite matrix.
The calculation is rather efficient since symmetrie and decomposition
of positive definite matrices is exploited.
</p>
<p>
Due to symmetry, the matrix&nbsp;<strong>A</strong> is uniquely defined by
a&nbsp;triangle, i.e. the upper or the lower triangular matrix.
It is assumed, that the input to describe&nbsp;<strong>A</strong> is either
a&nbsp;Cholesky factor or part of matrix&nbsp;<strong>A</strong> itself.
This is defined by the user with the boolean inputs <code>isTriangular</code>
and <code>upper</code> which is true when&nbsp;<strong>A</strong> is already
Cholesky factor and when&nbsp;<strong>A</strong> is upper triangular,
respectively.
</p>
<p>
Considering the Cholesky decomposition
</p>

<blockquote><pre>
       T
<strong>A</strong> = <strong>L</strong>*<strong>L</strong>
</pre></blockquote>

<p>
with lower triangular matrix&nbsp;<strong>L</strong> the equation above
could be rewritten as
</p>

<blockquote><pre>
     T
<strong>X</strong>*<strong>L</strong>*<strong>L</strong> = <strong>B</strong>
</pre></blockquote>

<p>
which is solved with BLAS function <em>dtrmm</em> applied to a&nbsp;upper
triangular matrix and subsequently to a&nbsp;lower triangular matrix.
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
end solveSymRight;

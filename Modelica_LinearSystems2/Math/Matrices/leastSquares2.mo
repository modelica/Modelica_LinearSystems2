within Modelica_LinearSystems2.Math.Matrices;
function leastSquares2
  "Solve overdetermined or underdetermined real system of linear equations A*X=B in a least squares sense (A may be rank deficient)"
  extends Modelica.Icons.Function;
  input Real A[:, :] "Matrix A";
  input Real B[size(A, 1),:] "Matrix B";
  input Real rcond=100*Modelica.Constants.eps
    "Reciprocal condition number to estimate rank of A";
  output Real X[size(A, 2), size(B,2)]
    "Matrix X such that min|A*X-B|^2 if size(A,1) >= size(A,2) or min|X|^2 and A*X=B, if size(A,1) < size(A,2)";
  output Integer rank "Rank of A";
protected
  Integer info;
  Real XX[max(size(A,1),size(A,2)), size(B,2)];
algorithm
  (XX,info,rank) := LAPACK.dgelsx(A, B, rcond);
  X := XX[1:size(A,2), :];
  assert(info == 0, "Solving an overdetermined or underdetermined linear system of
equations with function \"Matrices.leastSquares2\" failed.");
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.leastSquares2 instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
X = Matrices.<b>leastSquares2</b>(A,B);
</pre></blockquote>


<h4>Description</h4>
<p>
Returns a solution of equation A*X = B in a least
square sense (A may be rank deficient):
</p>
<pre>
  minimize | A*X - B |
</pre>

<p>
Several different cases can be distinguished (note, <b>rank</b> is an
output argument of this function):
</p>

<p>
<b>size(A,1) = size(A,2)</b>
</p>

<p> A solution is returned for a regular, as well as a singular matrix A:
</p>

<ul>
<li> <b>rank</b> = size(A,1):<br>
     A is <b>regular</b> and the returned solution X fulfills the equation
     A*X = B uniquely.</li>

<li> <b>rank</b> &lt; size(A,1):<br>
     A is <b>singular</b> and no unique solution for equation A*X = B exists.
     <ul>
     <li>  If an infinite number of solutions exists, the one is selected that fulfills
           the equation and at the same time has the minimum norm |x| for all solution
           vectors that fulfill the equation.</li>
     <li>  If no solution exists, X is selected such that |A*X - B| is as small as
           possible (but A*X - B is not zero).</li>
     </ul>
</ul>

<p>
<b>size(A,1) &gt; size(A,2):</b>
</p>

<p>
The equation A*X = B has no unique solution. The solution X is selected such that
|A*X - B| is as small as possible. If rank = size(A,2), this minimum norm solution is
unique. If rank &lt; size(A,2), there are an infinite number of solutions leading to the
same minimum value of |A*X - B|. From these infinite number of solutions, the one with the
minimum norm |X| is selected. This gives a unique solution that minimizes both
|A*X - B| and |X|.
</p>

<p>
<b>size(A,1) &lt; size(A,2):</b>
</p>

<ul>
<li> <b>rank</b> = size(A,1):<br>
     There are an infinite number of solutions that fulfill the equation A*X = B.
     From this infinite number, the unique solution is selected that minimizes |X|.
     </li>

<li> <b>rank</b> &lt; size(A,1):<br>
     There is either no solution of equation A*X = B, or there are again an infinite
     number of solutions. The unique solution X is returned that minimizes
      both |A*X - B| and |X|.</li>
</ul>


<h4>Note</h4>
<p>
The solution is computed with the LAPACK function &quot;dgelsx&quot;,
i.e., QR or LQ factorization of A with column pivoting.
</p>


<h4>Algorithmic details</h4>
<p>
The function first computes a QR factorization with column pivoting:
</p>

<pre>
      A * P = Q * [ R11 R12 ]
                  [  0  R22 ]
</pre>

<p>
with R11 defined as the largest leading submatrix whose estimated
condition number is less than 1/rcond.  The order of R11, <b>rank</b>,
is the effective rank of A.
</p>

<p>
Then, R22 is considered to be negligible, and R12 is annihilated
by orthogonal transformations from the right, arriving at the
complete orthogonal factorization:
</p>

<pre>
     A * P = Q * [ T11 0 ] * Z
                 [  0  0 ]
</pre>

<p>
The minimum-norm solution is then
</p>

<pre>
     X = P * Z' [ inv(T11)*Q1'*B ]
                [        0       ]
</pre>

<p>
where Q1 consists of the first &quot;rank&quot; columns of Q.
</p>


<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.leastSquares\">Matrices.leastSquares</a>
(same as leastSquares2, but with a right hand side vector),
<a href=\"modelica://Modelica.Math.Matrices.solve2\">Matrices.solve2</a>
(for square, regular matrices A)
</p>
</html>"));
end leastSquares2;

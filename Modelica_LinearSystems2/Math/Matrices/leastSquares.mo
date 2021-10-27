within Modelica_LinearSystems2.Math.Matrices;
function leastSquares
  "Solve overdetermined or underdetermined real system of linear equations A*x=b in a least squares sense (A may be rank deficient)"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Matrix A";
  input Real b[size(A, 1)] "Vector b";
  output Real x[size(A, 2)]
    "Vector x such that min|A*x-b|^2 if size(A,1) >= size(A,2) or min|x|^2 and A*x=b, if size(A,1) < size(A,2)";

protected
  Integer rank;
algorithm
  (x,rank) :=Modelica.Math.Matrices.leastSquares(A,b);
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
x = Matrices.<strong>leastSquares</strong>(A,b);
</pre></blockquote>
<h4>Description</h4>
<p>
A linear system of equations A*x = b has no solutions or infinitely
many solutions if A is not square. Function &quot;leastSquares&quot; returns
a solution in a least squarse sense:
</p>
<pre>
  size(A,1) &gt; size(A,2):  returns x such that |A*x - b|^2 is a minimum
  size(A,1) = size(A,2):  returns x such that A*x = b
  size(A,1) &lt; size(A,2):  returns x such that |x|^2 is a minimum for all
                          vectors x that fulfill A*x = b
</pre>

<h4>Note</h4>
<p>
The solution is computed with the LAPACK function &quot;dgelsx&quot;,
i.e., QR or LQ factorization of A with column pivoting.
If A does not have full rank,
the solution is not unique and from the infinitely many solutions
the one is selected that minimizes both |x|^2 and |A*x - b|^2.
</p>
</html>"));
end leastSquares;

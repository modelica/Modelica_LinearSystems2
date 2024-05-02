within Modelica_LinearSystems2.Math.Matrices;
function QR
  "QR decomposition of a rectangular matrix without column pivoting (A = Q*R). Return the full square Q-matrix"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Rectangular matrix with size(A,1) >= size(A,2)";
  output Real Q[size(A, 1),size(A, 2)]
    "Rectangular matrix with orthonormal columns such that Q*R=A[:,p]";
  output Real R[min(size(A, 1), size(A, 2)),size(A, 2)]
    "Square upper triangular matrix";

  output Real tau[min(size(A, 1), size(A, 2))] "Scalar factors of the elementary reflectors";
  output Real Q2[size(A,1),size(A, 1)] "Orthogonal matrix defined as the product of elementary reflectors";
protected
  Integer nrow=size(A, 1);
  Integer ncol=size(A, 2);
  Integer minrowcol=min(nrow, ncol);

algorithm
  if minrowcol > 0 then

    (Q,tau) := Modelica.Math.Matrices.LAPACK.dgeqrf(A);

  // determine R
    R := zeros(minrowcol, ncol);
    for i in 1:minrowcol loop
      for j in i:ncol loop
        R[i, j] := Q[i, j];
      end for;
    end for;

    //Q2 := Modelica_LinearSystems2.Math.Matrices.LAPACK.dorgqr(Q, tau);// would return the economic size Q matrix

    Q2 := Modelica_LinearSystems2.Math.Matrices.Internal.multiplyWithOrthogonalQ_qr(identity(size(A, 1)), Q, tau, "L", "N");

  else
    Q := fill(1, size(A, 1), size(A, 2));
    R := fill(0, min(size(A, 1), size(A, 2)), size(A, 2));
  end if;

  annotation ( Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(Q, R, tau, Q2) = Matrices.<strong>QR</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function returns the QR decomposition of
a rectangular matrix <strong>A</strong> (the number of columns of <strong>A</strong>
must be less than or equal to the number of rows):
</p>
<blockquote>
<strong>Q</strong> * <strong>R</strong> = <strong>A</strong>
</blockquote>
<p>
where <strong>Q</strong> is a rectangular matrix that has orthonormal columns and
has the same size as A (<strong>Q</strong><sup>T</sup><strong>Q</strong>=<strong>I</strong>),
<strong>R</strong> is a square, upper triangular matrix.
Matrix <strong>R</strong> has the following important properties:
</p>
<ul>
<li> The absolute value of a diagonal element of <strong>R</strong> is the largest
     value in this row, i.e.,
     abs(R[i,i]) &ge; abs(R[i,j]).</li>
<li> The diagonal elements of <strong>R</strong> are sorted according to size, such that
     the largest absolute value is abs(R[1,1]) and
     abs(R[i,i]) &ge; abs(R[j,j]) with i &lt; j. </li>
</ul>
<p>
This means that if abs(R[i,i]) &le; &epsilon; then abs(R[j,k]) &le; &epsilon;
for j &ge; i, i.e., the i-th row up to the last row of <strong>R</strong> have
small elements and can be treated as being zero.
This allows to, e.g., estimate the row-rank
of <strong>R</strong> (which is the same row-rank as <strong>A</strong>). Furthermore,
<strong>R</strong> can be partitioned in two parts
</p>
<blockquote><pre>
<strong>A</strong> = <strong>Q</strong> * [<strong>R</strong><sub>1</sub>, <strong>R</strong><sub>2</sub>;
          <strong>0</strong>,  <strong>0</strong>]
</pre></blockquote>
<p>
where <strong>R</strong><sub>1</sub> is a regular, upper triangular matrix.
</p>
<p>
Note, the solution is computed with the LAPACK function \"dgeqrf\",
i.e., by a&nbsp;QR factorization without column pivoting.
If <strong>Q</strong> is not needed, the function may be
called as: <code>(,R) = QR(A)</code>.
</p>
<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1, 2, 3;
                 3, 4, 5;
                 2, 1, 4];
  Real R[3,3];
<strong>algorithm</strong>
  (,R) := Matrices.QR(A);
  // R = [-3.74.., -4.27.., -6.94..;
           0.0   , -1.64.., -0.17..;
           0.0   ,  0.0   , -1.29..]
</pre></blockquote>
</html>"));
end QR;

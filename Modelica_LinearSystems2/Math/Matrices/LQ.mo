within Modelica_LinearSystems2.Math.Matrices;
function LQ
  "LQ decomposition of a rectangular matrix without column pivoting (A = L*Q)"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:] "Rectangular matrix";
  output Real L[size(A, 1),size(A,1)]
    "Rectangular matrix containing the lower triangular matrix";
  output Real Q[size(A,1),size(A, 2)]
    "Rectangular matrix with orthonormal columns such that L*Q=A";

protected
  Integer nrow=size(A, 1);
  Integer ncol=size(A, 2);
  Real tau[min(nrow, ncol)];
  Real LQ[nrow,ncol];
  Integer info;

algorithm
  assert(nrow <= ncol, "\nInput matrix A[" + String(nrow) + "," + String(ncol) + "] has more rows as columns.
  This is not allowed when calling Matrices.LQ(A).");
  if nrow > 0 then
    (LQ,tau,info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgelqf(A);
    assert(info == 0, "LQ-computation failed in function \"LQ\"");
  // determine R
    L :=  Matrices.triangle(LQ[:, 1:nrow], false);

    (Q, info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dorglq(LQ, tau);
    assert(info==0,"Computation of matrix Q failed ind function \"LQ\"");

  else
    L := fill(0, nrow, nrow);
    Q := fill(1, nrow, ncol);

  end if;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(Q,R,p) = Matrices.<strong>QR</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function returns the QR decomposition of
a rectangular matrix <strong>A</strong> (the number of columns of <strong>A</strong>
must be less than or equal to the number of rows):
</p>
<blockquote>
  <strong>Q</strong>*<strong>R</strong> = <strong>A</strong>[:,<strong>p</strong>]
</blockquote>
<p>
where <strong>Q</strong> is a rectangular matrix that has orthonormal columns and
has the same size as A (<strong>Q</strong><sup>T</sup><strong>Q</strong>=<strong>I</strong>),
<strong>R</strong> is a square, upper triangular matrix and <strong>p</strong> is a permutation
vector. Matrix <strong>R</strong> has the following important properties:
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
<strong>A</strong>[:,<strong>p</strong>] = <strong>Q</strong> * [<strong>R</strong><sub>1</sub>, <strong>R</strong><sub>2</sub>;
              <strong>0</strong>,  <strong>0</strong>]
</pre></blockquote>
<p>
where <strong>R</strong><sub>1</sub> is a regular, upper triangular matrix.
</p>

<h4>Note</h4>
<p>
The solution is computed with the LAPACK functions \"dgeqp3\"
and \"dorgqr\", i.e., by Housholder transformations with
column pivoting. If <strong>Q</strong> is not needed, the function may be
called as: <code>(,R,p) = QR(A)</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real R[3,3];
<strong>algorithm</strong>
  (,R) := Matrices.QR(A);  // R = [-7.07.., -4.24.., -3.67..;
                                    0     , -1.73.., -0.23..;
                                    0     ,  0     ,  0.65..];
</pre></blockquote>
</html>"));
end LQ;

within Modelica_LinearSystems2.Math.Matrices;
function eigenValuesAsRealMatrix
  "Return eigenvalues for a real, nonsymmetric matrix in a Real representation (computation with optional balancing)"

  extends Modelica.Icons.Function;
  input Real A[:, size(A, 1)] "Matrix";
  input Boolean balance=true
    "=true, if A is balanced (pre-scaled) before computation of the eigen values";
  output Real eigenvalues[size(A, 1), 2]
    "Eigenvalues of matrix A (Re: first column, Im: second column)";
protected
  Integer info;
algorithm
  if size(A, 1) > 0 then
    if balance then
      (eigenvalues[:, 1],eigenvalues[:, 2],info) :=
         Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx_eigenValues(A);
    else
      (eigenvalues[:, 1],eigenvalues[:, 2],info) :=
         Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev_eigenValues(A);
    end if;
    assert(info == 0, "Calculating the eigen values with function
\"eigenvaluesAsRealMatrix\" is not possible, since the
numerical algorithm does not converge.");
  end if;
  annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
eigenvalues = Matrices.<strong>eigenValuesAsRealMatrix</strong>(A);
</pre></blockquote>
<h4>Description</h4>
<p>
This function call returns the eigenvalues of a square matrix
<strong>A</strong>. The first column of &quot;eigenvalues&quot; contains the real and the
second column contains the imaginary part of the eigenvalues.
Before calculating the eigenvalues, matrix A is permuted and scaled (balanced)
to improve the computation. For details see the
<a href=\"http://www.netlib.org/lapack/lug/node94.html\">lapack documentation</a>.
</p>
<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real eval[3,2];
<strong>algorithm</strong>
  eval := Matrices.eigenValuesAsRealMatrix(A);  // eval = [ 8.0  , 0;
                                                //         -0.618, 0;
                                                //          1.618, 0];
</pre></blockquote>
<p>
i.e., matrix A has the 3 real eigenvalues 8.0, -0.618, 1.618.
</p>
</html>"));
end eigenValuesAsRealMatrix;

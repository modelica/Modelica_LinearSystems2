within Modelica_LinearSystems2.Math.Matrices;
function eigenValues
  "Compute eigenvalues and eigenvectors for a real, nonsymmetric matrix"

  extends Modelica.Icons.Function;

  input Real A[:, size(A, 1)] "Matrix";
  output Real eigenvalues[size(A, 1), 2]
    "Eigenvalues of matrix A (Re: first column, Im: second column)";

  output Real leftEigenvectors[size(A,1), size(A,2)]
    "Real-valued eigenvector matrix";
 output Real rightEigenvectors[size(A,1), size(A,2)]
    "Real-valued eigenvector matrix";

protected
  Integer info;
  Boolean onlyEigenvalues = false;
algorithm
if size(A,1) > 0 then
  if onlyEigenvalues then
      (eigenvalues[:, 1],eigenvalues[:, 2],info) :=
        Modelica.Math.Matrices.LAPACK.dgeev_eigenValues(A);
     rightEigenvectors :=zeros(size(A, 1), size(A, 1));
     leftEigenvectors :=zeros(size(A, 1), size(A, 1));
  else
      (eigenvalues[:, 1],eigenvalues[:, 2],leftEigenvectors,rightEigenvectors,info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeev(A);
  end if;
  assert(info == 0, "Calculating the eigenvalues with function
\"Matrices.eigenvalues\" is not possible, since the
numerical algorithm does not converge.");
end if;
  annotation (
    Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
                eigenvalues = Matrices.<b>eigenValues</b>(A);
(eigenvalues, eigenvectors) = Matrices.<b>eigenValues</b>(A);
</pre></blockquote>
<h4>Description</h4>
<p>
This function call returns the eigenvalues and 
optionally the (right) eigenvectors of a square matrix 
<b>A</b>. The first column of \"eigenvalues\" contains the real and the
second column contains the imaginary part of the eigenvalues.
If the i-th eigenvalue has no imaginary part, then eigenvectors[:,i] is
the corresponding real eigenvector. If the i-th eigenvalue
has an imaginary part, then eigenvalues[i+1,:] is the conjugate complex
eigenvalue and eigenvectors[:,i] is the real and eigenvectors[:,i+1] is the
imaginary part of the eigenvector of the i-th eigenvalue.
With function 
<a href=\"modelica://Modelica.Math.Matrices.eigenValueMatrix\">Matrices.eigenValueMatrix</a>,
a real block diagonal matrix is constructed from the eigenvalues 
such that 
</p>
<blockquote>
<pre>
A = eigenvectors * eigenValueMatrix(eigenvalues) * inv(eigenvectors)
</pre>
</blockquote>
<p>
provided the eigenvector matrix \"eigenvectors\" can be inverted
(an inversion is possible, if all eigenvalues are different
and no eigenvalue is zero).
</p>
<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3; 
                 3,4,5;
                 2,1,4];
  Real eval;
<b>algorithm</b>
  eval := Matrices.eigenValues(A);  // eval = [-0.618, 0; 
                                    //          8.0  , 0;
                                    //          1.618, 0];
</pre>
</blockquote>
<p>
i.e., matrix A has the 3 real eigenvalues -0.618, 8, 1.618.
</p>
</pre></blockquote>
<h4>See also</h4>
<a href=\"modelica://Modelica.Math.Matrices.eigenValueMatrix\">Matrices.eigenValueMatrix</a>,
<a href=\"modelica://Modelica.Math.Matrices.singularValues\">Matrices.singularValues</a>
</HTML>
"));
end eigenValues;

within Modelica_LinearSystems2.ComplexMathAdds;
function eigenVectors
  "Calculate the right eigenvectors of a linear state space system and write them columnwise in a matrix."
  extends Modelica.Icons.Function;

  import Modelica.ComplexMath.j;

  input Real A[:,size(A, 1)] "Real square matrix";
  output Complex eigvec[size(A, 1),size(A, 2)] "Eigenvectors of the system";
  output Complex eigval[size(A, 1)]=fill(Complex(0), size(A, 1))
    "Eigenvalues of the system";
protected
  Integer info;
  Real eigvecRe[size(A, 1),size(A, 2)];
  Real eigvalRe[size(A, 1)]=fill(0, size(A, 1));
  Real eigvalIm[size(A, 1)]=fill(0, size(A, 1));
  Integer n=size(A, 1);
  Integer i;
algorithm
  if size(A, 1) > 0 then

    (eigvalRe,eigvalIm,eigvecRe,info) := Modelica.Math.Matrices.LAPACK.dgeev(A);
    for i in 1:size(A, 1) loop
      eigval[i].re := eigvalRe[i];
      eigval[i].im := eigvalIm[i];
    end for;

    assert(info == 0, "Calculating the eigen values is not possible"
      + " since the numerical algorithm does not converge.");

    i := 1;
    while i <= n loop
      if abs(eigvalIm[i]) > 0 then
        for ii in 1:n loop
          eigvec[ii, i] := eigvecRe[ii, i] + j*eigvecRe[ii, i + 1];
          eigvec[ii, i + 1] := eigvecRe[ii, i] - j*eigvecRe[ii, i + 1];
        end for;
        i := i + 2;
      else
        for ii in 1:n loop
          eigvec[ii, i] := Complex(1)*eigvecRe[ii, i];
        end for;
        i := i + 1;
      end if;
    end while;

  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(eigvec, eigval) = <strong>eigenVectors</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvectors and optionally (onlyEigenvectors=false) the
eigenvalues of a state space system. The output <code>eigenvectors</code>
is a matrix with the same dimension as matrix <strong>ss.A</strong>.
Just like in <a href=\"modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>,
if the i-th eigenvalue has an imaginary part, then <code>eigenvectors</code>[:,i]
is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the
eigenvector of the i-th eigenvalue.
The eigenvalues are returned as a complex vector <code>eigenvalues</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1,1;-1,-1],
    B=[1;1],
    C=[1,1],
    D=[0]);

  Real eigenvectors[2,2];
  Complex eigenvalues[2];

<strong>algorithm</strong>
  (eigenvectors, eigenvalues) = eigenVectors(ss, true);
// eigenvectors = [0.707, 0; 0, 0.707]
// eigenvalues = {-1 + 1j, -1 - 1j}

          |0.707 |         | 0.707 |
i.e. v1 = |      |,   v2 = |       |
          |0.707i|         |-0.707i|
</pre></blockquote>
</html>"));
end eigenVectors;

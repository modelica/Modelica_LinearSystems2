within Modelica_LinearSystems2.Math.Matrices.Internal;
function reorderRSFc
  "Reorder a real Schur factorization for pole assignment design for continuous systems"
  extends Modelica.Icons.Function;

  input Real T[:,:] "Upper quasi-triangular matrix in Schur canonical form";
  input Real Q[:,size(T, 2)] "Matrix of Schur vectors";
  input Real alphaReal[size(T, 1)]
    "Real part of eigenvalue = alphaReal + i*alphaImag";
  input Real alphaImag[size(T, 1)]
    "Imaginary part of eigenvalue = alphaReal + i*alphaImag";
  input Real alpha
    "Maximum admissible value for real parts of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";

  output Real To[size(T, 1),size(T, 2)] "Reordered Schur form";
  output Real Qo[size(T, 1),size(T, 2)] "Reordered Schur vectors";
  output Real wr[size(T, 2)] "Reordered eigenvalues, real part";
  output Real wi[size(T, 2)] "Reordered eigenvalues, imaginary part";

protected
  Integer n=size(T, 2);
  Boolean select[n]=fill(false, n);
  Integer i;
algorithm
  for i in 1:n loop
    if alphaReal[i] < alpha then
      select[i] := true;
    end if;
  end for;

  (To, Qo, wr, wi) := Modelica.Math.Matrices.LAPACK.dtrsen("E", "V", select, T, Q);

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(To, Qo, wr, wi) = Matrices.Internal.<strong>reorderRSFc</strong>(T, Q, alphaReal, alphaImag, alpha)
</pre></blockquote>

<h4>Description</h4>
<p>
Reorder real Schur form according to <code>alpha</code>, as used e.g. in
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Design.assignPolesMI\">assignPolesMI</a>
for <em>continuous</em> systems.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.Internal.reorderRSF\">reorderRSF</a>
or
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.Internal.reorderRSFd\">reorderRSFd</a>
(for discrete systems).
</p>
</html>"));
end reorderRSFc;

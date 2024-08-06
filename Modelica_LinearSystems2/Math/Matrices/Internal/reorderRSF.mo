within Modelica_LinearSystems2.Math.Matrices.Internal;
function reorderRSF
  "Reorder a real Schur factorization according to a given pattern of the eigenvalues"
  extends Modelica.Icons.Function;

  input Boolean iscontinuous;
  input Real T[:,:] "Upper quasi-triangular matrix in Schur canonical form";
  input Real Q[:,size(T, 2)] "Matrix of Schur vectors";
  input Real alphaReal[size(T, 1)]
    "Real part of eigenvalue = alphaReal + i*alphaImag";
  input Real alphaImag[size(T, 1)]
    "Imaginary part of eigenvalue = alphaReal + i*alphaImag";

  output Real To[size(T, 1),size(T, 2)] "Reordered Schur form";
  output Real Qo[size(T, 1),size(T, 2)] "Reordered Schur vectors";
  output Real wr[size(T, 2)] "Reordered eigenvalues, real part";
  output Real wi[size(T, 2)] "Reordered eigenvalues, imaginary part";

protected
  Integer n=size(T, 2);
  Boolean select[size(T, 2)]=fill(false, size(T, 2));
  Integer i;
algorithm
  if iscontinuous then
    for i in 1:n loop
      if alphaReal[i] < 0 then
        select[i] := true;
      end if;
    end for;
  else
    for i in 1:n loop
      if alphaReal[i]^2 + alphaImag[i]^2 < 1 then
        select[i] := true;
      end if;
    end for;
  end if;

  (To, Qo, wr, wi) := Modelica.Math.Matrices.LAPACK.dtrsen("E", "V", select, T, Q);

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(To, Qo, wr, wi) = Matrices.Internal.<strong>reorderRSF</strong>(iscontinuous, T, Q, alphaReal, alphaImag)
</pre></blockquote>

<h4>Description</h4>
<p>
Reorder real Schur form used for pole assignment design for continuous or discrete systems.
For a&nbsp;continuous system (<code>iscontinuous&nbsp;= true</code>), those eigenvalues
of&nbsp;A will not be modified by the eigenvalue assignment algorithm which real part is negative.
For a&nbsp;discrete system (<code>iscontinuous&nbsp;= false</code>), eigenvalues of moduli less
then one will not be modified.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.Internal.reorderRSFc\">reorderRSFc</a>
for continuous systems or
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.Internal.reorderRSFd\">reorderRSFd</a>
for discrete systems.
</p>
</html>"));
end reorderRSF;

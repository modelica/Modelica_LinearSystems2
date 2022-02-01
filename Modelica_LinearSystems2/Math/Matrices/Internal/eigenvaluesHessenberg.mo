within Modelica_LinearSystems2.Math.Matrices.Internal;
function eigenvaluesHessenberg
  "Compute eigenvalues of an upper Hessenberg form matrix"

  input Real H[:,size(H, 1)];
  output Real alphaReal[size(H, 1)]
    "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
  output Real alphaImag[size(H, 1)]
    "Imaginary part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
  output Integer info=0;
protected
  Integer n=size(H, 1);
  Integer ilo=1;
  Integer ihi=n;
  Real Z[n,n]=fill(
        0,
        n,
        n);

algorithm
  if size(H, 1) > 0 then
    (alphaReal,alphaImag,info) := Modelica.Math.Matrices.LAPACK.dhseqr(H);
  else
    alphaReal := fill(0, size(H, 1));
    alphaImag := fill(0, size(H, 1));
  end if;

  annotation (Documentation(info="<html>
This function uses DHSEQR Lapack-routine to calculate the eigenvalues of an upper Hessenberg form <strong>H</strong>.
Therefore, <strong>H</strong> is reduced to Schur form <strong>T</strong>. The eigenvalues are obtained from the diagonal of <strong>T</strong>.

<p>
See Modelica_LinearSystems2.Math.Matrices.LAPACK.dhseqr for details
</p>
</html>"));
end eigenvaluesHessenberg;

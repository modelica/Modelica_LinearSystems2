within Modelica_LinearSystems2.Math.Matrices.Internal;
function eigenvaluesHessenberg
  "Compute eigenvalues of an upper Hessenberg form matrix"
  import Modelica_LinearSystems2.Math.Matrices.Internal;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
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
  Integer lwork=Internal.dhseqr_workdim(H);
  Real work[lwork];
  Real Z[n,n]=fill(
        0,
        n,
        n);

  annotation (Documentation(info="<html>
This function uses DHSEQR Lapack-routine to calculate the eigenvalues of an upper Hessenberg form <b>H</b>.
Therefore, <b>H</b> is reduced to Schur form <b>T</b>. The eigenvalues are obtained from the diagonal of <b>T</b>.

<p>
See Modelica_LinearSystems2.Math.Matrices.LAPACK.dhseqr for details
</p>
</html>
"));
algorithm
  if size(H, 1) > 0 then
    (alphaReal,alphaImag,info) := LAPACK.dhseqr(H, lwork);
  else
    alphaReal := fill(0, size(H, 1));
    alphaImag := fill(0, size(H, 1));
  end if;

end eigenvaluesHessenberg;

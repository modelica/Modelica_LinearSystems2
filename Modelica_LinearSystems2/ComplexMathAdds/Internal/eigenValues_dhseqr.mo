within Modelica_LinearSystems2.ComplexMathAdds.Internal;
function eigenValues_dhseqr
  "Compute eingenvalues of a upper Hessenberg matrix using lapack routine DHSEQR"
  extends Modelica.Icons.Function;

  import
    Modelica_LinearSystems2.Math.Matrices.Internal.eigenvaluesHessenberg;

  input Real H[:,size(H, 1)] "Real upper Hessenberg matrix";
  output Complex zeros[size(H, 1)]
    "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

protected
  Integer nx=size(H, 1) "Number of states";
  Real alphaReal[nx];
  Real alphaImag[nx];
  Integer info;

algorithm
  (alphaReal,alphaImag,info) := eigenvaluesHessenberg(H);
  assert(info == 0,
    "Failed to compute eigenvalues with function Internal.eigenValues_dhseqr(..)");

  for i in 1:nx loop
    zeros[i].re := alphaReal[i];
    zeros[i].im := alphaImag[i];
  end for;

end eigenValues_dhseqr;

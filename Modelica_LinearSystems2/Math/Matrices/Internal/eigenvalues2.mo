within Modelica_LinearSystems2.Math.Matrices.Internal;
function eigenvalues2 "Compute eigenvalues and unnormalized eigenvectors"
  extends Modelica.Icons.Function;

  import Modelica.Math.Matrices.LAPACK;

  input Real A[:,size(A, 1)];

  output Real alphaReal[size(A, 1)]
    "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
  output Real alphaImag[size(A, 1)]
    "Imaginary part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
  output Real lEigenVectors[size(A, 1),size(A, 1)]
    "Left eigenvectors of matrix A";
  output Real rEigenVectors[size(A, 1),size(A, 1)]
    "Right eigenvectors of matrix A";
  output Integer info=0;
protected
  Integer n=size(A, 1);
  Real H[size(A, 1),size(A, 2)] "Upper Hessenberg form";
  Real V[size(A, 1),size(A, 2)]
    "V=[v1,v2,..vn-1,0] with vi are vectors which define the elementary reflectors";
  Real tau[size(A, 1) - 1];
  Real Q[size(A, 1),size(A, 2)];
  Real Zo[:,:];
  Real Ho[:,:];

algorithm
  (H,V,tau) := Modelica.Math.Matrices.Utilities.toUpperHessenberg(A,1,n);
  Q := LAPACK.dorghr(V, 1, n, tau);

  if size(H, 1) > 0 then
   (alphaReal,alphaImag,info,Ho,Zo) := LAPACK.dhseqr(H, false, "V", Q);
  else
    alphaReal := fill(0, size(H, 1));
    alphaImag := fill(0, size(H, 1));
    Zo := fill(0, size(H, 1), size(H, 2));
    Ho := fill(0, size(H, 1), size(H, 2));
  end if;

  (lEigenVectors,rEigenVectors) := LAPACK.dtrevc(Ho, "B", "B", Zo);

  annotation (Documentation(info="<html>
</html>"));
end eigenvalues2;

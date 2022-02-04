within Modelica_LinearSystems2.Math.Matrices;
function toUpperHessenberg
  "Transform a real general matrix A to upper Hessenberg form H by an orthogonal similarity transformation:  Q' * A * Q = H"

  input Real A[:,size(A, 1)] "Square matrix A";
  input Integer ilo=1
    "Lowest index where the original matrix had been Hessenbergform";
  input Integer ihi=size(A, 1)
    "Highest index where the original matrix had been Hessenbergform";
  output Real H[size(A, 1),size(A, 2)] "Upper Hessenberg form";
  output Real V[size(A, 1),size(A, 2)]
    "V=[v1,v2,..vn-1,0] with vi are vectors which define the elementary reflectors";

  output Real tau[max(0,size(A, 1) - 1)]
    "Scalar factors of the elementary reflectors";
  output Integer info;

protected
  Integer n=size(A, 1);
  Real Aout[size(A, 1),size(A, 2)];
  Integer i;

algorithm
  if n>0 then
  (Aout,tau,info) := Modelica.Math.Matrices.LAPACK.dgehrd(
    A,
    ilo,
    ihi);
  H[1:2, 1:ihi] := Aout[1:2, 1:ihi];
  H[1:2, ihi + 1:n] := A[1:2, ihi + 1:n];

  for i in 3:n loop
    H[i, i - 1:ihi] := Aout[i, i - 1:ihi];
    H[i, ihi + 1:n] := A[i, ihi + 1:n];
  end for;

  for i in 1:min(n - 2, ihi) loop
    V[i + 1, i] := 1.0;
    V[i + 2:n, i] := Aout[i + 2:n, i];

  end for;
  V[n, n - 1] := 1;
  end if;

  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.Utilities.toUpperHessenberg instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
                H = Matrices.<strong>toUpperHessenberg</strong>(A);
(H, V, tau, info) = Matrices.<strong>toUpperHessenberg</strong>(A,ilo, ihi);
</pre></blockquote>

<h4>Description</h4>
<p>
Function <strong>toUpperHessenberg</strong> computes a upper Hessenberg form <strong>H</strong>
of a matrix <strong>A</strong> by orthogonal similarity transformation:
<strong>Q</strong>' * <strong>A</strong> * <strong>Q</strong> = <strong>H</strong>. It calls LAPACK function DGEHRD.
See <a href=\"Modelica://Modelica.Math.Matrices.LAPACK.dgehrd\">LAPACK.dgehrd</a>
for more information about the additional outputs V, tau, info and
inputs ilo, ihi for more information.
</p>

<h4>Example</h4>
<blockquote><pre>
A  = [1, 2,  3;
      6, 5,  4;
      1, 0,  0];

H = toUpperHessenberg(A);
</pre></blockquote>
<p>
results in:
</p>
<blockquote><pre>
H = [1.0,  -2.466,  2.630;
    -6.083, 5.514, -3.081;
     0.0,   0.919, -0.514]
</pre></blockquote>
</html>"));
end toUpperHessenberg;

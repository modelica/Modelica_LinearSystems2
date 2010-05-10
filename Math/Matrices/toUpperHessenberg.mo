within Modelica_LinearSystems2.Math.Matrices;
function toUpperHessenberg
  "transform a real general matrix A to upper Hessenberg form H by an orthogonal similarity transformation:  Q' * A * Q = H"
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,size(A, 1)] "Square matrix A";
  input Integer ilo=1
    "lowest index where the original matrix had been Hessenbergform";
  input Integer ihi=size(A, 1)
    "highest index where the original matrix had been Hessenbergform";
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
  (Aout,tau,info) := LAPACK.dgehrd(
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

  annotation (Documentation(info="<html>
  
   <h4>Syntax</h4>
<blockquote><pre>
         H = Matrices.Utilities<b>toUpperHessenberg</b>(A);
         (H, V, tau, info) = Matrices.Utilities<b>toUpperHessenberg</b>(A,ilo, ihi);
</pre></blockquote>
<h4>Description</h4>
Function <b>toUpperHessenberg</b> computes a upper Hessenberg form <b>H</b> of a matrix <b>A</b> by orthogonal similarity transformation:  <b>Q</b>' * <b>A</b> * <b>Q</b> = <b>H</b>.
It calls LAPACK function DGEHRD. See Matrices.Lapack.dgehrd for more information about the additional outputs V, tau, info and inputs ilo, ihi for more information.
<p>


<h4>Example</h4>
<blockquote><pre>
 A  = [1, 2,  3;
       6, 5,  4;
       1, 0,  0]; 

 H = toUpperHessenberg(A);

  results in:
  
 H = [1.0,  -2.466,  2.630;
     -6.083, 5.514, -3.081;
      0.0,   0.919, -0.514]
      
</pre></blockquote>

</html>"));
end toUpperHessenberg;

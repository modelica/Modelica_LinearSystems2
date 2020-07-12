within Modelica_LinearSystems2.ComplexMathAdds;
function eigenValues
  "Compute eingenvalues of a matrix A, using lapack routine dgeevx"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,size(A, 1)] "Real matrix";
  output Complex eigval[size(A, 1)]
    "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

protected
  Integer nx=size(A, 1) "Number of states";

  Real alphaReal[nx];
  Real alphaImag[nx];

  Integer info;

algorithm
 // Compute eigenvalues
  if size(A, 1) > 0 then
     (alphaReal,alphaImag,,,,info) := LAPACK.dgeevx(A);
     assert(info == 0,
       "Failed to compute eigenvalues with function eigenValues_dgeevx(..)");

     for i in 1:nx loop
       eigval[i].re := alphaReal[i];
       eigval[i].im := alphaImag[i];
     end for;
  end if;
  annotation (Documentation(info="<html>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
z that make the following matrix singular:
</p>
<pre>
    | <b>A</b> <b>B</b> |     | <b>I</b> <b>0</b> |
    |     | - z*|     |
    | <b>C</b> <b>D</b> |     | <b>0</b> <b>0</b> |
</pre>
<p>
where <b>I</b> is the identity matrix of the same size as <b>A</b>
and <b>0</b> are zero matrices of appropriate dimensions.
</p>
<p>
Currently, there is the restriction that the number of
inputs and the number of outputs must be identical.
</p>
</html>"));
end eigenValues;

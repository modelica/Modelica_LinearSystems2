within Modelica_LinearSystems2.ComplexMathAdds;
function eigenValues
  "Compute eingenvalues of a matrix A, using lapack routine dgeevx"
  extends Modelica.Icons.Function;

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
     (alphaReal,alphaImag,,,,info) := Modelica.Math.Matrices.LAPACK.dgeevx(A);
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
   der(<strong>x</strong>) = <strong>A</strong>*<strong>x</strong> + <strong>B</strong>*<strong>u</strong>
        <strong>y</strong> = <strong>C</strong>*<strong>x</strong> + <strong>D</strong>*<strong>u</strong>
</pre>
<p>
The invariant zeros of this system are defined as the variables
z that make the following matrix singular:
</p>
<pre>
    | <strong>A</strong> <strong>B</strong> |     | <strong>I</strong> <strong>0</strong> |
    |     | - z*|     |
    | <strong>C</strong> <strong>D</strong> |     | <strong>0</strong> <strong>0</strong> |
</pre>
<p>
where <strong>I</strong> is the identity matrix of the same size as <strong>A</strong>
and <strong>0</strong> are zero matrices of appropriate dimensions.
</p>
<p>
Currently, there is the restriction that the number of
inputs and the number of outputs must be identical.
</p>
</html>"));
end eigenValues;

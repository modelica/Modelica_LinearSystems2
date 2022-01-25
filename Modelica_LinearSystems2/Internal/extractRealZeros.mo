within Modelica_LinearSystems2.Internal;
function extractRealZeros
  "Extract real and conjugate complex elements from a complex vector"

  input Complex complexVector[:]
    "Zeros of a polynomial with real coefficients, e.g., eigen values of a real matrix";
  input Integer numberOfRealZeros
    "Number of real zeros of ComplexVector determined with function numberOfRealZeros";
  input String name="complexVector"
    "Name of complexVector to be used in error message";
  output Real realZeros[numberOfRealZeros] "Real zeros of complexVector";
  output Complex complexZeros[:]=fill(Complex(0), integer((size(complexVector, 1) -
      numberOfRealZeros)/2))
    "Complex zeros without the corresponding conjugate complex pair element";

protected
  Integer n=size(complexVector, 1);
  Integer i;
  Integer i1;
  Integer i2;
algorithm
  i1 := 0;
  i2 := 0;
  i := 1;
  while i <= n loop
    if complexVector[i].im == 0.0 then
      i1 := i1 + 1;
      realZeros[i1] := complexVector[i].re;
      i := i + 1;
    else
        // check that the next two zeros are a conjugate complex pair
      assert(i < n, "Argument " + name +
        " does not define a real valued polynomial\n" + name + "[" + String(n)
         + "] is complex without complex conjugate.");

      assert(complexVector[i].re == complexVector[i + 1].re,
        "No conjugate complex pair (checked real parts)\n" + "  " + name + "["
         + String(i) + "] = " + String(complexVector[i]) + "\n" + "  " + name
         + "[" + String(i + 1) + "] = " + String(complexVector[i]) + "\n" +
        "and the real parts of these two complex numbers should be identical\n"
         + "since conjugate complex pairs required.");
      assert(complexVector[i].im == -complexVector[i + 1].im,
        "No conjugate complex pair (checked the imaginary parts)\n" + "  " +
        name + "[" + String(i) + "] = " + String(complexVector[i]) + "\n" +
        "  " + name + "[" + String(i + 1) + "] = " + String(complexVector[i + 1])
         + "\n" +
        "and the imaginary parts of these two complex numbers should be identical\n"
         + "with opposite sign, since a conjugate complex pair is required.");
      i2 := i2 + 1;
        // Store the zero with the positive imaginary part
      if complexVector[i].im >= 0 then
        complexZeros[i2] := complexVector[i];
      else
        complexZeros[i2] := complexVector[i + 1];
      end if;
      i := i + 2;
    end if;
  end while;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
                realZeros = Matrices.<strong>extractRealZeros</strong>(complexVector, numberOfRealRoots);
(realZeros, complexZeros) = Matrices.<strong>extractRealZeros</strong>(
                              complexVector,
                              numberOfRealRoots,
                              name=&quot;complexVector&quot;);
</pre></blockquote>

<h4>Description</h4>
<p>
This function extracts the real zeros from the Complex vector
&quot;complexVector&quot;. It is required that all elements of complexVector
define either a real zero (complexVector[i].im=0) or a conjugate complex zero
pair
(complexVector[i].re == complexVector[i+1].re and
complexVector[i].im == -complexVector[i+1].im).
The second argument &quot;numberOfRealZeros&quot; is determined by a function
call of Internal.numberOfRealZeros(). The optional input argument
&quot;name&quot; is used as name of &quot;complexVector&quot; in error messages.
</p>
<p>
The function returns the real elements of complexVector in
vector &quot;realZeros&quot; and the real and imaginary part of a conjugate
complex zero pair in matrix &quot;complexZeros[:]&quot;.

<h4>Example</h4>
<blockquote><pre>
  // c = {0, 1+2*j, 1-2*j, 2, -3, -1-1*j, -1+1*j};
  Integer n = numberOfRealZeros(c);
  Real realZeros[n];
  Real complexZeros[:] = fill(Complex(0), integer((size(c,1)-n)/2));
algorithm
  (realZeros, complexZeros) := extractRealZeros(c, n);
           -> realZeros    = {0, 2, (-3)};
              complexZeros = { 1+2*j,
                              -1+1*j}
</pre></blockquote>
</html>"));
end extractRealZeros;

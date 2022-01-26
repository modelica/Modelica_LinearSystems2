within Modelica_LinearSystems2.Internal;
function reorderEigenvalue
  "Extract real and conjugate complex elements from a complex vector as part of a record Eigenvalue"
  import Modelica_LinearSystems2.Internal.Eigenvalue;
  input Eigenvalue EigenvalueVector[:]
    "Zeros of a polynomial with real coefficients, e.g., eigen values of a real matrix";
  input String name="EigenvalueVector"
    "Name of complexVector to be used in error message";
  output Eigenvalue reorderedEigenvalues[size(EigenvalueVector, 1)]
    "Reordered zeros";
  output Integer nRealEigenvalues
    "Number of real zeros (EigenvalueVector[1:nRealEigenvalues] are the real zeros)";

protected
  Integer n=size(EigenvalueVector, 1);
  Integer i;
  Integer jr;
  Integer jc;
  Complex complexVector[size(EigenvalueVector, 1)];
algorithm
  nRealEigenvalues:=numberOfRealZeros(complexVector);

  i := 1;
  jr := 1;
  jc := nRealEigenvalues + 1;
  while i <= n loop
    if EigenvalueVector[i].ev.im == 0.0 then

      reorderedEigenvalues[jr] := EigenvalueVector[i];
      reorderedEigenvalues[jr].ev.im := 0.0;
      i := i + 1;
      jr := jr + 1;
    else
// check that the next two zeros are a conjugate complex pair
      assert(i < n, "Argument " + name +
        " does not define a real valued polynomial\n" + name + "[" + String(n)
         + "] is complex without complex conjugate.");
      assert(abs(EigenvalueVector[i].ev.re - EigenvalueVector[i+1].ev.re) < max(Modelica.Constants.eps,
        abs(EigenvalueVector[i+1].ev.re)*100*Modelica.Constants.eps),
        "No conjugate complex pair (checked the real parts)\n" + "  " + name +
        "[" + String(i) + "] = " + String(EigenvalueVector[i].ev) + "\n" + "  " +
        name + "[" + String(i + 1) + "] = " + String(EigenvalueVector[i].ev) + "\n"
         +
        "and the real parts of these two complex numbers should be identical\n"
         + "since conjugate complex pairs required.");
      assert(abs(EigenvalueVector[i].ev.im + EigenvalueVector[i+1].ev.im) < max(Modelica.Constants.eps,
        abs(EigenvalueVector[i+1].ev.im)*100*Modelica.Constants.eps),
        "No conjugate complex pair (checked the imaginary parts)\n" + "  " +
        name + "[" + String(i) + "] = " + String(EigenvalueVector[i].ev) + "\n" +
        "  " + name + "[" + String(i + 1) + "] = " + String(EigenvalueVector[i + 1].ev)
         + "\n" +
        "and the imaginary parts of these two complex numbers should be identical\n"
         + "with opposite sign, since a conjugate complex pair is required.");

      // Store the zero with the positive imaginary part

      if EigenvalueVector[i].ev.im >= 0 then
        reorderedEigenvalues[jc] := EigenvalueVector[i];
        reorderedEigenvalues[jc + 1] := EigenvalueVector[i+1];
      else
        reorderedEigenvalues[jc] := EigenvalueVector[i+1];
        reorderedEigenvalues[jc + 1] := EigenvalueVector[i];
      end if;
      i := i + 2;
      jc := jc + 2;
    end if;
  end while;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
               reorderedZero = Matrices.<strong>reorderZeros</strong>(complexVector);
(reorderedZeros, nRealZeros) = Matrices.<strong>reorderZeros</strong>(
                                 complexVector,
                                 name=&quot;complexVector&quot;);
</pre></blockquote>

<h4>Description</h4>
<p>
Function <strong>reorderZeros</strong>(..) reorders the zeros from the
Complex vector &quot;complexVector&quot; such that the returned Complex vector
reorderedZeros contains first all real Zeros and afterwards the conjugate
complex zero pairs. It is required that all elements
of complexVector define either a real zero (complexVector[i].im=0)
or a conjugate complex zero pair
(complexVector[i].re == complexVector[i+1].re and
complexVector[i].im == -complexVector[i+1].im).
The optional input argument
&quot;name&quot; is used as name of &quot;complexVector&quot; in error messages.
</p>
<p>
The function returns the vector element reordered, as well as
the number of real zeros (nRealZeros).
</p>

<h4>Example</h4>
<blockquote><pre>
  // c = {0, 1+2*j, 1-2*j, 2, -3, -1-1*j, -1+1*j};
  Real complexZeros[:] = fill(Complex(0), integer((size(c,1)-n)/2));
<strong>algorithm</strong>
  (reorderedZeros, nRealZeros) := reorderZeros(c);
  // reorderedZeros = {0, 2, (-3), 1+2*j, 1-2*j, -1+1*j, -1-1*j}
  // nRealZeros     = 3
</pre></blockquote>
</html>"));
end reorderEigenvalue;

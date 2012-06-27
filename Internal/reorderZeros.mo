within Modelica_LinearSystems2.Internal;
function reorderZeros
  "Extract real and conjugate complex elements from a complex vector"
  import Modelica_LinearSystems2.Math.Complex;

  input Complex complexVector[:]
    "Zeros of a polynomial with real coefficients, e.g., eigen values of a real matrix";
  input String name="complexVector"
    "Name of complexVector to be used in error message";
  output Complex reorderedZeros[size(complexVector, 1)] "Reordered zeros";
  output Integer nRealZeros=numberOfRealZeros(complexVector)
    "Number of real zeros (reorderedZeros[1:nRealZeros] are the real zeros)";

protected
  Integer n=size(complexVector, 1);
  Integer i;
  Integer jr;
  Integer jc;
algorithm
  i := 1;
  jr := 1;
  jc := nRealZeros + 1;
  while i <= n loop
    if abs(complexVector[i].im) < Modelica.Constants.eps then
      reorderedZeros[jr].re := complexVector[i].re;
      reorderedZeros[jr].im := 0.0;
      i := i + 1;
      jr := jr + 1;
    else
// check that the next two zeros are a conjugate complex pair
      assert(i < n, "Argument " + name +
        " does not define a real valued polynomial\n" + name + "[" + String(n)
         + "] is complex without complex conjugate.");
      assert(abs(complexVector[i].re - complexVector[i + 1].re) < max(Modelica.Constants.eps,
        abs(complexVector[i + 1].re)*100*Modelica.Constants.eps),
        "No conjugate complex pair (checked the real parts)\n" + "  " + name +
        "[" + String(i) + "] = " + String(complexVector[i]) + "\n" + "  " +
        name + "[" + String(i + 1) + "] = " + String(complexVector[i]) + "\n"
         +
        "and the real parts of these two complex numbers should be identical\n"
         + "since conjugate complex pairs required.");
      assert(abs(complexVector[i].im + complexVector[i + 1].im) < max(Modelica.Constants.eps,
        abs(complexVector[i + 1].im)*100*Modelica.Constants.eps),
        "No conjugate complex pair (checked the imaginary parts)\n" + "  " +
        name + "[" + String(i) + "] = " + String(complexVector[i]) + "\n" +
        "  " + name + "[" + String(i + 1) + "] = " + String(complexVector[i + 1])
         + "\n" +
        "and the imaginary parts of these two complex numbers should be identical\n"
         + "with opposite sign, since a conjugate complex pair is required.");

        // Store the zero with the positive imaginary part
      if complexVector[i].im >= 0 then
        reorderedZeros[jc] := complexVector[i];
        reorderedZeros[jc + 1] := complexVector[i + 1];
      else
        reorderedZeros[jc] := complexVector[i + 1];
        reorderedZeros[jc + 1] := complexVector[i];
      end if;
      i := i + 2;
      jc := jc + 2;
    end if;
  end while;
  annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
               reorderedZero = Internal.<b>reorderZeros</b>(complexVector);
(reorderedZeros, nRealZeros) = Internal.<b>reorderZeros</b>(complexVector, 
                                                     name=\"complexVector\");
</pre></blockquote>
<h4>Description</h4>
<p>
Function <b>reorderZeros</b>(..) reorders the zeros from the
Complex vector \"complexVector\" such that the returned Complex vector
reorderedZeros contains first all real Zeros and afterwards the conjugate
complex zero pairs. It is required that all elements
of complexVector define either a real zero (complexVector[i].im=0) 
or a conjugate complex zero pair
(complexVector[i].re == complexVector[i+1].re and
complexVector[i].im == -complexVector[i+1].im). 
The optional input argument
\"name\" is used as name of \"complexVector\" in error messages.
</p>
<p>
The function returns the vector element reordered, as well as
the number of real zeros (nRealZeros).
<h4>Example</h4>
<blockquote><pre>
    
  // c = {0; 1+2j; 1-2j; 2; -3; -1-j; -1+j};
    Real complexZeros[:] = fill(Complex(0), integer((size(c,1)-n)/2));
  algorithm
  (reorderedZeros, nRealZeros) := reorderZeros(c);
      -> reorderedZeros = {0, 2, (-3), 1+2j, 1-2j, -1+j, -1-j}
         nRealZeros     = 3
 
</pre></blockquote>
</HTML>"));
end reorderZeros;

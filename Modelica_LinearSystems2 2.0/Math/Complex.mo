within Modelica_LinearSystems2.Math;
record Complex "Record defining a Complex number"
  extends Modelica.Icons.Record;
  Real re "Real part of complex number" annotation(Dialog);
  Real im "Imaginary part of complex number" annotation(Dialog);

  annotation (
    Documentation(info="<html>
<p>
This record defines a complex number consisting of a real
and an imaginary part. Additionally, various utility functions
are provided in the record to operate on instances of this record.
Example (note: \"j\" in the comments is defined as j=sqrt(-1)):
</p>
<pre>
   <b>import</b> Modelica.Utilities.Streams;
   <b>import</b> Modelica_LinearSystems2.Math.Complex;
   Complex c1=Complex(re=2, im=3) \"= 2 + 3j\";
   Complex c2=Complex(3,4) \"= 3 + 4j\";
 <b>algorithm</b>
   c3 := Complex.'+'(c1, c2) \"= c1 + c2\";
   Streams.print(\"c3 = \" + Complex.'String'(c3));
   Streams.print(\"c3 = \" + Complex.'String'(c3,\"i\"));
   // This gives the following print-out:
   c3 = 5 + 7j
   c3 = 5 + 7i
</pre>
<p>
The utility functions are written in such a way that it
is convenient to work with them, once operator overloading
is provided in Modelica and Modelica tools. Example:
</p>
<pre>
   // Assume operator overloading is available (this is not yet the case):
   Complex j  = Complex.j();
   Complex c4 = -2 + 5*j;
   // A Modelica tool will transform the statement of c4 into
   Complex c4 = Complex.'+'( Complex.(-2,0), Complex.'*'(Complex(5,0),j)));
</pre>
<p>
The utility functions are implemented so that they can be easily
inlined by a tool. In such a case, the above statement will
not lead to any overhead.
</p>
</html>"));

 encapsulated package Examples
    "Library demonstrating the usage of complex numbers"

    import Modelica;
    extends Modelica.Icons.Library;

    function addTwoComplexNumbers "Show how to add 2 complex number"
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica.Utilities.Streams;

      output Boolean ok;
    protected
      Complex j = Complex.j();
      Complex c1=2+3*j;
      Complex c2=3+4*j;
      Complex c3;
    algorithm

      Streams.print("c1 = " + String(c1));
      Streams.print("c2 = " + String(c2));

      c3 := c1 + c2;
      Streams.print("c3 = c1 + c2 = " + String(c3));
      ok := true;
    end addTwoComplexNumbers;

 end Examples;

encapsulated package Vectors
 function print
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2.Math.Complex;

  input String name="" "Name of complex vector";
  input Complex c[:] "Complex vector to be printed";

  input String fileName=""
        "Print to fileName; empty fileName prints to the log window";

 algorithm
  print(name + " =", fileName);
  for i in 1:size(c, 1) loop
     print("   " + String(c[i]), fileName);
  end for;
 end print;

function length "Return length of a complex vector"
      import Modelica_LinearSystems2.Math.Complex;
  input Complex v[:] "Vector";
  output Real result "Length of vector v";

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
Vectors.<b>length</b>(v);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Vectors.<b>length</b>(v)</code>\" returns the
<b>Euclidean length</b> \"<code>sqrt(v*v)</code>\" of vector v.
The function call is equivalent to Vectors.norm(v). The advantage of
length(v) over norm(v)\"is that function length(..) is implemented
in one statement and therefore the function is usually automatically
inlined. Further symbolic processing is therefore possible, which is
not the case with function norm(..).
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  v = {2, -4, -2, -1};
  <b>length</b>(v);  // = 5
</pre></blockquote>
<h4><font color=\"#008000\">See also</font></h4>
<a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.norm\">Vectors.norm</a>
</html>"));
algorithm
  result := sqrt(sum(v[i].re^2 + v[i].im^2 for i in 1:size(v,1)));
end length;

function norm "Returns the norm of a complex vector"
      import Modelica;
      import Modelica_LinearSystems2.Math.Complex;
  input Complex v[:] "Vector";
  input Real p(min=1) = 2
        "Type of p-norm (often used: 1, 2, or Modelica.Constants.inf)";
  output Real result "p-norm of vector v";

  annotation (Documentation(info="<HTML>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
Vectors.<b>norm</b>(v);
Vectors.<b>norm</b>(v,p=2);   // 1 &le; p &le; &#8734;
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Vectors.<b>norm</b>(v)</code>\" returns the
<b>Euclidean norm</b> \"<code>sqrt(v*v)</code>\" of vector v.
With the optional
second argument \"p\", any other p-norm can be computed:
</p>
<center>
<IMG SRC=\"../Extras/Images/vectorNorm.png\" ALT=\"function Vectors.norm\">
</center>
<p>
Besides the Euclidean norm (p=2), also the 1-norm and the
infinity-norm are sometimes used:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>1-norm</b></td>
      <td>= sum(abs(v))</td>
      <td><b>norm</b>(v,1)</td>
  </tr>
  <tr><td><b>2-norm</b></td>
      <td>= sqrt(v*v)</td>
      <td><b>norm</b>(v) or <b>norm</b>(v,2)</td>
  </tr>
  <tr><td><b>infinity-norm</b></td>
      <td>= max(abs(v))</td>
      <td><b>norm</b>(v,Modelica.Constants.<b>inf</b>)</td>
  </tr>
</table>
<p>
Note, for any vector norm the following inequality holds:
</p>
<blockquote><pre>
<b>norm</b>(v1+v2,p) &le; <b>norm</b>(v1,p) + <b>norm</b>(v2,p)
</pre></blockquote>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  v = {2, -4, -2, -1};
  <b>norm</b>(v,1);    // = 9
  <b>norm</b>(v,2);    // = 5
  <b>norm</b>(v);      // = 5
  <b>norm</b>(v,10.5); // = 4.00052597412635
  <b>norm</b>(v,Modelica.Constants.inf);  // = 4
</pre></blockquote>
<h4><font color=\"#008000\">See also</font></h4>
<a href=\"Modelica:Modelica.Math.Matrices.norm\">Matrices.norm</a>
</HTML>"));
algorithm
  if p == 2 then
    result:= sqrt(sum(v[i].re^2 + v[i].im^2 for i in 1:size(v,1)));
  elseif p == Modelica.Constants.inf then
    result:= Complex.'max'(v);
  elseif p == 1 then
    result:= sum(Complex.'abs'(v[i]) for i in 1:size(v,1));
  else
    result:=(sum(Complex.'abs'(v[i])^p for i in 1:size(v, 1)))^(1/p);
  end if;
end norm;

function normalize
      "Return normalized complex vector such that length = 1 and prevent zero-division for zero vector"
      import Modelica;
      import Modelica_LinearSystems2.Math.Complex;

  input Complex v[:] "Vector";
  input Real eps = 100*Modelica.Constants.eps "if |v| < eps then result = v";
  output Complex result[size(v, 1)] "Input vector v normalized to length=1";

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
Vectors.<b>normalize</b>(v);
Vectors.<b>normalize</b>(v,eps=100*Modelica.Constants.eps);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Vectors.<b>normalize</b>(v)</code>\" returns the
<b>unit vector</b> \"<code>v/length(v)</code>\" of vector v.
If length(v) is close to zero (more precisely, if length(v) &lt; eps),
v is returned in order to avoid
a division by zero. For many applications this is useful, because
often the unit vector <b>e</b> = <b>v</b>/length(<b>v</b>) is used to compute
a vector x*<b>e</b>, where the scalar x is in the order of length(<b>v</b>),
i.e., x*<b>e</b> is small, when length(<b>v</b>) is small and then
it is fine to replace <b>e</b> by <b>v</b> to avoid a division by zero.
</p>
<p>
Since the function is implemented in one statement,
it is usually inlined and therefore symbolic processing is
possible.
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  <b>normalize</b>({1,2,3});  // = {0.267, 0.534, 0.802}
  <b>normalize</b>({0,0,0});  // = {0,0,0}
</pre></blockquote>
<h4><font color=\"#008000\">See also</font></h4>
<a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.length\">Vectors.length</a>
</html>"));
    protected
  Real length_v = Complex.Vectors.length(v);
  Complex j = Complex.j();
algorithm
  if length_v >= eps then
     for i in 1:size(v,1) loop
         result[i] :=v[i].re/length_v + (v[i].im/length_v)*j;
     end for;
  else
     result :=v;
  end if;
end normalize;

function sortComplex "Sort elements of complex vector"
      import Modelica_LinearSystems2.Math.Complex;
  input Complex v[:] "Vector to be sorted";
  input Boolean ascending = true
        "= true if ascending order, otherwise descending order";
  input Boolean sortFrequency=true
        "= true, if sorting is first for imaginary then for real value; = false, if sorting is for absolute value";
  output Complex sorted_v[size(v,1)] = v "Sorted vector";
  output Integer indices[size(v,1)] = 1:size(v,1) "sorted_v = v[indices]";

  annotation (Documentation(info="<HTML>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
           sorted_v = Vectors.<b>sort</b>(v);
(sorted_v, indices) = Vectors.<b>sort</b>(v, ascending=true);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>sort</b>(..) sorts a Real vector v
in ascending order and returns the result in sorted_v.
If the optional argument \"ascending\" is <b>false</b>, the vector
is sorted in descending order. In the optional second
output argument the indices of the sorted vector with respect
to the original vector are given, such that sorted_v = v[indices].
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  (v2, i2) := Vectors.sort({-1, 8, 3, 6, 2});
       -> v2 = {-1, 2, 3, 6, 8}
          i2 = {1, 5, 3, 4, 2}
</pre></blockquote>
</HTML>"));
  /* shellsort algorithm; should be improved later */
    protected
  Integer gap;
  Integer i;
  Integer j;
  Complex wv;
  Integer wi;
  Integer nv = size(v,1);
  Boolean swap;
  Integer k1;
  Integer k2;
algorithm
  gap := div(nv,2);

  while gap > 0 loop
     i := gap;
     while i < nv loop
        j := i-gap;
        if j>=0 then
           k1 := j+1;
           k2 := j + gap + 1;
           if sortFrequency then
              if ascending then
                 swap := abs(sorted_v[k1].im) >  abs(sorted_v[k2].im) or
                         abs(sorted_v[k1].im) == abs(sorted_v[k2].im) and
                         (sorted_v[k1].re  > sorted_v[k2].re or
                          sorted_v[k1].re  == sorted_v[k2].re and sorted_v[k1].im < sorted_v[k2].im);
              else
                 swap := abs(sorted_v[k1].im) <  abs(sorted_v[k2].im) or
                         abs(sorted_v[k1].im) == abs(sorted_v[k2].im) and
                         (sorted_v[k1].re  < sorted_v[k2].re or
                          sorted_v[k1].re  == sorted_v[k2].re and sorted_v[k1].im < sorted_v[k2].im);
              end if;
           else
              if ascending then
                 swap := Complex.'abs'(sorted_v[k1]) > Complex.'abs'(sorted_v[k2]);
              else
                 swap := Complex.'abs'(sorted_v[k1]) < Complex.'abs'(sorted_v[k2]);
              end if;
           end if;
        else
           swap := false;
        end if;

        while swap loop
           wv := sorted_v[j+1];
           wi := indices[j+1];
           sorted_v[j+1] := sorted_v[j+gap+1];
           sorted_v[j+gap+1] := wv;
           indices[j+1] := indices[j+gap+1];
           indices[j+gap+1] := wi;
           j := j - gap;
           if j >= 0 then
              k1 := j+1;
              k2 := j + gap + 1;
              if sortFrequency then
                 if ascending then
                    swap := abs(sorted_v[k1].im) >  abs(sorted_v[k2].im) or
                            abs(sorted_v[k1].im) == abs(sorted_v[k2].im) and
                            (sorted_v[k1].re  > sorted_v[k2].re or
                             sorted_v[k1].re  == sorted_v[k2].re and sorted_v[k1].im < sorted_v[k2].im);
                 else
                    swap := abs(sorted_v[k1].im) <  abs(sorted_v[k2].im) or
                            abs(sorted_v[k1].im) == abs(sorted_v[k2].im) and
                            (sorted_v[k1].re  < sorted_v[k2].re or
                             sorted_v[k1].re  == sorted_v[k2].re and sorted_v[k1].im < sorted_v[k2].im);
                 end if;
              else
                 if ascending then
                    swap := Complex.'abs'(sorted_v[k1]) > Complex.'abs'(sorted_v[k2]);
                 else
                    swap := Complex.'abs'(sorted_v[k1]) < Complex.'abs'(sorted_v[k2]);
                 end if;
              end if;
           else
              swap := false;
           end if;
        end while;
        i := i + 1;
     end while;
     gap := div(gap,2);
  end while;
end sortComplex;
end Vectors;

  encapsulated operator 'constructor'
    function fromReal
      import Modelica_LinearSystems2.Math.Complex;
      input Real re "Real part of complex number";
      input Real im=0 "Imaginary part of complex number";
      output Complex result "Complex number";
    algorithm
      result.re :=re;
      result.im :=im;
    end fromReal;
  end 'constructor';

  encapsulated operator '-' "Unary and binary minus"
    function negate "Unary minus (multiply complex number by -1)"
      import Modelica_LinearSystems2.Math.Complex; // changed to Modelica_LinearSystems2
      input Complex c1 "Complex number";
      output Complex c2 "= -c1";
    algorithm
      c2 := Complex(-c1.re, -c1.im);
    end negate;

    function subtract "Subtract two complex numbers"
      import Modelica_LinearSystems2.Math.Complex;// changed to Modelica_LinearSystems2
      input Complex c1 "Complex number 1";
      input Complex c2 "Complex number 2";
      output Complex c3 "= c1 - c2";
    algorithm
      c3 := Complex(c1.re - c2.re, c1.im - c2.im);
    end subtract;
  end '-';

  encapsulated operator function '+' "Add two complex numbers"
    import Modelica_LinearSystems2.Math.Complex;

      input Complex c1 "Complex number 1";
      input Complex c2 "Complex number 2";
      output Complex c3 "= c1 + c2";

  algorithm
    c3 := Complex(c1.re + c2.re, c1.im + c2.im);
    //    end add;
  end '+';

  encapsulated operator function '*' "Multiply two complex numbers"
    import Modelica_LinearSystems2.Math.Complex;

      input Complex c1 "Complex number 1";
      input Complex c2 "Complex number 2";
      output Complex c3 "= c1*c2";
  algorithm
      c3 := Complex(c1.re*c2.re - c1.im*c2.im, c1.re*c2.im + c1.im*c2.re);
    //    end multiply;
  end '*';

  encapsulated operator function '/' "Divide two complex numbers"
    import Modelica_LinearSystems2.Math.Complex;

      input Complex c1 "Complex number 1";
      input Complex c2 "Complex number 2";
      output Complex c3 "= c1/c2";
  algorithm
      c3 := Complex((c1.re*c2.re + c1.im*c2.im)/(c2.re^2 + c2.im^2), (-c1.re*c2.im + c1.im*c2.re)/(c2.re^2 + c2.im^2));
    //    end divide;
  end '/';

  encapsulated operator function '=='
    "Test whether two complex numbers are identical"
    import Modelica_LinearSystems2.Math.Complex;
      input Complex c1 "Complex number 1";
      input Complex c2 "Complex number 2";
      output Boolean result "c1 == c2";
  algorithm
      result := c1.re == c2.re and c1.im == c2.im;
    //    end equals;
  end '==';

  encapsulated operator function 'String'
    "Transform Complex number into a String representation"
    import Modelica_LinearSystems2.Math.Complex;
      input Complex c
      "Complex number to be transformed in a String representation";
      input String name="j"
      "Name of variable representing sqrt(-1) in the string";
      input Integer significantDigits=6
      "Number of significant digits that are shown";
      output String s="";
  algorithm
      s := String(c.re, significantDigits=significantDigits);
      if c.im <> 0 then
        if c.im > 0 then
          s := s + " + ";
        else
          s := s + " - ";
        end if;
        s := s + String(abs(c.im), significantDigits=significantDigits) + "*" + name;
      end if;
    //    end toString;
  end 'String';

  encapsulated function j "Returns sqrt(-1)"
    import Modelica_LinearSystems2.Math.Complex;

    output Complex c "= sqrt(-1)";
  algorithm
    c := Complex(0,1);
  end j;

  encapsulated function 'abs' "Absolute value of complex number"
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c "Complex number";
    output Real result "= abs(c)";
  algorithm
    result := (c.re^2 + c.im^2)^0.5; //changed from sqrt
  end 'abs';

  encapsulated function 'sqrt' "Square root of complex number"
    import Modelica.Math;
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c1 "Complex number";
    output Complex c2 "= sqrt(c1)";
  algorithm
    c2 := Complex(sqrt(Complex.'abs'(c1))*Math.cos(Complex.arg(c1)/2), sqrt(
      Complex.'abs'(c1))*Math.sin(Complex.arg(c1)/2));
  end 'sqrt';

encapsulated function 'max' "Return maximum element of complex vector"
    import Modelica_LinearSystems2.Math.Complex;

  input Complex v[:] "Vector";
  output Real result "Element of v with largest absolute value";
  output Integer index "v[index] has the largest absolute value";

  annotation (Documentation(info="<html>

</html>"));
  protected
  Real absv_i;
algorithm
  if size(v,1) > 0 then
     result := Complex.'abs'(v[1]);
     index  := 1;
     for i in 2:size(v,1) loop
        absv_i :=Complex.'abs'(v[i]);
        if absv_i > result then
           result := absv_i;
           index := i;
        end if;
     end for;
  else
     result := 0;
     index  := 0;
  end if;
end 'max';

  encapsulated function exp "Exponential of complex number"
    import Modelica_LinearSystems2.Math.Complex;
    import Modelica.Math;

    input Complex c1 "Complex number";
    output Complex c2 "= exp(c1)";
  algorithm
     c2 := Complex(Math.exp(c1.re)*Math.cos(c1.im), Math.exp(c1.re)*Math.sin(c1.im));
  end exp;

  encapsulated function log "Logarithm of complex number"

    import Modelica;
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c1 "Complex number";
    output Complex c2 "= log(c1)";
  algorithm
    c2 := Complex(Modelica.Math.log(Complex.'abs'(c1)), Complex.arg(c1));
  end log;

  encapsulated function sin "Sine of complex number"

    import Modelica_LinearSystems2.Math.Complex;

    input Complex c1 "Complex number";
    output Complex c2 "sin(c1)";

  algorithm
     c2 := (Complex.exp(Complex(-c1.im, +c1.re)) - Complex.exp(Complex(+c1.im, -c1.re)))/
      Complex(0, 2);

  end sin;

  encapsulated function cos "Cosine of complex number"
    import Modelica_LinearSystems2.Math.Complex;
    import Modelica_LinearSystems2.Math.Complex.exp;

    input Complex c1 "Complex number";
    output Complex c2 "= cos(c1)";

  algorithm
   c2 := (exp(Complex(-c1.im, +c1.re)) + exp(Complex(+c1.im, -c1.re)))/2;
  end cos;

  encapsulated function arg "Phase angle of complex number"
    import Modelica;
    import Modelica_LinearSystems2.Math;
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c "Complex number";
    input Modelica.SIunits.Angle phi0=0
      "phase angle phi shall be in the range: -pi < phi-phi0 < pi";
    output Modelica.SIunits.Angle phi "= phase angle of c";

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
   Complex.<b>arg</b>(c);
   Complex.<b>arg</b>(c, phi0=0);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Complex.<b>arg</b>(c)</code>\" returns the
phase angle phi of the Complex number c in the range
-pi &lt; phi &lt; pi.<br>
The function call \"<code>Complex.<b>arg</b>(c,phi0)</code>\" returns the
phase angle phi of the Complex number c in the range
-pi &lt; phi - phi0 &lt; pi.
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  Complex.<b>arg</b>( Complex(1,0.5), 4*pi );  // = 4*pi+pi/4 = 13.351...
</pre></blockquote>
</html>"));
  algorithm
    phi := Modelica.Math.atan3(
        c.im,
        c.re,
        phi0);
  end arg;

  encapsulated function conj "Conjugate of complex number"
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c1 "Complex number";
    output Complex c2 "= c1.re - j*c1.im";
  algorithm
    c2 := Complex(c1.re, -c1.im);
  end conj;

  encapsulated function real "Real part of complex number"
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c "Complex number";
    output Real r "= c.re ";
  algorithm
    r := c.re;
  end real;

  encapsulated function imag "imaginary part of complex number"
    import Modelica_LinearSystems2.Math.Complex;

    input Complex c "Complex number";
    output Real r "= c.im ";
  algorithm
    r := c.im;
  end imag;

  encapsulated function eigenValues
    "compute eingenvalues of a matrix A, using lapack routine dgeevx"

    import Modelica_LinearSystems2.Math.Complex;
    import Modelica_LinearSystems2.Math.Matrices.LAPACK;
    import Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx;

    input Real A[:,size(A, 1)] "Real matrix";
    output Complex eigval[size(A, 1)]
      "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

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
  protected
    Integer nx=size(A, 1) "Number of states";

    Real alphaReal[nx];
    Real alphaImag[nx];

    Integer info;

  algorithm
   // Compute eigenvalues

    (alphaReal,alphaImag,,,,info) := dgeevx(A);
    assert(info == 0,
      "Failed to compute eigenvalues with function eigenValues_dgeevx(..)");

    for i in 1:nx loop
      eigval[i].re := alphaReal[i];
      eigval[i].im := alphaImag[i];
    end for;

  end eigenValues;

  encapsulated function frequency
    "Frequency and damping of conjugated complex pole pair"
    import Modelica;
    import Modelica_LinearSystems2.Math;
    import Modelica_LinearSystems2.Math.Complex;
    input Complex c "Complex number";
    output Modelica.SIunits.Frequency f "Frequency of c (= c.im in Hz)";
    output Real damping "Damping of c (= c.re/c.im)"
    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
   Complex.<b>arg</b>(c);
   Complex.<b>arg</b>(c, phi0=0);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Complex.<b>arg</b>(c)</code>\" returns the
phase angle phi of the Complex number c in the range
-pi &lt; phi &lt; pi.<br>
The function call \"<code>Complex.<b>arg</b>(c,phi0)</code>\" returns the
phase angle phi of the Complex number c in the range
-pi &lt; phi - phi0 &lt; pi.
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  Complex.<b>arg</b>( Complex(1,0.5), 4*pi );  // = 4*pi+pi/4 = 13.351...
</pre></blockquote>
</html>"));

  protected
    Real abs_ev=(c.re^2 + c.im^2)^0.5;
  algorithm
    f := if abs(c.im) > 10*Modelica.Constants.eps then abs_ev/(2*Modelica.Constants.pi) else 0;
    damping := if abs(c.im) > 10*Modelica.Constants.eps then if abs_ev > Modelica.Constants.eps then -c.re/abs_ev else 0.0 else
      1.0;
  end frequency;

encapsulated package Internal
 function eigenValues_dhseqr
      "compute eingenvalues of a upper Hessenberg matrix using lapack routine DHSEQR"

      import Modelica;
      import Modelica_LinearSystems2.Math.Complex;
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
end Internal;

end Complex;

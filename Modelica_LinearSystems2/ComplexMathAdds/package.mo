within Modelica_LinearSystems2;
package ComplexMathAdds "Library of additional complex mathematical functions and of functions operating on complex vectors and matrices"
  extends Modelica.Icons.Package;

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
end ComplexMathAdds;

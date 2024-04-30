within Modelica_LinearSystems2.ComplexMathAdds.Vectors;
function print "Convert a complex vector into a string and save it in a file"
  extends Modelica.Icons.Function;
  import Modelica.Utilities.Streams.print;

  input String name="" "Name of complex vector";
  input Complex c[:] "Complex vector to be printed";

  input String fileName=""
    "Print to fileName; empty fileName prints to the log window";

algorithm
  print(name + " =", fileName);
  for i in 1:size(c, 1) loop
    print("   " + String(c[i]), fileName);
  end for;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.<strong>print</strong>(name, c, fileName);
</pre></blockquote>

<h4>Description</h4>
<p>
Convert a&nbsp;complex vector into a&nbsp;string representation and save it line by line
in a&nbsp;file.
If <code>fileName</code> stays empty, the string will be printed to the log window.
</p>

<h4>Example</h4>
<blockquote><pre>
c1 = {Complex(1,3), Complex(2,7), Complex(-3,1)};

Vectors.print(\"c_vec\",c1)
// c_vec =
//    1 + 3*j
//    2 + 7*j
//    -3 + 1*j
</pre></blockquote>
</html>"));
end print;

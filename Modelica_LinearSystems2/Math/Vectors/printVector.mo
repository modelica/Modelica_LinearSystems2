within Modelica_LinearSystems2.Math.Vectors;
function printVector "Print vector"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica.Utilities.Strings;

  input Real v[:] "Vector of real numbers to be printed";
  input Integer significantDigits=6
    "Number of significant digits that are shown";
  input String name="v" "Independent variable name used for printing";
  output String s="" "String containing v";
protected
  String blanks=Strings.repeat(significantDigits);
  String space=Strings.repeat(8);
  String space2=Strings.repeat(3);
  Integer r=size(v, 1);

algorithm
  if r == 0 then
    s := name + " = []";
  else
    s := "\n" + name + " = \n";
    for i in 1:r loop
      s := s + space;

      if v[i] >= 0 then
        s := s + " ";
      end if;
      s := s + String(v[i], significantDigits=significantDigits) +
        Strings.repeat(significantDigits + 8 - Strings.length(String(abs(v[i]))));

      s := s + "\n";
    end for;

  end if;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
s = Vectors.<b>printVector</b>(v,significantDigits,name);
</pre></blockquote>

<h4>Description</h4>
<p>
This function returns string <code>s</code> containing vector of real numbers <code>v</code> called <code>name</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
v = {2,0,-342.1};
<b>printVector</b>(v, 6, &quot;myVec&quot;);
//  = &quot;
// myVec = 
//          2          
//          0          
//         -342      
// &quot;
</pre></blockquote>
</html>"));
end printVector;

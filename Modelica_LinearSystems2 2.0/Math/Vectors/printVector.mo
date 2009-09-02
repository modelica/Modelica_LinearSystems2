within Modelica_LinearSystems2.Math.Vectors;
function printVector "print vector"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica.Utilities.Strings;

  input Real v[:];
  input Integer significantDigits=6
    "Number of significant digits that are shown";
  input String name="v" "Independent variable name used for printing";
  output String s="";
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
end printVector;

within Modelica_LinearSystems2.ComplexMathAdds.Matrices;
function print "Print matrix"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Strings;

  input Complex M[:,:];
  input Integer significantDigits=6
    "Number of significant digits that are shown";
  input String name="M" "Independent variable name used for printing";
  output String s="";
protected
  String blanks=Strings.repeat(significantDigits);
  String space=Strings.repeat(8);
  String space2=Strings.repeat(3);
  Integer r=size(M, 1);
  Integer c=size(M, 2);

algorithm
  if r == 0 or c == 0 then
    s := name + " = []";
  else
    s := "\n" + name + " = \n";
    for i in 1:r loop
      s := s + space;
      for j in 1:c loop
        if Modelica.ComplexMath.abs(M[i, j]) >= 0 then
          s := s + " ";
        end if;
        s := s + String(M[i, j], significantDigits=significantDigits) +
          Strings.repeat(significantDigits + 8 - Strings.length(String(Modelica.ComplexMath.abs(M[i,j]))));

      end for;
      s := s + "\n";
    end for;

  end if;
end print;

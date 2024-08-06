within Modelica_LinearSystems2.Internal;
function printComplexVector "Print elements of vector of complex numbers"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;

  input String name "Name of vector c";
  input Complex c[:] "Vector of complex numbers";

algorithm
  print(name + " =");
  for i in 1:size(c, 1) loop
    print("   " + String(c[i]));
  end for;
end printComplexVector;

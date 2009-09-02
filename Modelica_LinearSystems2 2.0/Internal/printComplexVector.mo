within Modelica_LinearSystems2.Internal;
function printComplexVector
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Complex;

  input String name;
  input Complex c[:];

algorithm
  print(name + " =");
  for i in 1:size(c, 1) loop
    print("   " + String(c[i]));
  end for;
end printComplexVector;

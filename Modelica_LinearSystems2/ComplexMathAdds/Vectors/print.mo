within Modelica_LinearSystems2.ComplexMathAdds.Vectors;
function print "Print vector"
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
end print;

within Modelica_LinearSystems2.ComplexMathAdds.Vectors;
function multiply "Scalar product of two complex vectors"
  extends Modelica.Icons.Function;

  input Complex v1[:];
  input Complex v2[size(v1,1)];
  output Complex result=Complex(0);

algorithm
  for i in 1:size(v1,1) loop
    result := result + v1[i]*v2[i];
  end for;

end multiply;

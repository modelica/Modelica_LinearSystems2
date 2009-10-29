within Modelica_LinearSystems2.Tests.Design;
function data_Kautsky2 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[5,5],
  redeclare Real B[5,2],
  redeclare Complex assignedPoles[5]);

protected
 Complex j = Complex.j();

algorithm
  data.A:=[-0.1094,  0.0628, 0.0,     0.0,  0.0;
                1.306,  -2.132,   0.9807,  0.0,  0.0;
                0.0,  1.595,   -3.149,  1.547,  0.0;
                0.0,  0.0355,  2.632,  -4.257, 1.855;
                0.0,  0.00227,  0.0,  0.1636,  -0.1625];
  data.B:=[0.0,  0.0;  0.0638,  0.0;  0.0838,  -0.1396;  0.1004,  -0.206;  0.0063,  -0.0128];
  data.assignedPoles:=cat(1,Complex(1)*{-0.2,-0.5,-1.0},{-1.0+j, -1.0-j});
  data.K:=fill(0, 0, 0);
end data_Kautsky2;

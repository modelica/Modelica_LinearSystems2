within Modelica_LinearSystems2.Tests.Design;
function data_Byers3 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[4,4],
  redeclare Real B[4,2],
  redeclare Complex assignedPoles[4]);

algorithm
  data.A:=[-65.0,  65.0,  -19.5,  19.5;
                 0.1,  -0.1,    0.0,   0.0;
                 1.0,   0.0,   -0.5,  -1.0;
                 0.0,   0.0,    0.4,   0.0];
  data.B:=[65.0,  0.0; 0.0,  0.0; 0.0,  0.0;0.0, 0.4];
  data.assignedPoles:=Complex(1)*{-1.0, -2.0, -3.0, -4.0};
  data.K:=fill(0, 0, 0);
end data_Byers3;

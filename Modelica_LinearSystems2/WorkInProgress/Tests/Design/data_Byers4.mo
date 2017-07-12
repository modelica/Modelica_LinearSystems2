within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Byers4 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[3,3],
  redeclare Real B[3,2],
  redeclare Complex assignedPoles[3]);

algorithm
  data.A:=[0.0,   1.0,   0.0;
               0.0,   0.0,   1.0;
              -6.0, -11.0,  -6.0];
  data.B:=[1.0,  1.0;  0.0,  1.0; 1.0,  1.0];
  data.assignedPoles:=Complex(1)*{-1.0, -2.0, -3.0};
  data.K:=fill(0, 0, 0);
end data_Byers4;

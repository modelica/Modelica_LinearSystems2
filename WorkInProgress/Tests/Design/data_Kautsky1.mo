within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Kautsky1 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[4,4],
  redeclare Real B[4,2],
  redeclare Complex assignedPoles[4]);

algorithm
  data.A:=[1.38,-0.2077,6.715,-5.676; -0.5814,-4.29,0.0,0.675; 1.067,4.273,-6.654,
    5.893; 0.048,4.273,1.343,-2.104];
  data.B:=[0.0,0.0; 5.679,0.0; 1.136,-3.146; 1.136,0.0];
  data.assignedPoles:=Complex(1)*{-0.2,-0.5,-5.05657,-8.66589};
  data.K:=fill(0,0,0);
end data_Kautsky1;

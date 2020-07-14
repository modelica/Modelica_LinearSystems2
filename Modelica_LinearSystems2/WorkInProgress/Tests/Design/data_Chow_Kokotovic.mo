within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Chow_Kokotovic "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;
  input Real d=1e-6;
  output DesignData data(
  redeclare Real A[4,4],
  redeclare Real B[4,1],
  redeclare Complex assignedPoles[4]);

algorithm
  data.A:=[0, 0.4, 0, 0; 0, 0, 0.345, 0; 0, -0.524/d, -0.465/d, 0.262/d; 0, 0, 0, -1/d];
  data.B:=[0; 0; 0; 1/d];
  data.assignedPoles:=Complex(1)*{-1.0, -1.0, -3.0, -4.0};
  data.K:=fill(0, 0, 0);
end data_Chow_Kokotovic;

within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Byers6 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[4,4],
  redeclare Real B[4,2],
  redeclare Complex assignedPoles[4]);

algorithm
  data.A:=[5.8765,  9.3456,  4.5634,  9.3520;
                 6.6526,  0.5867,  3.5829,  0.6534;
                 0.0,     9.6738,  7.4876,  4.7654;
                 0.0,     0.0,     6.6784,  2.5678];
  data.B:=[3.9878,  0.5432;  0.0,  2.765;  0.0,  0.0; 0.0,  0.0];

  data.assignedPoles:={Complex(-29.4986), Complex(-10.0922), Complex(2.5201,6.89), Complex(2.5201,-6.89)};
  data.K:=fill(0, 0, 0);
end data_Byers6;

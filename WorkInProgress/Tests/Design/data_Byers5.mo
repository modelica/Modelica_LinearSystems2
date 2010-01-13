within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Byers5 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;

  output DesignData data(
  redeclare Real A[5,5],
  redeclare Real B[5,2],
  redeclare Complex assignedPoles[5]);

algorithm
  data.A:=diagonal({1, 10, 0.1, 0.1, 10})*[-1.29e-1,  0.0,  3.96e-2, 2.5e-2,   1.91e-2;
                                                3.29e-3,  0.0, -7.79e-5, 1.22e-4, -6.21e-1;
                                                7.18e-2,  0.0, -1.0e-1,  8.87e-4,  -3.85;
                                                4.11e-2,  0.0,  0.0,    -8.22e-2,   0.0;
                                                3.51e-4,  0.0,  3.5e-5,  4.26e-5,  -7.43e-2]*diagonal({1, 0.1, 10, 10, 0.1});
  data.B:=diagonal({1, 10, 0.1, 0.1, 10})*[0.0,  1.39e-3; 0.0,  3.59e-5; 0.0,  -9.89e-3; 2.49e-5,  0.0; 0.0,  -5.34e-6]*diagonal({10000, 100});
  data.assignedPoles:=Complex(1)*{-0.01,  -0.02,  -0.03,  -0.04,  -0.05};
  data.K:=fill(0, 0, 0);
end data_Byers5;

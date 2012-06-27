within Modelica_LinearSystems2.WorkInProgress;
function LinearSystems2Tests

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.StateSpace;

protected
  StateSpace sys1 "State system 1";
  StateSpace sys2 "State system 1";

algorithm
  sys1 := Modelica_LinearSystems2.StateSpace.'constructor'.fromABCDMatrices([-0.4371, 1; -0.3593, -0.0970], [-0.1220; -0.6353], [1, 0; 0,1; 31.5183, 0; 32.0534, -0.1450], [0;0;8.799;7.8492], {"dm"}, {"alpha", "q", "Gz_cdg", "Gz_sen"}, {"alpha", "q"});
  sys2 := Modelica_LinearSystems2.StateSpace.'constructor'.fromABCDMatrices([-0.4371, 1; 0.3593, -0.0970], [-0.1220; -0.6353], [1, 0; 0,1; 31.5183, 0; 32.0534, -0.1450], [0;0;8.799;7.8492], {"dm"}, {"alpha", "q", "Gz_cdg", "Gz_sen"}, {"alpha", "q"});

  print("page 10");
  Modelica_LinearSystems2.StateSpace.Transformation.extract(sys1, {4},{1});

  print("page 11");
  Modelica_LinearSystems2.StateSpace.Analysis.eigenVectors(sys1,true);
  Modelica_LinearSystems2.StateSpace.Analysis.eigenVectors(sys2,true);

  print("page 14");
  // Modelica_LinearSystems2.StateSpace.Analysis.zerosAndPoles(sys1);

  print("page 16");
  Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunctionMIMO(sys1);
  Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunctionMIMO(sys2);

end LinearSystems2Tests;

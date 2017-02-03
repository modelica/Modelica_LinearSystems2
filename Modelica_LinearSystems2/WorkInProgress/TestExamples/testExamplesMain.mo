within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesMain "Test all examples from package "

algorithm
  Modelica_LinearSystems2.Examples.StateSpace();
  Modelica_LinearSystems2.Examples.ZerosAndPoles();
  Modelica_LinearSystems2.Examples.TransferFunction();
  Modelica_LinearSystems2.Examples.DiscreteStateSpace();
  Modelica_LinearSystems2.Examples.DiscreteZerosAndPoles();
  Modelica_LinearSystems2.Examples.DiscreteTransferFunction();

  Modelica.Utilities.Streams.print("testExamplesMain done!");
end testExamplesMain;

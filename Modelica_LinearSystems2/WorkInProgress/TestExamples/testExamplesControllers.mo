within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesControllers "Test all examples from package Controllers"

algorithm
  simulateModel("Modelica_LinearSystems2.Controller.Examples.FirstExample");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.SimpleControlledDrive");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.Discretization1");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.Discretization2");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.DiscretizationSeries");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.Interpolator");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.DoublePendulum");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.InverseDoublePendulum");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.InverseDoublePendulumWithObserver");
  simulateModel("Modelica_LinearSystems2.Controller.Examples.MixingUnit");

  Modelica.Utilities.Streams.print("testExamplesControllers done!");
end testExamplesControllers;

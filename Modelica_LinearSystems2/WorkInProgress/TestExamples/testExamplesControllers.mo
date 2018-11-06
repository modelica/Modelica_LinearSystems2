within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesControllers "Test all examples from package Controllers"
  import DymolaCommands.SimulatorAPI.simulateModel;

  output Boolean ok "True, if all OK";
algorithm
  ok := false;
  ok := simulateModel("Modelica_LinearSystems2.Controllers.Examples.FirstExample");
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.SimpleControlledDrive")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.Discretization1")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.Discretization2")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.DiscretizationSeries")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.Interpolator")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.DoublePendulum")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.InverseDoublePendulum")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.InverseDoublePendulumWithObserver")
        else false;
  ok := if ok then simulateModel("Modelica_LinearSystems2.Controllers.Examples.MixingUnit")
        else false;

  if ok then
    Modelica.Utilities.Streams.print("testExamplesControllers done!");
  else
    Modelica.Utilities.Streams.print("testExamplesControllers FAILED!");
  end if;
end testExamplesControllers;

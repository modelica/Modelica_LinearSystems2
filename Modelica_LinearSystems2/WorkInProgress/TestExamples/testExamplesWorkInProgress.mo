within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesWorkInProgress "Test all examples from package WorkInProgress"
  import DymolaCommands.SimulatorAPI.simulateModel;

algorithm
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Controller.Examples.TestComponents");
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Controller.Examples.ZerosAndPolesBlock");
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Controller.Examples.limIntegrator");
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designCraneController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designCraneControllerWithObserver();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInverseDoublePendulumController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInverseDoublePendulumControllerWithObserver();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInversePendulumController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designStateSpaceController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.analysis2();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.plotBodeSISODiscrete();
  Modelica_LinearSystems2.WorkInProgress.TransferFunction.Examples.plotBodeDiscrete();
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Tests.Examples.SimpleDrive_SISO");
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Tests.Examples.Motor");
  simulateModel("Modelica_LinearSystems2.WorkInProgress.Tests.Examples.SimpleStateSpaceSystem");

  Modelica.Utilities.Streams.print("testExamplesWorkInProgress done!");
  annotation (
    __Dymola_interactive=true);
end testExamplesWorkInProgress;

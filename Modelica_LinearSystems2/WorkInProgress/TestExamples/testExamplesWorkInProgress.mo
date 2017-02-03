within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesWorkInProgress "Test all examples from package WorkInProgress"

algorithm
  Modelica_LinearSystems2.WorkInProgress.Controller.Examples.TestComponents();
  Modelica_LinearSystems2.WorkInProgress.Controller.Examples.ZerosAndPolesBlock();
  Modelica_LinearSystems2.WorkInProgress.Controller.Examples.limIntegrator();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designCraneController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designCraneControllerWithObserver();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInverseDoublePendulumController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInverseDoublePendulumControllerWithObserver();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designInversePendulumController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.designStateSpaceController();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.analysis2();
  Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples.plotBodeSISODiscrete();
  Modelica_LinearSystems2.WorkInProgress.TransferFunction.Examples.plotBodeDiscrete();
  Modelica_LinearSystems2.WorkInProgress.Tests.Examples.SimpleDrive_SISO();
  Modelica_LinearSystems2.WorkInProgress.Tests.Examples.Motor();
  Modelica_LinearSystems2.WorkInProgress.Tests.Examples.SimpleStateSpaceSystem();

  Modelica.Utilities.Streams.print("testExamplesWorkInProgress done!");
end testExamplesWorkInProgress;

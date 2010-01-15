within Modelica_LinearSystems2.WorkInProgress.Tests.Controller;
model DiscretizationInitialState
  extends Modelica_LinearSystems2.Controller.Examples.Discretization(
      sampleClock(initType=Modelica_LinearSystems2.Controller.Types.Init.InitialState));
        annotation (interactive=true);
end DiscretizationInitialState;

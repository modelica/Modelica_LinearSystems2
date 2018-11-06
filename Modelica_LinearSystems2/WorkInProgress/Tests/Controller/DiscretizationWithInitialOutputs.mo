within Modelica_LinearSystems2.WorkInProgress.Tests.Controller;
model DiscretizationWithInitialOutputs
  parameter Real y_start = step.offset;
  extends Modelica_LinearSystems2.Controllers.Examples.Discretization1(
    sampleClock(initType=Modelica_LinearSystems2.Controllers.Types.Init.InitialOutput),
    continuous(y_start=y_start),
    explicitEuler(y_start=y_start),
    implicitEuler(y_start=y_start),
    trapezoid(y_start=y_start),
    impulseExact(y_start=y_start),
    stepExact(y_start=y_start),
    rampExact(y_start=y_start));

end DiscretizationWithInitialOutputs;

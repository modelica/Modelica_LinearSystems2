within Modelica_LinearSystems2.Tests;
package Controller
  model DiscretizationInitialState
    extends Modelica_LinearSystems2.Controller.Examples.Discretization(
        sampleClock(initType=Modelica_LinearSystems2.Controller.Types.Init.InitialState));
  end DiscretizationInitialState;

  model DiscretizationWithInitialOutputs
    parameter Real y_start = step.offset;
    extends Modelica_LinearSystems2.Controller.Examples.Discretization(
      sampleClock(initType=Modelica_LinearSystems2.Controller.Types.Init.InitialOutput),
      continuous(y_start=y_start),
      explicitEuler(y_start=y_start),
      implicitEuler(y_start=y_start),
      trapezoid(y_start=y_start),
      impulseExact(y_start=y_start),
      stepExact(y_start=y_start),
      rampExact(y_start=y_start));
  end DiscretizationWithInitialOutputs;
end Controller;

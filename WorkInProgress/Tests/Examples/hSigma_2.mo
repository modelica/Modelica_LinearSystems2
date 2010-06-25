within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function hSigma_2 "Pendular output function 2"
  extends Modelica_LinearSystems2.DiscreteStateSpace.Internal.hBase;

//   input Real x[:] "State at instant k";
//   input Real u[:] "Input at instant k";
//   input Modelica.SIunits.Time Ts "Sample time";
//
//   output Real y[2] "Output at instant k";

algorithm
    y := {x[3] + 5*sin(x[1]), 5*x[1]};
end hSigma_2;

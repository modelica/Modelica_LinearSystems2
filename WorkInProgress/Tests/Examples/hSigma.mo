within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function hSigma "Pendular output function"
  extends Modelica_LinearSystems2.DiscreteStateSpace.Internal.hBase;

//   input Real x[:] "State at instant k";
//   input Real u[:] "Input at instant k";
//   input Modelica.SIunits.Time Ts "Sample time";
//
//   output Real y[2] "Output at instant k";

algorithm
    y := {x[3] + 10*sin(x[1]), x[1]};
end hSigma;

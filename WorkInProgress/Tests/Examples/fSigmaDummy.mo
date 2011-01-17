within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function fSigmaDummy "Dummy function for the discretetized state function"

  import Modelica_LinearSystems2.DiscreteStateSpace.Internal;
  extends Modelica.Icons.Function;

  extends Internal.fBase;

algorithm
  x_new := x;
end fSigmaDummy;

within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
function hSigmaDummy "Dummy function for the discretetized output function"
  import Modelica_LinearSystems2;
  extends Modelica_LinearSystems2.DiscreteStateSpace.Internal.hBase;

algorithm
    y := fill(1,ny);
end hSigmaDummy;

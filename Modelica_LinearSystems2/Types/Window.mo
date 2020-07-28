within Modelica_LinearSystems2.Types;
type Window = enumeration(
    Rectangle "Rectangle",
    Bartlett "Bartlett",
    Hann "Hann",
    Hamming "Hamming",
    Blackman "Blackman",
    Kaiser "Kaiser") "Enumeration of window type"
    annotation (
      Evaluate=true,
      obsolete = "Obsolete enumeration - use Modelica_LinearSystems2.Utilities.Types.Window instead");

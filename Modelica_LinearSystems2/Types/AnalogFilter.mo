within Modelica_LinearSystems2.Types;
type AnalogFilter = enumeration(
    CriticalDamping "Filter with critical damping",
    Bessel "Bessel filter",
    Butterworth "Butterworth filter",
    Chebyshev "Chebyshev filter")
  "Enumeration defining the method of filtering"
    annotation (
      Evaluate=true,
      obsolete = "Obsolete enumeration - use Modelica_LinearSystems2.Utilities.Types.AnalogFilter instead");

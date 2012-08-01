within Modelica_LinearSystems2.Controller.Types;
type AnalogFilter = enumeration(
    CriticalDamping "Filter with critical damping",
    Bessel "Bessel filter",
    Butterworth "Butterworth filter",
    Chebyshev "Chebyshev filter")
  "Enumeration defining the method of filtering"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));

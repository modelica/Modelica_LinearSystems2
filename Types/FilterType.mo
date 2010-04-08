within Modelica_LinearSystems2.Types;
type FilterType = enumeration(
    LowPass "Low pass filter",
    HighPass "High pass filter",
    BandPass "Band pass filter",
    BandStop "Band stop / notch filter")
  "Enumeration of analog filter types (high pass or low pass)"
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

within Modelica_LinearSystems2.Controllers.Types;
type Window = enumeration(
    Rectangle "Rectangle",
    Bartlett "Bartlett",
    Hann "Hann",
    Hamming "Hamming",
    Blackman "Blackman",
    Kaiser "Kaiser") "Enumeration of window type for FIR filter"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));

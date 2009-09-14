within Modelica_LinearSystems2.Types;
type Window = enumeration(
    Rectangle "Rectangle",
    Bartlett "Bartlett",
    Hann "Hann",
    Hamming "Hamming",
    Blackman "Blackman",
    Kaiser "Kaiser") "Enumeration of window type" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

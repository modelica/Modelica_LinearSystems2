within Modelica_LinearSystems2.Controller.Types;
type FIRspec = enumeration(
    MeanValue "Mean value FIR filter",
    Window "FIR filter defined by window",
    Coefficients "FIR filter defined by coefficient vector")
  "Enumeration of FIR specification" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

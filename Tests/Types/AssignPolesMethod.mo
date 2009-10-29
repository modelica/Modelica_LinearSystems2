within Modelica_LinearSystems2.Tests.Types;
type AssignPolesMethod = enumeration(
    Schur "Schur method for pole assignment",
    KNV "Robust pole assignment accordimg to MATLAB's place-algorithm ")
  "Enumeration of algorithms for pole aassignment" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

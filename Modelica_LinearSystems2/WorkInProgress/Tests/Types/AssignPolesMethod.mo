within Modelica_LinearSystems2.WorkInProgress.Tests.Types;
type AssignPolesMethod = enumeration(
    Schur "Schur, Schur method for pole assignment",
    KNV "KNV, Robust pole assignment accordimg to MATLAB's place-algorithm ")
  "Enumeration of algorithms for pole assignment"
    annotation (Evaluate=true, Documentation(info="<html>

</html>"));

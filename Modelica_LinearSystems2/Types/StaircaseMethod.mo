within Modelica_LinearSystems2.Types;
type StaircaseMethod = enumeration(
    QR "Apply staircase algorithm based on QR factorization",
    SVD "Apply staircase algorithm based on SVD")
  "Enumeration of methods for staircase algorithm"
    annotation (
      Evaluate=true,
      obsolete = "Obsolete enumeration - use Modelica_LinearSystems2.Utilities.Types.StaircaseMethod instead");

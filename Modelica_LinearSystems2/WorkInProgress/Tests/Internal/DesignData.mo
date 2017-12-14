within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
record DesignData
  "Contains the system matrix A, the input matrix B, the assigned Poles and optionally the ideal feedback matrix"
  import Complex;
  extends Modelica.Icons.Record;

  Real A[:,size(A,1)] "System matrix";
  Real B[size(A,1),:] "Input matrix";
  Complex assignedPoles[:] "Assigned poles";
  Real K[:,:] "Feedback gain matrix";

end DesignData;

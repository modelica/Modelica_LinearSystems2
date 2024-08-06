within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
record DesignData_Chow_Kokotovic
  "Contains the system matrix A, the input matrix B, the assigned Poles and optionally the ideal feedback matrix"
  import Complex;
  import
    Modelica_LinearSystems2.WorkInProgress.Tests.Design.DesignData_Chow_Kokotovic;
  extends Modelica.Icons.Record;
  Real A[:,size(A,1)] "System matrix";
  Real B[size(A,1),:] "Input matrix";
  Complex assignedPoles[:] "Assigned poles";
  Real K[:,:] "Feedback gain matrix";

  encapsulated operator 'constructor'
    "Default constructors for a DesignData record"
    import Modelica_LinearSystems2;

  function parameters "Default constructor for a DesignData record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Complex;
      import
        Modelica_LinearSystems2.WorkInProgress.Tests.Design.DesignData_Chow_Kokotovic;

    input Real d=1e-6;
    output DesignData_Chow_Kokotovic data(
      redeclare Real A[4,4],
      redeclare Real B[4,1],
      redeclare Complex assignedPoles[4]);

  algorithm
    data.A := [0,0.4,0,0; 0,0,0.345,0; 0,-0.524/d,-0.465/d,0.262/d; 0,0,0,-1/d];
    data.B := [0; 0; 0; 1/d];
    data.assignedPoles := Complex(1)*{-1.0,-2.0,-3.0,-4.0};
    data.K := fill(0, 0, 0);

  end parameters;

  end 'constructor';

end DesignData_Chow_Kokotovic;

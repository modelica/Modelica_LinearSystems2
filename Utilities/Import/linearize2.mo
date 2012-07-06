within Modelica_LinearSystems2.Utilities.Import;
function linearize2
  "Linearize a model after simulation up to a given time and return only the A matrix"
  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica.SIunits.Time t_linearize= 0
    "Simulate until t_linearize and then linearize" annotation(Dialog);
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options it t_linearize > 0" annotation(Dialog);

protected
  String fileName="dslin";
  String fileName2=fileName+".mat";

  // Simulate until t_linearize and then linearize at this time instant
  Boolean OK1 = if t_linearize <= 0.0 then true else
                   simulateModel(problem=modelName, startTime=0, stopTime=t_linearize,
                                 method=simulationSetup.method,
                                 tolerance=simulationSetup.tolerance,
                                 fixedstepsize=simulationSetup.fixedStepSize);
  Boolean OK2 = if t_linearize <= 0.0 then true else importInitial("dsfinal.txt");
  Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

  // Read linear system from file
  Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
  Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
  Integer nx=integer(nxMat[1, 1]);
  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
public
  output Real A[nx,nx] =  ABCD[1:nx, 1:nx] "A-matrix";
algorithm

   annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This function is identical to function
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Import.linearize\">linearize</a>
but returns only the A-matrix.
</p>
</html>"));
end linearize2;

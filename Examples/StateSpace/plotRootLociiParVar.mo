within Modelica_LinearSystems2.Examples.StateSpace;
function plotRootLociiParVar
  "Plot root locii of Modelica model with variations of one model parameter"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;

  input String modelName="Modelica_LinearSystems2.Utilities.Plants.DoublePendulum"
    "Name of the Modelica model"
    annotation(Dialog(translatedModel));
  input String modelParameter="length" "Parameter to be varied"
    annotation(Dialog(__Dymola_importDsin(onlyStart=true,fields(name=initialName))));
  input Integer nVariations(min=1) = 1 "Number of variations modelParameter";
  input Real minP = 1 "Minimum value of modelParameter";
  input Real maxP = 2 "Maximum value of modelParameter";

  input Real T_linearize=0
    "Simulate until T_linearize and then linearize the model"
    annotation (Dialog(group="Simulation setup"));
  input String method="Dassl" "Integration method"
    annotation (Dialog(group="Simulation setup"));
  output String modelName2 "Name of the model with parameter and its value";
  //input Modelica_LinearSystems2.StateSpace resultX;
protected
  String fileName="dslin" "Name of the result file";
  String fileName2=fileName+".mat" "Name of the result file with extension";
  Boolean ok "True, if all calls are ok";
  // StateSpace resultX;
  StateSpace resultX=StateSpace.Import.fromModel(modelName, 0, fileName, method);

//   Real nxMat[1,1];
//   Integer nx;
// //  Integer ABCDsizes[2]=readMatrixSize(fileName, "ABCD");
//   Real Amat[nx,nx];
//   String xuyName[nx];
//  Real eigenvalues[size(A, 1), 2]
  Real dp "Step of parameter equidistand grid";
  Real parValue "Value of parameter at a loop";
algorithm
  dp := (maxP-minP)/max(nVariations-1,1);
  Modelica.Utilities.Streams.print("------------ Step 1");
  for i in 1:nVariations loop
    parValue := minP+(i-1)*dp;
    modelName2 := modelName+"("+modelParameter+"="+String(parValue)+")";
    ok := simulateModel(problem=modelName2, startTime=0, stopTime=T_linearize, method=method);
    ok := importInitial("dsfinal.txt");
    ok := linearizeModel(problem=modelName2, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize+1, method=method);
  Modelica.Utilities.Streams.print("------------ Step 2:");

    // Read matrix Amat from result file after linearization
    resultX:=Modelica_LinearSystems2.StateSpace.Internal.read_dslin(fileName);

    // Plot root locii
    Modelica_LinearSystems2.Examples.StateSpace.plotEigenvalues(
      resultX.A,
      if i == 1 then true else false,
      "Root locii of "+modelName+"; parVar: "+modelParameter,
       String(parValue),
      i,
      {0,0,255});
    Modelica.Utilities.Streams.print("------------ Step 3");
  end for;
end plotRootLociiParVar;

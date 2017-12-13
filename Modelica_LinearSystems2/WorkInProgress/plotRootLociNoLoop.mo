within Modelica_LinearSystems2.WorkInProgress;
function plotRootLociNoLoop
  "Plot root loci of nonlinear Modelica model by linearizing the model for variations of one model parameter"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

  input String modelName="Modelica.Mechanics.Rotational.Examples.First"
    "Name of the Modelica model"
    annotation(Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.WorkInProgress.Internal.ModelParameters modelParams[:]=
    {Modelica_LinearSystems2.WorkInProgress.Internal.ModelParameters(
    parName="Jload",
    parMin=1,
    parMax=6,
    nVar=10)}
    annotation (Dialog(__Dymola_label="Model parameters", __Dymola_importDsin(
        button="Select model parameter" "Select parameters to be optimized",
        onlyStart=true,
        fields(
          parName=initialName,
          parValue=initialValue.value,
          parMin=initialValue.minimum,
          parMax=initialValue.maximum,
          parUnit=initialValue.unit))));
  input Modelica_LinearSystems2.WorkInProgress.Internal.LinearizationOptions simulationOptions=
    Modelica_LinearSystems2.WorkInProgress.Internal.LinearizationOptions(
    method="Dassl",
    tolerance=0.0001)
    annotation (Dialog(__Dymola_label="Simulation setup"));

//   input Modelica_LinearSystems2.WorkInProgress.Internal.ModelSetup modelSetup
//     "Setup of the model to be analysed";

//   input String modelParameter="Jload" "Parameter to be varied"
//     annotation(Dialog(__Dymola_importDsin(onlyStart=true,fields(name=initialName))));
//   input Integer nVariations(min=2) = 10 "Number of variations modelParameter";
//   input Real p_min = 1 "Minimum value of modelParameter";
//   input Real p_max = 6 "Maximum value of modelParameter";

//   input Modelica.SIunits.Time T_linearize=0
//     "Simulate until T_linearize and then linearize the model"
//     annotation (Dialog(group="Simulation setup"));
//   input String method="Dassl" "Integration method"
//     annotation (Dialog(group="Simulation setup"),
//     choices(choice="Dassl" "Dassl",
//     choice="Rkfix2" "Rkfix2",
//     choice="Rkfix3" "Rkfix3",
//     choice="Rkfix4" "Rkfix4"));

  input Integer position[4]={5, 5, 600, 450} "Window Position"
    annotation (Dialog(group="Plot settings"));
  input Boolean useLegend = true "Use legend"
    annotation (Dialog(group="Plot settings"));
  input Boolean grid = true "Add grid"
    annotation (Dialog(group="Plot settings"));
  input MarkerStyles markerStyle=MarkerStyles.Square "Style of marker"
    annotation (Dialog(group="Plot settings"));
  input Integer markerColorMin[3]={0,0,255}
    "Color of marker for minimum parameter value"
    annotation(Dialog(group="Plot settings", colorSelector=true));
  input Integer markerColorMax[3]={255,0,0}
    "Color of marker for maximum parameter value"
    annotation(Dialog(group="Plot settings", colorSelector=true));
protected
  String fileName="dslin" "Name of the result file";
  String fileName2 "Name of the result file with extension";
  Boolean ok "True, if all calls are ok";
  Real dp "Step of parameter equidistand grid";
  Real parValue "Value of parameter in a loop";
  Real parValues[:] "Vector of all parameter values";
  Integer nVarMin = max(modelParams[1].nVar,2)
    "Minimum number of variations (at least 2 are necessary)";
  String method = simulationOptions.method "Integration method";

  StateSpace sSpace=StateSpace.Import.fromModel(modelName, 0, fileName, method);
  Real color[nVarMin,3] "Color of markers in a loop";
algorithm
  dp := (modelParams[1].parMax-modelParams[1].parMin)/(nVarMin-1);
  color := [linspace(markerColorMin[1],markerColorMax[1],nVarMin),linspace(markerColorMin[2],markerColorMax[2],nVarMin),linspace(markerColorMin[3],markerColorMax[3],nVarMin)];

  parValues:=linspace(modelParams[1].parMin,modelParams[1].parMax,nVarMin);
  ok := translateModel(modelName);
  assert(ok, "Translation of model " + modelName + " failed.");
  ok:=simulateMultiExtendedModel(
    problem=modelName,
    startTime=0,
    stopTime=0,
    initialNames={modelParams[1].parName,"linearize:"},
    initialValues=transpose({parValues,linspace(1,nVarMin,nVarMin)}),
    finalNames=fill("",0),
    method=method,
    tolerance=simulationOptions.tolerance,
    fixedstepsize=simulationOptions.fixedStepSize);

  // input String problem := "" "Name of model, e.g. Modelica.Mechanics.Rotational.Components.Clutch";
  // input Real startTime := 0.0 "Start of simulation";
  // input Real stopTime := 1.0 "End of simulation";
  // input Integer numberOfIntervals := 0 "Number of output points";
  // input Real outputInterval := 0.0 "Distance between output points";
  // input String method := "Dassl" "Integration method";
  // input Real tolerance := 0.0001 "Tolerance of integration";
  // input Real fixedstepsize := 0.0 "Fixed step size for Euler";
  // input String resultFile := "dsres" "Where to store result";
  // output Boolean result "true if successful";
  // input String initialNames[:] "Parameters and start-values to set";
  // input Real initialValues[:, size(initialNames, 1)] "Parameter values";
  // input String finalNames[:] "Variables at end-point";
  // output Real finalValues[size(initialValues, 1), size(finalNames, 1)] "Values at end-point";

  for i in 1:nVarMin loop
    // Read matrix Amat from result file after linearization
    fileName2 := fileName+String(i);
    sSpace:=Modelica_LinearSystems2.StateSpace.Internal.read_dslin(fileName2);

    // Plot root loci
    Modelica_LinearSystems2.WorkInProgress.plotEigenvalues(
      sSpace.A,
      if i == 1 then true else false,
      position,
      "Root loci of " + modelName + "; parVar: " + modelParams[1].parName,
      useLegend,
      String(parValues[i]),
      grid,
      markerStyle,
      integer(color[i, :]));
  end for;
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>plotRootLoci(modelName, modelParams, simulationOptions, position, useLegend, grid, markerStyle, markerColorMin, markerColorMax)</pre></blockquote>

<h4>Description</h4>
<p>Calling this function with default setup the following root loci plot will be generated.</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
end plotRootLociNoLoop;

within Modelica_LinearSystems2.WorkInProgress;
function plotRootLoci
  "Plot root loci of nonlinear Modelica model by linearizing the model for variations of one model parameter"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

  input String modelName="Modelica.Mechanics.Rotational.Examples.First"
    "Name of the Modelica model"
    annotation(Dialog(__Dymola_translatedModel));
  input Boolean simulate = false
    "Linearize model after simulation (time-consuming!), otherwise linearization of all parameter variants at once"
    annotation (Dialog(__Dymola_compact=false),
      choices(checkBox=true));

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

  input Integer position[4]={5, 5, 600, 450} "Window Position"
    annotation (Dialog(group="Plot settings"));
  input Boolean useLegend = true "Use legend"
    annotation (Dialog(group="Plot settings", __Dymola_compact=true, __Dymola_descriptionLabel = true),
      choices(checkBox=true));
  input Boolean grid = true "Use grid"
    annotation (Dialog(group="Plot settings", __Dymola_compact=true, __Dymola_descriptionLabel = true),
      choices(checkBox=true));
  input MarkerStyles markerStyle=MarkerStyles.Square "Style of marker"
    annotation (Dialog(group="Plot settings"));
  input Integer markerColorMin[3]={0,0,255}
    "Color of marker for minimum parameter value"
    annotation(Dialog(group="Plot settings", colorSelector=true));
  input Integer markerColorMax[3]={255,0,0}
    "Color of marker for maximum parameter value"
    annotation(Dialog(group="Plot settings", colorSelector=true));

  input Boolean deleteResult = false "Delete result files of linearization"
    annotation (Dialog(__Dymola_compact=false),
      choices(checkBox=true));
protected
  String fileName="dslin" "Name of the result file";
  String fileName2=fileName+".mat" "Name of the result file with extension";
  String modelName2 "Name of the model with parameter and its value";
  Boolean ok "True, if all calls are ok";
//  Real dp "Step of parameter equidistand grid";
//  Real parValue "Value of parameter in a loop";
  Real parValues[:] "Vector of all parameter values";
  Integer nVarMin = max(modelParams[1].nVar,2)
    "Minimum number of variations (at least 2 are necessary)";
  String method = simulationOptions.method "Integration method";

  StateSpace sSpace=StateSpace.Import.fromModel(modelName, 0, fileName, method);
  Real color[nVarMin,3] "Color of markers in a loop";
algorithm
  parValues:=linspace(modelParams[1].parMin,modelParams[1].parMax,nVarMin);
  color := [linspace(markerColorMin[1],markerColorMax[1],nVarMin),linspace(markerColorMin[2],markerColorMax[2],nVarMin),linspace(markerColorMin[3],markerColorMax[3],nVarMin)];

  if not simulate then // and simulationOptions.t_linearize==0
    // Linearization of all parameter variants at once
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
  else
    // Stepwise simulation and linearization (time-consuming)
    for i in 1:nVarMin loop
      modelName2 := modelName+"("+modelParams[1].parName+"="+String(parValues[i])+")";
      fileName2 := fileName+String(i);
      Modelica.Utilities.Streams.print("  ...linearizing "+modelName2);

      ok := simulateModel(
        problem=modelName2,
        startTime=0,
        stopTime=simulationOptions.t_linearize,
        method=method,
        numberOfIntervals=simulationOptions.numberOfIntervals,
        outputInterval=simulationOptions.outputInterval,
        tolerance = simulationOptions.tolerance,
        fixedstepsize = simulationOptions.fixedStepSize);
      ok := importInitial("dsfinal.txt");
      ok := linearizeModel(
        problem=modelName2,
        resultFile=fileName2,
        startTime=simulationOptions.t_linearize,
        stopTime=simulationOptions.t_linearize+1,
        method=method,
        tolerance = simulationOptions.tolerance,
        fixedstepsize = simulationOptions.fixedStepSize);
    end for;
  end if;

  for i in 1:nVarMin loop
    // Read matrix Amat from result file after linearization
    fileName2 := fileName+String(i);
    sSpace:=Modelica_LinearSystems2.StateSpace.Internal.read_dslin(fileName2);

    // Plot root locus
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

    if deleteResult then
      Modelica.Utilities.Files.removeFile(Modelica.Utilities.Files.fullPathName(fileName2+".mat"));
      Modelica.Utilities.Streams.print("  ... "+fileName2+".mat deleted");
    end if;
  end for;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
plotRootLoci(modelName, modelParams, simulationOptions, position, useLegend, grid, markerStyle, markerColorMin, markerColorMax)
</pre></blockquote>

<h4>Description</h4>
<p>
This function linearizes Modelica model, calculates its eigenvalues and
plots root loci for variants of selected model parameter <code>modelParams[1].parName</code>.
Optionally, the model will be simulated before linearization.
This is reasonable e.g. if the model is initially not at equilibrium and it
contains nonlinearities.
The end time of the simulation is <code>simulationOptions.t_linearize</code>.
The linearization will be performed for the states of model at that time.
Therefore, this process can generally be very time-consuming.
</p>
<p>
On the contrary, if the linearization time <code>simulationOptions.t_linearize&nbsp;=&nbsp;0</code>
the linearization is performed at once for various parameter values.
This is efficient and sufficient for most cases.
</p>

<p>Calling this function with default setup the following root loci plot will be generated.</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
end plotRootLoci;

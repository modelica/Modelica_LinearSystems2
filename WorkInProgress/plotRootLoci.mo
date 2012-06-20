within Modelica_LinearSystems2.WorkInProgress;
function plotRootLoci
  "Plot root loci of nonlinear Modelica model by linearizing the model for variations of one model parameter"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Utilities.Types.MarkerStyles;

  input String modelName="Modelica.Mechanics.Rotational.Examples.First"
    "Name of the Modelica model"
    annotation(Dialog(translatedModel));
  input String modelParameter="Jload" "Parameter to be varied"
    annotation(Dialog(__Dymola_importDsin(onlyStart=true,fields(name=initialName))));
  input Integer nVariations(min=2) = 10 "Number of variations modelParameter";
  input Real p_min = 1 "Minimum value of modelParameter";
  input Real p_max = 6 "Maximum value of modelParameter";

  input Modelica.SIunits.Time T_linearize=0
    "Simulate until T_linearize and then linearize the model"
    annotation (Dialog(group="Simulation setup"));
  input String method="Dassl" "Integration method"
    annotation (Dialog(group="Simulation setup"),
    choices(choice="Dassl" "Dassl",
    choice="Rkfix2" "Rkfix2",
    choice="Rkfix3" "Rkfix3",
    choice="Rkfix4" "Rkfix4"));

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
  String fileName2=fileName+".mat" "Name of the result file with extension";
  String modelName2 "Name of the model with parameter and its value";
  Boolean ok "True, if all calls are ok";
  Real dp "Step of parameter equidistand grid";
  Real parValue "Value of parameter in a loop";
  Integer nVarMin = max(nVariations,2)
    "Minimum number of variations (at least 2 are necessary)";

  StateSpace sSpace=StateSpace.Import.fromModel(modelName, 0, fileName, method);
public
  output Real color[nVarMin,3] "Color of markers in a loop";
algorithm
  dp := (p_max-p_min)/(nVarMin-1);
  color := [linspace(markerColorMin[1],markerColorMax[1],nVarMin),linspace(markerColorMin[2],markerColorMax[2],nVarMin),linspace(markerColorMin[3],markerColorMax[3],nVarMin)];
  for i in 1:nVarMin loop
    parValue := p_min+(i-1)*dp;
    modelName2 := modelName+"("+modelParameter+"="+String(parValue)+")";
    ok := simulateModel(problem=modelName2, startTime=0, stopTime=T_linearize, method=method);
    ok := importInitial("dsfinal.txt");
    ok := linearizeModel(problem=modelName2, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize+1, method=method);

    // Read matrix Amat from result file after linearization
    sSpace:=Modelica_LinearSystems2.StateSpace.Internal.read_dslin(fileName);

    // Plot root loci
    Modelica_LinearSystems2.Examples.StateSpace.plotEigenvalues(sSpace.A, if i == 1 then true else false, position, "Root loci of " + modelName + "; parVar: " + modelParameter, useLegend, String(parValue), grid, markerStyle, integer(color[i, :]));
  end for;
  annotation (interactive=true, Documentation(info="<html>
<p><h4>Syntax</h4></p>
<blockquote><pre>plotRootLoci(modelName, modelParameter, nVariations, p_min, p_max, T_linearize, method, position, useLegend, grid, markerStyle, markerColorMin, markerColorMax)</pre></blockquote>
<p><h4>Description</h4></p>
<p>Calling this function with default setup the following root loci plot will be generated.</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
end plotRootLoci;

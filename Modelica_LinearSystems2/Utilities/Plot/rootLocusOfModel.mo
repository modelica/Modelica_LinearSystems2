within Modelica_LinearSystems2.Utilities.Plot;
function rootLocusOfModel
  "Compute and plot the root locus of one parameter of a model (= eigen values of the model that is linearized for every parameter value)"
  extends Modelica.Icons.Function;

  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.ParameterVariation modelParam[:]
    "Model parameter to be varied (exactly one) and values for other parameters";
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options" annotation(Dialog(enable=not linearizeAtInitial));
  input Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram diagram=
    Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram() annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
    Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);

protected
  Real Re[:,:]
    "Real values of eigenvalues Re[i,j], where i are the different parameter values and j the eigenvalue numbers";
  Real Im[:,:]
    "Imaginary values of eigenvalues Im[i,j], where i are the different parameter values and j the eigenvalue numbers";
  Real s[:]
    "The different parameter values s[i] associated with Re[i,j] and Im[i,j]";
  String paramName;
  String paramUnit;
  String heading;
  String pName;
  Boolean reorder=diagram.linePattern <> Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None;
algorithm
  // Compute the root locus data
  (Re, Im, s, paramName, paramUnit) :=
    Modelica_LinearSystems2.Utilities.Import.rootLocusOfModel(
    modelName,
    modelParam,
    simulationSetup,
    reorder);
  if paramUnit == "" or paramUnit == " " then
     pName :=paramName;
  else
     pName :=paramName + " [" + paramUnit + "]";
  end if;

  // Plot the data
  if diagram.heading == "" then
     heading :="Root locus of " + modelName + " over " + pName;
  else
     heading :=diagram.heading;
  end if;

  Modelica_LinearSystems2.Utilities.Plot.parameterizedCurves(
     diagram=Modelica_LinearSystems2.Utilities.Plot.Records.ParametrizedCurves(
       X=Re, Y=Im, s=s,
       xName=diagram.ReName,
       yName=diagram.ImName,
       sName=pName,
       heading=heading,
       xLabel=diagram.xLabel,
       yLabel=diagram.yLabel,
       labelWithS=diagram.labelWithParam,
       heightRatio=diagram.heightRatio,
       grid=diagram.grid,
       logX=diagram.logX,
       logY=diagram.logY,
       uniformScaling=diagram.uniformScaling,
       curveProperties={Modelica_LinearSystems2.Utilities.Plot.Records.CurveProperties(
                        lineColor=diagram.lineColor,
                        linePattern=diagram.linePattern,
                        lineSymbol=diagram.lineSymbol,
                        lineThickness=diagram.lineThickness)}),
      device=device);

  annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Plot.<b>rootLocusOfModel</b>(modelName, modelParam, simulationSetup, diagram, device)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes and plots a root locus of a selected Modelica model by
varying one parameter and performing an eigen value analysis for every parameter value.
Other parameters of the model can be set to a specific value. An equidistant or
a logarithmic gridding of the parameter to be varied can be selected and then the gridding
is performed between the given minimum and maximum value.
</p>

<h4>Example</h4>
<p>
Execute function <a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfDrive\">Examples.rootLocusOfDrive</a>
that is defined to analyse an example mechanical model from the Modelica standard library:
</p>
<p>
Calling the function
</p>
<blockquote><pre>
Utilities.Plot.<b>rootLocusOfModel</b>(
  modelName = &quot;Modelica.Mechanics.Rotational.Examples.First&quot;,
  modelParam={
    Modelica_LinearSystems2.Records.ParameterVariation(
      Name=&quot;Jload&quot;,
      grid=Modelica_LinearSystems2.Types.Grid.Equidistant
      Min=1,
      Max=6,
      nPoints=30)});
</pre></blockquote>
<p>
yields the following diagram (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same parameter value are
marked with a red square):
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/RootLocusOfModel.png\"/></p>
</html>"));
end rootLocusOfModel;

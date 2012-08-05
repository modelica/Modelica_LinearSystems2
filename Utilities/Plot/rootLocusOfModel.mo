within Modelica_LinearSystems2.Utilities.Plot;
function rootLocusOfModel
  "Compute and plot the root locus of one parameter of a model (= eigen values of the model that is linearized for every parameter value)"
  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.ParameterVariation modelParam[:]
    "Model parameter to be varied and modified values for other parameters";
  input Boolean linearizeAtInitial=true
    "= true, if linearization at inital time; otherwise simulate until t_linearize"
     annotation (choices(__Dymola_checkBox=true));
  input Modelica.SIunits.Time t_linearize= 0
    "Simulate until t_linearize and then linearize, if linearizeAtInitial == false"
                                                                                    annotation(Dialog(enable=not linearizeAtInitial));
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options it t_linearize > 0, if linearizeAtInitial == false" annotation(Dialog(enable=not linearizeAtInitial));
  input Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram diagram annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device
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
    linearizeAtInitial,
    t_linearize,
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

  annotation (__Dymola_interactive=true);
end rootLocusOfModel;

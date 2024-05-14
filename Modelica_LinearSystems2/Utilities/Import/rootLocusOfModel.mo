within Modelica_LinearSystems2.Utilities.Import;
function rootLocusOfModel
  "Return the root locus of one parameter (= eigen values of the model that is linearized for every parameter value)"
  extends Modelica.Icons.Function;

  import Simulator = DymolaCommands.SimulatorAPI;
  import Modelica_LinearSystems2.Utilities.Types.Grid;

  input String modelName "Name of the Modelica model" annotation (Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.ParameterVariation modelParam[:]
    "Model parameter to be varied (exactly one) and values for other parameters";
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization
    simulationSetup = Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options";
  input Boolean reorder=false
    "True, if eigen values shall be reordered so that they are closest to the previous ones";
  output Real Re[:, :]
    "Re[nx,np] Real values of eigenvalues Re[j,i], where i are the different parameter values and j the eigenvalue numbers";
  output Real Im[:, :]
    "Im[nx,np] Imaginary values of eigenvalues Im[j,i], where i are the different parameter values and j the eigenvalue numbers";
  output Real s[:]
    "s[np] The different parameter values s[i] associated with Re[i,j] and Im[i,j]";
  output String paramName "Name of the parameter that was varied";
  output String paramUnit "Unit of parameter paramName";
protected
  Integer nParam=size(modelParam, 1);
  Boolean OK;
  String parameterModel="_rootLocusOfOneParameter_model";
  String parameterModelFile=parameterModel + ".mo";
  String str;
  Integer index_p_var;
  Integer is[:] "File indices X of the dslinX.mat files";
  Integer np;

  String fileName="dslin";
  String fileName2;
  Integer xuy[3];
  Integer nx;
  Integer nu;
  Integer ny;
  Boolean newModel=false;
  Boolean first;
  Real Min;
  Real Max;
  Real logMin;
  Real logMax;
algorithm
  // Check that the system has eigen values
  // assert(nx > 0,"Model " + modelName + " does not has states. Therefore, root locus does not make sense.");

  // Determine the parameter to be varied and assign new parameter values to the model if necessary
  if nParam == 0 then
    // No parameter defined
    Modelica.Utilities.Streams.error(
      "No parameter defined that shall be varied for the root locus");

  elseif nParam == 1 then
    // Exactly one parameter defined
    assert(modelParam[1].grid == Grid.Equidistant or modelParam[1].grid == Grid.Logarithmic,
      "One parameter defined, but grid is not defined as Equidistant or Logarithmic");
    np := modelParam[1].nPoints;
    index_p_var := 1;
    OK := Simulator.translateModel(modelName);

  else
    // More then one parameter defined; find the parameter to be varied
    index_p_var := 0;
    for i in 1:nParam loop
      if modelParam[i].grid == Grid.OneValue then
        // do nothing
      elseif modelParam[i].grid == Grid.Equidistant or modelParam[i].grid ==
          Grid.Logarithmic then
        if index_p_var > 0 then
          Modelica.Utilities.Streams.error("Parameters " + modelParam[
            index_p_var].Name + " and " + modelParam[i].Name +
            " shall be varied,\n" +
            "but this is only possible for one parameter.\n" +
            " Therefore, change the definition of \"grid\".");
        end if;
        index_p_var := i;
      else
        assert(false, "Wrong definition of \"grid\" for parameter " +
          modelParam[i].Name);
      end if;
    end for;
    assert(index_p_var > 0,
      "No parameter defined that shall be varied for the root locus.");
    np := modelParam[index_p_var].nPoints;

    // Translate model and set the new parameter values
    OK := Simulator.translateModel(modelName);
    assert(OK, "Translation of model " + modelName + " failed.");
    for i in 1:nParam loop
      if i <> index_p_var then
        OK := SetVariable(modelParam[i].Name, modelParam[i].Value);
        assert(OK, "Setting parameter " + modelParam[i].Name + " = " + String(
          modelParam[i].Value) + " failed.");
      end if;
    end for;
  end if;

  // Parameter that is varied
  paramName := modelParam[index_p_var].Name;
  paramUnit := "";
  Modelica.Utilities.Streams.print("Perform parameter variation of " + paramName);

  // Check min/max values
  Min := modelParam[index_p_var].Min;
  Max := modelParam[index_p_var].Max;
  assert(Min > -1e99, "Minimum value not set for parameter to be varied: " +
    paramName);
  assert(Max < 1e99, "Maximum value not set for parameter to be varied: " +
    paramName);

  // Compute all parameter values
  if modelParam[index_p_var].grid == Grid.Logarithmic then
    // logarithmic spacing
    assert(Min*Max >= 0.0,
      "Since grid = Logarithmic for parameter to be varied: " + paramName +
      "\nThe Min and Max values need to have the same sign.");
    if Min < 0.0 then
      logMin := -log10(Min);
    elseif Min > 0.0 then
      logMin := log10(Min);
    else
      // Min = 0.0
      logMin := log10(1e-15);
    end if;

    if Max < 0.0 then
      logMax := -log10(Max);
    elseif Max > 0.0 then
      logMax := log10(Max);
    else
      // Max = 0.0
      logMax := -log10(1e-15);
    end if;

    s := linspace(
      logMin,
      logMax,
      np);
    for i in 1:size(s, 1) loop
      s[i] := 10^s[i];
    end for;
  else
    s := linspace(
      Min,
      Max,
      np);
  end if;

  is := 1:np;
  if simulationSetup.linearizeAtInitial then
    // Linearization of all parameter variants at once at the initial point
    OK := Simulator.simulateMultiExtendedModel(
      problem=modelName,
      startTime=0,
      stopTime=0,
      initialNames={paramName, "linearize:"},
      initialValues=[s, is],
      finalNames=fill("", 0),
      method=simulationSetup.method,
      tolerance=simulationSetup.tolerance,
      fixedstepsize=simulationSetup.fixedStepSize);
      // Disable the following assert since OK is always 'false' for initialNames[2]="linearize:".
      // This setting of initialNames[2] is, in turn, necessary to get file dslin1.mat
      // which is needed to read linear system sizes: nx, nu and ny below.
      // assert(OK,
      //       "Linearization with function simulateMultiExtendedModel failed (maybe some parameter values are not meaningful?).");
  else
    // Simulate always until t_linearize and only then linearize
    Modelica.Utilities.Streams.error("Option not yet implemented");
  end if;

  // Determine array dimensions of the first linearization point
  fileName2 := fileName + String(is[1]) + ".mat";
  xuy :=Streams.readSystemDimension(fileName2, "ABCD");
  nx :=xuy[1];
  nu :=xuy[2];
  ny :=xuy[3];

  // Read all matrices from file, compute eigenvalues and store them in output arrays
  (Re,Im) := Modelica_LinearSystems2.Internal.eigenValuesFromLinearization(
    is,
    nx,
    nu,
    ny,
    reorder);

  annotation (__Dymola_interactive=true);
end rootLocusOfModel;

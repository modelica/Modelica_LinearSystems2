within Modelica_LinearSystems2.Utilities.Import;
function linearize2
  "Linearize a model at the start time, or optionally after simulation up to a given time instant, and return it as StateSpace object"
  extends Modelica.Icons.Function;

  import Simulator = DymolaCommands.SimulatorAPI;

  input String modelName "Name of the Modelica model"
    annotation (Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.SetParameter modelParam[:]=fill(
      Modelica_LinearSystems2.Records.SetParameter(Name="", Value=0.0), 0)
    "Values of model parameters used for linearization";
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization
    simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options";
protected
  String fileName="dslin";

  function setModelParam
    extends Modelica.Icons.Function;

    input String modelName;
    input Modelica_LinearSystems2.Records.SetParameter modelParam[:];
    output Boolean OK;
  algorithm
    OK := DymolaCommands.SimulatorAPI.translateModel(modelName);
    assert(OK, "Translation of model " + modelName + " failed.");
    for i in 1:size(modelParam, 1) loop
      OK := SetVariable(modelParam[i].Name, modelParam[i].Value);
      assert(OK, "Setting parameter " + modelParam[i].Name + " = " + String(
        modelParam[i].Value) + " failed.");
    end for;

    annotation (__Dymola_interactive=true);
  end setModelParam;

  // Set model parameters
  Boolean OK1 = if size(modelParam, 1) == 0 then true else
    setModelParam(modelName, modelParam);

  // Simulate until t_linearize and then linearize at this time instant
  Boolean OK2 = if simulationSetup.linearizeAtInitial then true else
    Simulator.simulateModel(
      problem=modelName,
      startTime=simulationSetup.t_start,
      stopTime=simulationSetup.t_linearize,
      method=simulationSetup.method,
      tolerance=simulationSetup.tolerance,
      fixedstepsize=simulationSetup.fixedStepSize);
  Boolean OK3 = if simulationSetup.linearizeAtInitial then true else
    Simulator.importInitial("dsfinal.txt");
  Boolean OK4 = Simulator.linearizeModel(
      problem=modelName,
      resultFile=fileName,
      startTime=simulationSetup.t_linearize,
      stopTime=simulationSetup.t_linearize);

  // Model is already translated. Reset to the default initial conditions
  Boolean OK5 = Simulator.translateModel(problem=modelName);
public
  output Modelica_LinearSystems2.StateSpace ss=
    Modelica_LinearSystems2.StateSpace.Internal.read_dslin(fileName)
    "Linearized system as StateSpace object";
algorithm

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This function initializes a Modelica model and then simulates the model
until time instant \"t_linearize\".
If t_linearize=0, no simulation takes place (only initialization).
At the simulation stop time, the model is linearized in such a form that
</p>
<ul>
  <li>
    all top-level signals with prefix \"input\" are treated as inputs <strong>u</strong>(t)
    of the model,
  </li>
  <li>
    all top-level signals with prefix \"output\" are treated as outputs <strong>y</strong>(t)
    of the model,
  </li>
  <li>
    all variables that appear differentiated and that are selected as states at this time
    instant are treated as states <strong>x</strong> of the model.
  </li>
</ul>

<p>
Formally, the non-linear hybrid differential-algebraic equation system is therefore treated
as the following ordinary equation system at time instant t_linearize:
</p>
<blockquote><pre>
der(<strong>x</strong>) = <strong>f</strong>(<strong>x</strong>,<strong>u</strong>)
    <strong>y</strong>  = <strong>g</strong>(<strong>x</strong>,<strong>u</strong>) 
</pre></blockquote>
<p>
Taylor series expansion (linearization) of this model around the simulation stop time t_linearize:
</p>
<blockquote><pre>
<strong>u</strong>0 = <strong>u</strong>(t_linearize)
<strong>y</strong>0 = <strong>y</strong>(t_linearize)
<strong>x</strong>0 = <strong>x</strong>(t_linearize)
</pre></blockquote>
<p>
and neglecting higher order terms results in the following system:
</p>
<blockquote><pre>
der(<strong>x</strong>0 + d<strong>x</strong>) = <strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0) + der(<strong>f</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>f</strong>,<strong>u</strong>)*d<strong>u</strong>
    <strong>y</strong>0 + d<strong>y</strong>  = <strong>g</strong>(<strong>x</strong>0,<strong>u</strong>0) + der(<strong>g</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>g</strong>,<strong>u</strong>)*d<strong>u</strong>
</pre></blockquote>
<p>
where der(<strong>f</strong>,<strong>x</strong>) is the partial derivative
of&nbsp;<strong>f</strong> with respect to&nbsp;<strong>x</strong>, and the partial
derivatives are computed at the linearization point t_linearize. Re-ordering of terms
gives (note <strong>der</strong>(<strong>x</strong>0)&nbsp;=&nbsp;<strong>0</strong>):
</p>
<blockquote><pre>
der(d<strong>x</strong>) = der(<strong>f</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>f</strong>,<strong>u</strong>)*d<strong>u</strong> + <strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0)
    d<strong>y</strong>  = der(<strong>g</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>g</strong>,<strong>u</strong>)*d<strong>u</strong> + (<strong>g</strong>(<strong>x</strong>0,<strong>u</strong>0) - <strong>y</strong>0)
</pre></blockquote>
<p>
or
</p>
<blockquote><pre>
der(d<strong>x</strong>) = <strong>A</strong>*d<strong>x</strong> + <strong>B</strong>*d<strong>u</strong> + <strong>f</strong>0
    d<strong>y</strong>  = <strong>C</strong>*d<strong>x</strong> + <strong>D</strong>*d<strong>u</strong>
</pre></blockquote>
<p>
This function returns the matrices&nbsp;<strong>A</strong>, <strong>B</strong>,
<strong>C</strong>, <strong>D</strong> and assumes that the linearization point is
a&nbsp;steady-state point of the simulation (i.e.,
<strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0)&nbsp;=&nbsp;0).
Additionally, the full Modelica names of all inputs, outputs and states shall be
returned if possible (default is to return empty name strings).
</p>
</html>"));
end linearize2;

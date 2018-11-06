within Modelica_LinearSystems2.Utilities.Import;
function linearize2
  "Linearize a model at the start time, or optionally after simulation up to a given time instant, and return it as StateSpace object"
  extends Modelica.Icons.Function;
  import Modelica.Utilities.Streams;
  import Modelica.Utilities.Streams.print;
  import Simulator = DymolaCommands.SimulatorAPI;

  input String modelName "Name of the Modelica model"
    annotation (Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.SetParameter modelParam[:] = fill(
    Modelica_LinearSystems2.Records.SetParameter(Name="", Value=0.0), 0)
    "Values of model parameters used for linearization";
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization
    simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options" annotation (Dialog(enable=not linearizeAtInitial));
protected
  String fileName = "dslin";

  function setModelParam
    input String modelName;
    input Modelica_LinearSystems2.Records.SetParameter modelParam[:];
    output Boolean OK;
  algorithm
    OK := Simulator.translateModel(modelName);
    assert(OK, "Translation of model " + modelName + " failed.");
    for i in 1:size(modelParam, 1) loop
      OK := SetVariable(modelParam[i].Name, modelParam[i].Value);
      assert(OK, "Setting parameter " + modelParam[i].Name + " = " + String(
        modelParam[i].Value) + " failed.");
    end for;
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

  annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = Utilities.Import.linearize2(
  modelName,modelParam,simulationSetup)
</pre></blockquote>

<h4>Description</h4>
<p>
This function initializes a Modelica model and then simulates the model
until time instant &quot;t_linearize&quot;.
If t_linearize=0, no simulation takes place (only initialization).
At the simulation stop time, the model is linearized in such a form that
</p>
<ul>
<li>all top-level signals with prefix &quot;input&quot; are treated as inputs <b>u</b>(t) of the model ,</li>
<li>all top-level signals with prefix &quot;output&quot; are treated as outputs <b>y</b>(t) of the model,</li>
<li>all variables that appear differentiated and that are selected as states at this time instant are treated as states <b>x</b> of the model.</li>
</ul>
<p>
Formally, the non-linear hybrid differential-algebraic equation system is
therefore treated as the following ordinary equation system at time
instant t_linearize:
</p>
<blockquote><pre>
der(<b>x</b>) = <b>f</b>(<b>x</b>,<b>u</b>)
     <b>y</b> = <b>g</b>(<b>x</b>,<b>u</b>)
</pre></blockquote>
<p>
Taylor series expansion (linearization) of this model around the simulation stop time t_linearize:
</p>
<blockquote><pre>
<b>u</b>0 = <b>u</b>(t_linearize)
<b>y</b>0 = <b>y</b>(t_linearize)
<b>x</b>0 = <b>x</b>(t_linearize)
</pre></blockquote>
<p>and neglecting higher order terms results in the following system: </p>
<blockquote>
<pre>
der(<b>x</b>0+d<b>x</b>) = <b>f</b>(<b>x</b>0,<b>u</b>0) + der(<b>f</b>,<b>x</b>)*d<b>x</b> + der(<b>f</b>,<b>u</b>)*d<b>u</b>
   <b>y</b>0 + d<b>y</b> = <b>g</b>(<b>x</b>0,<b>u</b>0) + der(<b>g</b>,<b>x</b>)*d<b>x</b> + der(<b>g</b>,<b>u</b>)*d<b>u</b>
</pre>
</blockquote>
<p>where der(<b>f</b>,<b>x</b>) is the partial derivative of <b>f</b> with respect to <b>x</b>, and the partial derivatives are computed at the linearization point t_linearize. Re-ordering of terms gives (note <b>der</b>(<b>x</b>0) = <b>0</b>): </p>
<blockquote><pre>
der(d<b>x</b>) = der(<b>f</b>,<b>x</b>)*d<b>x</b> + der(<b>f</b>,<b>u</b>)*d<b>u</b> + <b>f</b>(<b>x</b>0,<b>u</b>0)
     d<b>y</b> = der(<b>g</b>,<b>x</b>)*d<b>x</b> + der(<b>g</b>,<b>u</b>)*d<b>u</b> + (<b>g</b>(<b>x</b>0,<b>u</b>0) - <b>y</b>0)
</pre></blockquote>
<p>
or
</p>
<blockquote><pre>
der(d<b>x</b>) = <b>A</b>*d<b>x</b> + <b>B</b>*d<b>u</b> + <b>f</b>0
     d<b>y</b> = <b>C</b>*d<b>x</b> + <b>D</b>*d<b>u</b>
</pre></blockquote>
<p>
This function returns the state-space representation of the linearized system assuming that
the linearization point is a steady-state point of the simulation
(i.e., <b>f</b>(<b>x</b>0,&nbsp;<b>u</b>0)&nbsp;=&nbsp;0).
</p>
</html>"));
end linearize2;

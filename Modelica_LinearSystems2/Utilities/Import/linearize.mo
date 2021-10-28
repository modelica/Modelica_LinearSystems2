within Modelica_LinearSystems2.Utilities.Import;
function linearize "Linearize a model after simulation up to a given time"
  import Simulator = DymolaCommands.SimulatorAPI;

  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica.Units.SI.Time t_linearize=0
    "Simulate until T_linearize and then linearize" annotation (Dialog);

protected
  String fileName="dslin";
  String fileName2=fileName+".mat";

  // Simulate until t_linearize and then linearize at this time instant
  Boolean OK1 = Simulator.simulateModel(problem=modelName, startTime=0, stopTime=t_linearize);
  Boolean OK2 = Simulator.importInitial("dsfinal.txt");
  Boolean OK3 = Simulator.linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

  // Read linear system from file
  Integer xuy[3]=Streams.readSystemDimension(fileName2, "ABCD");
  Integer nx = xuy[1];
  Integer nu = xuy[2];
  Integer ny = xuy[3];
  Real ABCD[nx + ny,nx + nu] = Modelica.Utilities.Streams.readRealMatrix(
    fileName2, "ABCD", nx + ny, nx + nu);
  String xuyName[nx + nu + ny] = DymolaCommands.MatrixIO.readStringMatrix(
    fileName2, "xuyName", nx + nu + ny);

  // Model is already translated. Reset to the default initial conditions
  Boolean OK4 = Simulator.translateModel(problem=modelName);
public
  output Real A[nx,nx] =  ABCD[1:nx, 1:nx] "A-matrix";
  output Real B[nx,nu] =  ABCD[1:nx, nx + 1:nx + nu] "B-matrix";
  output Real C[ny,nx] =  ABCD[nx + 1:nx + ny, 1:nx] "C-matrix";
  output Real D[ny,nu] =  ABCD[nx + 1:nx + ny, nx + 1:nx + nu] "D-matrix";
  output String inputNames[nu] =  xuyName[nx + 1:nx + nu]
    "Modelica names of inputs";
  output String outputNames[ny] =  xuyName[nx + nu + 1:nx + nu + ny]
    "Modelica names of outputs";
  output String stateNames[nx] =  xuyName[1:nx] "Modelica names of states";
algorithm

   annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>This function initializes a Modelica model and then simulates the model with its default experiment options until time instant \"t_linearize\". If t_linearize=0, no simulation takes place (only initialization). At the simulation stop time, the model is linearized in such a form that </p>
<p><ul>
<li>all top-level signals with prefix \"input\" are treated as inputs <strong>u</strong>(t) of the model ,</li>
<li>all top-level signals with prefix \"output\" are treated as outputs <strong>y</strong>(t) of the model,</li>
<li>all variables that appear differentiated and that are selected as states at this time instant are treated as states <strong>x</strong> of the model.</li>
</ul></p>
<p>Formally, the non-linear hybrid differential-algebraic equation system is therefore treated as the following ordinary equation system at time instant t_linearize: </p>
<pre>    der(<strong>x</strong>) = <strong>f</strong>(<strong>x</strong>,<strong>u</strong>)</pre>
<pre>         <strong>y</strong> = <strong>g</strong>(<strong>x</strong>,<strong>u</strong>) </pre>
<p>Taylor series expansion (linearization) of this model around the simulation stop time t_linearize: </p>
<pre>   <strong>u</strong>0 = <strong>u</strong>(t_linearize)</pre>
<pre>   <strong>y</strong>0 = <strong>y</strong>(t_linearize)</pre>
<pre>   <strong>x</strong>0 = <strong>x</strong>(t_linearize) </pre>
<p>and neglecting higher order terms results in the following system: </p>
<pre>   der(<strong>x</strong>0+d<strong>x</strong>) = <strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0) + der(<strong>f</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>f</strong>,<strong>u</strong>)*d<strong>u</strong></pre>
<pre>      <strong>y</strong>0 + d<strong>y</strong> = <strong>g</strong>(<strong>x</strong>0,<strong>u</strong>0) + der(<strong>g</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>g</strong>,<strong>u</strong>)*d<strong>u</strong></pre>
<p>where der(<strong>f</strong>,<strong>x</strong>) is the partial derivative of <strong>f</strong> with respect to <strong>x</strong>, and the partial derivatives are computed at the linearization point t_linearize. Re-ordering of terms gives (note <strong>der</strong>(<strong>x</strong>0) = <strong>0</strong>): </p>
<pre>   der(d<strong>x</strong>) = der(<strong>f</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>f</strong>,<strong>u</strong>)*d<strong>u</strong> + <strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0)</pre>
<pre>        d<strong>y</strong> = der(<strong>g</strong>,<strong>x</strong>)*d<strong>x</strong> + der(<strong>g</strong>,<strong>u</strong>)*d<strong>u</strong> + (<strong>g</strong>(<strong>x</strong>0,<strong>u</strong>0) - <strong>y</strong>0)</pre>
<p>or </p>
<pre>   der(d<strong>x</strong>) = <strong>A</strong>*d<strong>x</strong> + <strong>B</strong>*d<strong>u</strong> + <strong>f</strong>0</pre>
<pre>        d<strong>y</strong> = <strong>C</strong>*d<strong>x</strong> + <strong>D</strong>*d<strong>u</strong></pre>
<p>This function returns the matrices <strong>A</strong>, <strong>B</strong>, <strong>C</strong>, <strong>D</strong> and assumes that the linearization point is a steady-state point of the simulation (i.e., <strong>f</strong>(<strong>x</strong>0,<strong>u</strong>0) = 0). Additionally, the full Modelica names of all inputs, outputs and states shall be returned if possible (default is to return empty name strings).</p>
</html>"));
end linearize;

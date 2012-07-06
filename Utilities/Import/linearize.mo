within Modelica_LinearSystems2.Utilities.Import;
function linearize "Linearize a model after simulation up to a given time"
  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica.SIunits.Time t_linearize= 0
    "Simulate until t_linearize and then linearize" annotation(Dialog);
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options it t_linearize > 0" annotation(Dialog);
protected
  String fileName="dslin";
  String fileName2=fileName+".mat";

  // Simulate until t_linearize and then linearize at this time instant
  Boolean OK1 = if t_linearize <= 0.0 then true else
                   simulateModel(problem=modelName, startTime=0, stopTime=t_linearize,
                                 method=simulationSetup.method,
                                 tolerance=simulationSetup.tolerance,
                                 fixedstepsize=simulationSetup.fixedStepSize);
  Boolean OK2 = if t_linearize <= 0.0 then true else importInitial("dsfinal.txt");
  Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

  // Read linear system from file
  Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
  Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
  Integer nx=integer(nxMat[1, 1]);
  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
  String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);
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
<p>
This function initializes a Modelica model and then simulates the model 
until time instant \"t_linearize\". 
If t_linearize=0, no simulation takes place (only initialization). 
At the simulation stop time, the model is linearized in such a form that 
</p>
<ul>
<li>all top-level signals with prefix \"input\" are treated as inputs <b>u</b>(t) of the model ,</li>
<li>all top-level signals with prefix \"output\" are treated as outputs <b>y</b>(t) of the model,</li>
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
<blockquote><pre>
der(<b>x</b>0+d<b>x</b>) = <b>f</b>(<b>x</b>0,<b>u</b>0) + der(<b>f</b>,<b>x</b>)*d<b>x</b> + der(<b>f</b>,<b>u</b>)*d<b>u</b>
    <b>y</b>0 + d<b>y</b> = <b>g</b>(<b>x</b>0,<b>u</b>0) + der(<b>g</b>,<b>x</b>)*d<b>x</b> + der(<b>g</b>,<b>u</b>)*d<b>u</b>
</pre></blockquote>
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
This function returns the matrices <b>A</b>, <b>B</b>, <b>C</b>, <b>D</b> and assumes that the linearization point is a steady-state point of the simulation (i.e., <b>f</b>(<b>x</b>0,<b>u</b>0) = 0). Additionally, the full Modelica names of all inputs, outputs and states shall be returned if possible (default is to return empty name strings).
</p>
</html>"));
end linearize;

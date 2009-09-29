within Modelica_LinearSystems2.Tests;
encapsulated function timeResponsePlots
  "Plot the time response of the system. The response type is selectable"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Types.TimeResponse;

  import Modelica_LinearSystems2.Utilities.Plot;

  extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Impulse-, Step-, ramp-response"));

  annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>timeResponse</b>(ss);
or
StateSpace.Plot.<b>timeResponse</b>(ss, dt, tSpan,response, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>


<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>timeResponse</b> plots the time response of a state space system. The character of the time response if defined by the input <tt>response</tt>, i.e. Impulse, Step, Ramp, or Initial. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.plotImpulse\">plotImpulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[0.0,1.0,1.0],
D=[0.0])

Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
Modelica_LinearSystems2.StateSpace.Plot.timeResponse(ss, response=response)
// gives:
</pre></blockquote>

</p>
<p align=\"center\">
<img src=\"../Extras/Images/timeResponseSS.png\">
</p>
<p>
</p>


</html> "));

protected
 StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; -1,-1],
      B=[1;0],
      C=[1,0],
      D=[0]);

  Real dt=0 "Sample time [s]";
  Real tSpan=0 "Simulation time span [s]";

  Modelica_LinearSystems2.Types.TimeResponse response1=
      Modelica_LinearSystems2.Types.TimeResponse.Impulse;
  Modelica_LinearSystems2.Types.TimeResponse response2=
      Modelica_LinearSystems2.Types.TimeResponse.Step;
  Modelica_LinearSystems2.Types.TimeResponse response3=
      Modelica_LinearSystems2.Types.TimeResponse.Ramp;

  input Real x0[2]=zeros(2) "Initial state vector";

  input Boolean subPlots=true
    "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));
  Plot.Records.Curve curve1;
  Plot.Records.Curve curve2;
  Plot.Records.Curve curve3;
  Integer i1;
  Integer i2;
  Plot.Records.Diagram diagram2[1];

  Real y1[:,size(ss.C, 1),size(ss.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  Real t1[:] "Time vector: (number of samples)";
  Real x1[:,size(ss.A, 1),size(ss.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";
  Real y2[:,size(ss.C, 1),size(ss.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  Real t2[:] "Time vector: (number of samples)";
  Real x2[:,size(ss.A, 1),size(ss.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";
  Real y3[:,size(ss.C, 1),size(ss.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  Real t3[:] "Time vector: (number of samples)";
  Real x3[:,size(ss.A, 1),size(ss.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";
  String yNames[size(ss.C, 1)];
  String uNames[size(ss.B, 2)];

algorithm
  (y1,t1,x1) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
    sc=ss,
    dt=dt,
    tSpan=tSpan,
    response=response1,
    x0=x0);

  (y2,t2,x2) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
    sc=ss,
    dt=dt,
    tSpan=tSpan,
    response=response2,
    x0=x0);

  (y3,t3,x3) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
    sc=ss,
    dt=dt,
    tSpan=tSpan,
    response=response3,
    x0=x0);

// generate headings

  for i2 in 1:size(ss.B, 2) loop
    for i1 in 1:size(ss.C, 1) loop
      curve1 := Plot.Records.Curve(
        x=t1,
        y=y1[:, i1, i2],
        autoLine=true);

      curve2 := Plot.Records.Curve(
        x=t2,
        y=y2[:, i1, i2],
        autoLine=true);

      curve3 := Plot.Records.Curve(
        x=t3,
        y=y3[:, i1, i2],
        autoLine=true);

      diagram2[i1] := defaultDiagram;
      diagram2[i1].curve := {curve1,curve2,curve3};
      diagram2[i1].heading := defaultDiagram.heading + "  " + uNames[i2] + " -> " + yNames[i1];
      diagram2[i1].yLabel := yNames[i1];

    end for;

    if subPlots then
      Plot.diagramVector(diagram2, device);
    else
      for i1 in 1:size(ss.C, 1) loop
        Plot.diagram(diagram2[i1], device);
      end for;
    end if;
  end for;

end timeResponsePlots;

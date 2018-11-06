within Modelica_LinearSystems2.Examples.StateSpace;
function analysisStepResponse "Step response example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;

protected
  Modelica_LinearSystems2.StateSpace sc=Modelica_LinearSystems2.StateSpace(
    A=[-1,1; 0,-2],
    B=[1,0; 0,1],
    C=[1,0; 0,1],
    D=[0,0; 0,0]);

  Real t[:] "Time vector: (number of samples)";
  Real x_continuous[:,size(sc.A, 1),size(sc.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

public
  output Real y[:,size(sc.C, 1),size(sc.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inuputs)";

algorithm
  (y,t,x_continuous) :=
    Modelica_LinearSystems2.StateSpace.Analysis.stepResponse(
      sc=sc,
      dt=0.1,
      tSpan=5);

  Modelica_LinearSystems2.Utilities.Plot.diagramVector({
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y[:,1,1],
          legend="y1"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y[:,2,1],
          legend="y2")},
      heading="Step response to u1",
      xLabel="time [s]",
      yLabel="y1, y2"),
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y[:,1,2],
          legend="y1"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y[:,2,2],
          legend="y2")},
      heading="Step response to u2",
      xLabel="time [s]",
      yLabel="y1, y2")});

  annotation (
    Documentation(info="<html>
<p>
In this example, the <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.stepResponse\">step response</a>
of a state space system is computed and plotted.
</p>
</html>"));
end analysisStepResponse;

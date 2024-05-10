within Modelica_LinearSystems2.Examples.StateSpace;
function analysisImpulseResponse "Impulse response example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Utilities.Plot;

protected
  StateSpace sc = StateSpace(
    A=[-1,1; 0,-2],
    B=[4,0; 0,2],
    C=[2,0; 0,3],
    D=[0,0; 0,0]);

  Real t[:] "Time vector: (number of samples)";
  Real x_continuous[:,size(sc.A, 1),size(sc.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

public
  output Real y[:,size(sc.C, 1),size(sc.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inuputs)";

algorithm
  (y,t,x_continuous) := StateSpace.Analysis.impulseResponse(sc=sc,dt=0.1,tSpan=5);

  Plot.diagramVector({
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(x=t, y=y[:,1,1], legend="y1"),
        Plot.Records.Curve(x=t, y=y[:,2,1], legend="y2")},
      heading="Impulse response to u1",
      xLabel="time [s]",
      yLabel="y1, y2"),
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(x=t, y=y[:,1,2], legend="y1"),
        Plot.Records.Curve(x=t, y=y[:,2,2], legend="y2")},
      heading="Impulse response to u2",
      xLabel="time [s]",
      yLabel="y1, y2")});

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Computes and plots the impulse response of a state space system.
</p>
</html>"));
end analysisImpulseResponse;

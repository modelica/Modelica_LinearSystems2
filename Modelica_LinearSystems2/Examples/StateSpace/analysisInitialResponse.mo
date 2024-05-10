within Modelica_LinearSystems2.Examples.StateSpace;
function analysisInitialResponse "Initial response example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Utilities.Plot;

protected
  StateSpace sc=StateSpace(
    A=[-1,1; 0,-2],
    B=[1,0; 0,1],
    C=[1,0; 0,1],
    D=[0,0; 0,0]);

  Real t[:] "Time vector: (number of samples)";
  Real x_continuous[:,2,1]
    "State trajectories: (number of samples) x (number of states) x 1";
  Real x0[2]=ones(2) "Initial state vector";

public
  output Real y[:,size(sc.C, 1),1]
    "Output response: (number of samples) x (number of outputs) x 1";

algorithm
  (y,t,x_continuous) := StateSpace.Analysis.initialResponse(x0=x0,sc=sc,dt=0.1,tSpan=5);

  Plot.diagramVector({
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(
          x=t,
          y=y[:,1,1],
          legend="y1")},
      heading="Initial response y1",
      xLabel="time [s]"),
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(
          x=t,
          y=y[:,2,1],
          legend="y2")},
      heading="Initial response y2",
      xLabel="time [s]")});

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Computes and plots the initial response of a state space system by given initial values.
</p>
</html>"));
end analysisInitialResponse;

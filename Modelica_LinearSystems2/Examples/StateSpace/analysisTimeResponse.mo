within Modelica_LinearSystems2.Examples.StateSpace;
function analysisTimeResponse "Compute time response of a state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
  import Modelica_LinearSystems2.Utilities.Plot;

  input StateSpace sc=StateSpace(
    A=[-1,1; 0,-2],
    B=[1,0; 0,1],
    C=[1,0; 0,1],
    D=[0,0; 0,0]);
  input TimeResponse response=TimeResponse.Step;
  input Real dt=0.01;
  input Real span=5;

//   Modelica_LinearSystems2.DiscreteStateSpace sd=
//       Modelica_LinearSystems2.DiscreteStateSpace(sc, Ts);
protected
  Real x0[2]={0,0};
  Integer samples = integer(span/dt+1);

public
  output Real y[samples,2,2]
    "System response (dimension: (input samples) x (number of outputs))";
  output Real t[samples] "Time vector used for simulation";
  output Real x[samples,2,2]
    "State trajectories (dimension: (input samples) x (number of states)";
algorithm
  // (y,x) := Modelica_LinearSystems2.DiscreteStateSpace.timeResponse(
  //   sd,
  //   u,
  //   x0);

  (y,t,x) := StateSpace.Analysis.timeResponse(response,sc,dt, span);

  // Plot.diagram(
  //   Plot.Records.Diagram(
  //     curve={
  //       Plot.Records.Curve(x=t, y=y[:, 1, 1], legend="y1"),
  //       Plot.Records.Curve(x=t, y=y[:, 2, 2], legend="y2")},
  //     heading="Step response to synchronous step of u1 and u2",
  //     xLabel="time [s]",
  //     yLabel="y1, y2"));

  Plot.diagram(
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(x=t, y=y[:, 1, 1], legend="y1_1"),
        Plot.Records.Curve(x=t, y=y[:, 1, 2], legend="y1_2"),
        Plot.Records.Curve(x=t, y=y[:, 2, 1], legend="y2_1"),
        Plot.Records.Curve(x=t, y=y[:, 2, 2], legend="y2_2")},
      heading="Step response to synchronous step of u1 and u2",
      xLabel="time [s]",
      yLabel="y1, y2"));
  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Computes the time response of the system
StateSpace <em>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</em>,
sampled at <em>Ts=0.01</em> with inititial state <em>x0=[0;0]</em>
subject to the system input <em>u = ones(samples,2)</em>, (<em>samples</em> is set to 30).
</p>
</html>"));
end analysisTimeResponse;

within Modelica_LinearSystems2.Examples.StateSpace;
function analysisTimeResponse "Compute time response of a state space system"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

  input Modelica_LinearSystems2.StateSpace sc=Modelica_LinearSystems2.StateSpace(
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

  (y,t,x) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(response,sc,dt, span);

  // Modelica_LinearSystems2.Utilities.Plot.diagram(
  //   Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
  //   curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
  //     x=t,
  //     y=y[:, 1, 1],
  //     legend="y1"),Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
  //     x=t,
  //     y=y[:, 2, 2],
  //     legend="y2")},
  //   heading="Step response to synchronous step of u1 and u2",
  //   xLabel="time [s]",
  //   yLabel="y1, y2"));

  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
    curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=t,
      y=y[:, 1, 1],
      legend="y1_1"),Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=t,
      y=y[:, 1, 2],
      legend="y1_2"), Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=t,
      y=y[:, 2, 1],
      legend="y2_1"),Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=t,
      y=y[:, 2, 2],
      legend="y2_2")},
    heading="Step response to synchronous step of u1 and u2",
    xLabel="time [s]",
    yLabel="y1, y2"));
  annotation (
    Documentation(info="<html>
<p>
Computes the time response of the system
StateSpace <i>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</i>,
sampled at <i>Ts=0.01</i> with inititial state <i>x0=[0;0]</i>
subject to the system input <i>u = ones(samples,2)</i>, (<i>samples</i> is set to 30).
</p>
</html>"));
end analysisTimeResponse;

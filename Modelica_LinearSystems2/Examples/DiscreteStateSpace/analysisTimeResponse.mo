within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function analysisTimeResponse
  "Compute time response of a discrete state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Utilities.Plot;

  input Real u[:,2]=ones(300, 2);
protected
  DiscreteStateSpace dss = DiscreteStateSpace(
    A=[0.99005, 0.00985;
       0,       0.9802],
    B=[0.01,0;
       0,0.01],
    C=[1,0; 0,1],
    D=[0,0; 0,0],
    Ts=0.01,
    B2=[0,0; 0,0],
    method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact);

  Integer samples=size(u,1);
  Real x0[2]={0, 0};
public
  output Real y[samples,2]
    "System response (dimension: (input samples) x (number of outputs))";
  output Real t[samples] "Time vector used for simulation";
  output Real x[samples,2]
    "State trajectories (dimension: (input samples) x (number of states)";
algorithm
  t := 0:dss.Ts:(samples*dss.Ts - dss.Ts);
  (y,x) := DiscreteStateSpace.timeResponse(dss, u, x0);

  Plot.diagram(
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(x=t, y=y[:,1], legend="y1"),
        Plot.Records.Curve(x=t, y=y[:,2], legend="y2")},
      heading="Step response to synchronous step of u1 and u2",
      xLabel="time [s]",
      yLabel="y1, y2"));

  annotation (Documentation(info="<html>
<p>
Computes the time response of the system
</p>
<blockquote><pre>
dss = DiscreteStateSpace(
        A=[0.99005, 0.00985;
           0,       0.9802],
        B=[0.01,0;
           0,0.01],
        C=[1,0; 0,1],
        D=[0,0; 0,0],
        Ts=0.01,
        B2=[0,0; 0,0],
        method=Modelica_LinearSystems2.Types.Method.StepExact),
</pre></blockquote>
<p>
sampled at <em>Ts=0.01</em> with inititial state <em>x0=[0;0]</em>
subject to the system input <em>u = ones(300,2)</em> (which results with Ts=0.01 in 3 sec).
</p>
</html>"), __Dymola_interactive=true);
end analysisTimeResponse;

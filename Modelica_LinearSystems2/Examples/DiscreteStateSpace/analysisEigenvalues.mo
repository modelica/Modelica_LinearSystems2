within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function analysisEigenvalues
  "Example to compute the eigenvalues of a discrete state space system"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input StateSpace ss = StateSpace(
    A=[-1,1; -1,-1],
    B=[1; 1],
    C=[1,1],
    D=[0]);

  input Modelica.SIunits.Time Ts=0.1 "Sample time";
  input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact
    "Discretization method";
protected
  DiscreteStateSpace dss = DiscreteStateSpace(
    ss,
    Ts=0.1,
    method=method);

  Complex evDiscrete[:] = DiscreteStateSpace.Analysis.eigenValues(dss);
  Complex evContinuous[:] = StateSpace.Analysis.eigenValues(ss);

   //alternative calculation
  Complex ev1=Complex.exp(evContinuous[1]*Ts);
  Complex ev2=Complex.exp(evContinuous[2]*Ts);
  Complex evDiscrete2[2]={ev1,ev2};

algorithm
  Complex.Vectors.print("evDiscrete", evDiscrete);
  Complex.Vectors.print("evDiscrete2", evDiscrete2);
  Complex.Vectors.print("evContiuous", evContinuous);

  annotation (Documentation(info="<html>
<p>
This example shows the computation of the poles and zeros of state space system.
</p>
</html>"));
end analysisEigenvalues;

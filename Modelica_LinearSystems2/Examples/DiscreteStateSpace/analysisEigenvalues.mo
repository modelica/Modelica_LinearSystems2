within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function analysisEigenvalues
  "Compute the eigenvalues of a discrete state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.ComplexMathAdds.Vectors;

  input StateSpace ss = StateSpace(
    A=[-1,1; -1,-1],
    B=[1; 1],
    C=[1,1],
    D=[0]);

  input Modelica.Units.SI.Time Ts=0.1 "Sample time";
  input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact "Discretization method";
protected
  DiscreteStateSpace dss = DiscreteStateSpace(
    ss,
    Ts=Ts,
    method=method);

  Complex evDiscrete[:] = DiscreteStateSpace.Analysis.eigenValues(dss);
  Complex evContinuous[:] = StateSpace.Analysis.eigenValues(ss);

   //alternative calculation
  Complex ev1=Modelica.ComplexMath.exp(evContinuous[1]*Ts);
  Complex ev2=Modelica.ComplexMath.exp(evContinuous[2]*Ts);
  Complex evDiscrete2[2]={ev1,ev2};

algorithm
  Vectors.print("evDiscrete", evDiscrete);
  Vectors.print("evDiscrete2", evDiscrete2);
  Vectors.print("evContiuous", evContinuous);

  annotation (Documentation(info="<html>
<p>
This example computes a&nbsp;discrete state space system from a&nbsp;continuous state space system.
Consequently, the eigenvalues of both systems are computed and printed.
</p>
</html>"));
end analysisEigenvalues;

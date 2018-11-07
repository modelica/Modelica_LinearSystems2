within Modelica_LinearSystems2.Controllers;
block StateSpace "Continuous or discrete state space system block"

  parameter Modelica_LinearSystems2.StateSpace system = Modelica_LinearSystems2.StateSpace(
    A=fill(0, 0, 0),
    B=fill(0, 0, 1),
    C=fill(0, 1, 0),
    D=fill(0, 1, 1)) "Continuous linear time-invariant system"
    annotation(Dialog(enable=continuous));
  extends Interfaces.PartialSampledBlock;
  parameter Real x_start[nx]=zeros(nx) "Initial or guess values of states"
    annotation(Dialog(tab="Advanced options"));
  parameter Real y_start[ny]=zeros(ny)
    "Initial values of outputs (remaining states are in steady state if possible)"   annotation(Dialog(tab="Advanced options"));
  Modelica.Blocks.Interfaces.RealInput u[size(system.B, 2)]
    "Continuous or discrete input signals of block"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y[size(system.C, 1)]
    "Continuous or discrete output signals of block"
    annotation (Placement(transformation(extent={{100,-10},{120,10}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput x[nx](start=x_start)
    "State vector of continuous system";
  final parameter Integer nx=size(system.A, 1) "Number of states x" annotation(HideResult=true);
  final parameter Integer ny=size(system.C, 1) "Number of outputs y" annotation(HideResult=true);

  parameter Boolean withDelay=false
    "True, if a unit delay should be considered"
    annotation(Evaluate=true, HideResult=true,Dialog(tab="Advanced options"));
protected
  Internal.DiscreteStateSpace discretePart(
    system=system,
    withDelay=withDelay,
    methodType=methodType,
    sampleFactor=sampleFactor,
    init=init,
    x_start=x_start,
    y_start=y_start) if not continuous "Discretized state space system";

equation
  if continuous then
    der(x) = system.A*x + system.B*u;
    y = system.C*x + system.D*u;
  end if;

  connect(u, discretePart.u);
  connect(x, discretePart.x);
  connect(y, discretePart.y);

initial equation
  if continuous then
    if init == Types.Init.InitialState then
      x = x_start;
    elseif init == Types.Init.SteadyState then
      der(x) = zeros(nx);
    elseif init == Types.Init.InitialOutput and nx>0 then
      x = Modelica.Math.Matrices.equalityLeastSquares(system.A, -system.B*u, system.C, y_start - system.D*u);
//      y = y_start;
//      der(x[ny + 1:nx]) = zeros(nx - ny);
    end if;
  end if;

  annotation (
    defaultComponentName="stateSpace",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Text(
          extent={{0,0},{-100,100}},
          lineColor={0,0,127},
          textString="A"),
        Text(
          extent={{100,0},{0,100}},
          lineColor={0,0,127},
          textString="B"),
        Text(
          extent={{-100,0},{0,-100}},
          lineColor={0,0,127},
          textString="C"),
        Text(
          extent={{0,0},{100,-100}},
          lineColor={0,0,127},
          textString="D"),
        Text(
          extent={{-96,18},{100,-16}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}),
        Text(
          extent={{-56,40},{60,0}},
          lineColor={0,0,0},
          textString="x=Ax+Bu"),
        Text(
          extent={{-60,0},{60,-40}},
          lineColor={0,0,0},
          textString=" y=Cx+Du"),
        Line(points={{-100,0},{-60,0}}),
        Line(points={{60,0},{100,0}})}),
    Documentation(info="<html>
</html>"));
end StateSpace;

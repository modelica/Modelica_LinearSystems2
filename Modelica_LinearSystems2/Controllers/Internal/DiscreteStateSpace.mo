within Modelica_LinearSystems2.Controllers.Internal;
model DiscreteStateSpace
  "Discrete linear time invariant state space system block parameterized with a continuous state space system"
  extends Controllers.Icons.PartialBlockIcon(cont=false);

  import Modelica_LinearSystems2.Controllers.Types;
  import Modelica.Math.Matrices;

  parameter Modelica_LinearSystems2.StateSpace system = Modelica_LinearSystems2.StateSpace(
    A=fill(0, 0, 0),
    B=fill(0, 0, 1),
    C=fill(0, 1, 0),
    D=fill(0, 1, 1), uNames={"u"},yNames={"y"})
    "Continuous linear time-invariant system" annotation(HideResult=true);

  parameter Types.MethodWithGlobalDefault methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.UseSampleClockOption
    "Type of discretization" annotation(Evaluate=true, HideResult=true);
  final parameter Types.Method method=convertToMethod(methodType, sampleClock.methodType)
    "Discretization method (explicitEuler/implicitEuler/trapezoidal/stepExact/rampExact)"
    annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1) = 1
    "Factor so that sample time = sampleClock.sampleTime * sampleFactor"
    annotation(HideResult=true);

  parameter Types.Init init = Modelica_LinearSystems2.Controllers.Types.Init.SteadyState
    "Type of initialization (No init/SteadyState/InitialState/InitialOutput)"
    annotation(Evaluate=true, HideResult=true);

  parameter Real x_start[nx]=zeros(nx)
    "Initial value of continuous state x, if init=InitialState (otherwise guess value)"
    annotation(Evaluate=true,HideResult=true);
  parameter Real y_start[ny]=zeros(ny)
    "Initial value of continuous output y, if init=InitialOutput (otherwise guess value)"
    annotation(Evaluate=true,HideResult=true);
  parameter Boolean withDelay = false
    "True, if a unit delay should be considered";

  final parameter Modelica.SIunits.Time Ts=sampleClock.sampleTime*sampleFactor
    "Sample time" annotation(HideResult=false);
  final parameter Integer nx=size(system.A, 1) "Number of states"  annotation(HideResult=true);
  final parameter Integer nu=size(system.B, 2) "Number of inputs"  annotation(HideResult=true);
  final parameter Integer ny=size(system.C, 1) "Number of outputs"  annotation(HideResult=true);

  Modelica.Blocks.Interfaces.RealInput u[nu]
    "Continuous or discrete input signals of block"
    annotation (HideResult=true, Placement(transformation(
          extent={{-140,-20},{-100,20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y[ny](start=y_start)
    "Discrete output signals of block" annotation (HideResult=true,
      Placement(transformation(extent={{100,-10},{120,10}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput x[nx](start=x_start)
    "State vector of continuous system at sample times" annotation(HideResult=true);

protected
  outer SampleClock sampleClock "Global options" annotation(HideResult=true);
  parameter Modelica_LinearSystems2.DiscreteStateSpace discreteSystem=
    Modelica_LinearSystems2.DiscreteStateSpace(system, Ts, method);

  discrete Real xd[nx](start=x_start)
    "State vector of discrete system (pre(xd) = x - B2*u)";
  discrete Real new_xd[nx](start=x_start) "Next valued of xd"
    annotation(HideResult=true);

// Derived quantities
  discrete Real u_sampled[nu] "Sampled continuous input signal u";
  discrete Real pre_u_sampled[nu] "Sampled continuous input signal u";
  discrete Real y_sampled[ny] "Sampled continuous output"        annotation(HideResult=true);
  discrete Real x_sampled[nx] "Sampled continuous state"        annotation(HideResult=true);
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(HideResult=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(HideResult=true);

equation
  if sampleClock.blockType == Types.BlockType.Continuous then
    // no sampling in sampleClock
    sampleTrigger = sample(Ts, Ts);
    ticks = 0;
  else
    when sampleClock.sampleTrigger then
      ticks = if pre(ticks) < sampleFactor then pre(ticks) + 1 else 1;
    end when;
    sampleTrigger = sampleClock.sampleTrigger and ticks >= sampleFactor;
  end if;

  when {initial(),sampleTrigger} then
    u_sampled = u;
    pre_u_sampled = pre(u_sampled);
    if withDelay then
      new_xd = discreteSystem.B*pre_u_sampled + discreteSystem.A*xd;
      y_sampled = discreteSystem.C*xd + discreteSystem.D*pre_u_sampled;
      x_sampled = xd + discreteSystem.B2*pre_u_sampled;
    else
      new_xd = discreteSystem.B*u_sampled + discreteSystem.A*xd;
      y_sampled = discreteSystem.C*xd + discreteSystem.D*u_sampled;
      x_sampled = xd + discreteSystem.B2*u_sampled;
    end if;
    xd = pre(new_xd);
  end when;

  y = y_sampled;
  x = x_sampled;

initial equation
  pre(ticks) = 0;

  if init == Types.Init.InitialState then
    x = x_start;

  elseif init == Types.Init.SteadyState then
    if Matrices.isEqual(
        discreteSystem.A,
        identity(nx),
        100*Modelica.Constants.eps) then
        // block contains an integrator and is only possible to initialize with steady states steady state when u==0
      u = fill(0.0, nu);
      // xd[nu + 1:nx] = pre(x[nu + 1:nx]);
      xd[nu + 1:nx] = x[nu + 1:nx];

    else
      xd = discreteSystem.B*u_sampled + discreteSystem.A*xd;//new_xd;

    end if;

  elseif init == Types.Init.InitialOutput and nx>0 then
    // y=y_start;
    // xd[ny+1:nx]=[zeros(nx-ny,ny),identity(nx-ny)]*Modelica.Math.Matrices.solve(identity(nx)-discreteSystem.A,discreteSystem.B*u);
    xd = Modelica.Math.Matrices.equalityLeastSquares(identity(nx)-discreteSystem.A, -discreteSystem.B*u, discreteSystem.C, y_start - discreteSystem.D*u);
  end if;

  annotation (
    defaultComponentName="discreteStateSpace",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Text(extent={{-90,15},{-15,90}}, textString="A"),
        Text(extent={{15,15},{90,90}}, textString="B"),
        Text(extent={{-90,-15},{-15,-90}}, textString="C"),
        Text(extent={{15,-15},{90,-90}}, textString="D")}),
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
end DiscreteStateSpace;

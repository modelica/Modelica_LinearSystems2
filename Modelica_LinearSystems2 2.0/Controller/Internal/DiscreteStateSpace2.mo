within Modelica_LinearSystems2.Controller.Internal;
model DiscreteStateSpace2
  "Discrete linear time invariant state space system block parameterized with a continuous state space system"

  import Modelica_LinearSystems2.Controller.Types;
  import Modelica.Math.Matrices;

  extends Interfaces.PartialBlockIcon;

  parameter Real ABCD[:,:] "Continuous linear time-invariant system"
                                              annotation(Hide=true);

protected
      parameter Real A[:,:] = ABCD[1:end-1,1:end-1];
      parameter Real B[:,:] = ABCD[1:end-1,end:end];
      parameter Real C[:,:] = ABCD[end:end,1:end-1];
      parameter Real D[:,:] = matrix(ABCD[end,end]);
public
  final parameter Integer nx=size(A, 1) "Number of states"  annotation(Hide=true);
  final parameter Integer nu=size(B, 2) "Number of inputs"  annotation(Hide=true);
  final parameter Integer ny=size(C, 1) "Number of outputs"  annotation(Hide=true);
  parameter Types.MethodWithGlobalDefault methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.UseSampleClockOption
    "Type of discretization" annotation(Evaluate=true, Hide=true);
  final parameter Types.Method method=if methodType == Types.MethodWithGlobalDefault.UseSampleClockOption then
            sampleClock.methodType else methodType
    "Discretization method (explicitEuler/implicitEuler/trapezoidal/stepExact/rampExact)"
     annotation(Evaluate=true, Hide=false);

  parameter Integer sampleFactor(min=1) = 1
    "sample time=sampleClock.sampleTime*sampleFactor"
     annotation(Hide=true);
  parameter Types.Init init=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (No init/InitialState/SteadyState/Output)"
    annotation(Evaluate=true, Hide=true);
  parameter Real x_start[nx]=zeros(nx)
    "Initial value of continuous state x, if init=InitialState (otherwise guess value)"
    annotation(Evaluate=true,Hide=true);
  parameter Real y_start[ny]=zeros(ny)
    "Initial value of continuous output y, if init=InitialOutput (otherwise guess value)"
    annotation(Evaluate=true,Hide=true);
  parameter Boolean withDelay = false
    "is true if a unit delay should be considered";

  final parameter Modelica.SIunits.Time Ts=sampleClock.sampleTime*sampleFactor
    "Sample time" annotation(Hide=false);

  Modelica.Blocks.Interfaces.RealInput u[nu]
    "Continuous or discrete input signals of block"
    annotation (                          Hide=true, Placement(transformation(
          extent={{-140,-20},{-100,20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y[ny](start=y_start)
    "Discrete output signals of block" annotation (                        Hide=true,
      Placement(transformation(extent={{100,-10},{120,10}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput x[nx](start=x_start)
    "State vector of continuous system at sample times"                                                        annotation(Hide=true);
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
    Window(
      x=0.27,
      y=0.13,
      width=0.55,
      height=0.75),
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
    Documentation(info="<HTML>
</HTML>"));

protected
  outer SampleClock sampleClock "Global options"                       annotation(Hide=true);
  parameter Modelica_LinearSystems2.DiscreteStateSpace discreteSystem=
      Modelica_LinearSystems2.DiscreteStateSpace(
      A, B, C, D,
      Ts,
      method);

 discrete Real xd[nx](start=x_start)
    "State vector of discrete system (pre(xd) = x - B2*u)";
  discrete Real new_xd[nx](start=x_start) "Next valued of xd"
                        annotation(Hide=true);

// Derived quantities
  discrete Real u_sampled[nu] "Sampled continuous input signal u";
  discrete Real pre_u_sampled[nu] "Sampled continuous input signal u";
  discrete Real y_sampled[ny] "Sampled continuous output"        annotation(Hide=true);
  discrete Real x_sampled[nx] "Sampled continuous state"        annotation(Hide=true);
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(Hide=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(Hide=true);

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
      xd[nu + 1:nx] = pre(x[nu + 1:nx]);
    else

      xd = discreteSystem.B*u_sampled + discreteSystem.A*xd;//new_xd;

    end if;

  elseif init == Types.Init.InitialOutput then
  y=y_start;
  xd[ny+1:nx]=[zeros(nx-ny,ny),identity(nx-ny)]*Modelica.Math.Matrices.solve(identity(nx)-discreteSystem.A,discreteSystem.B*u);

  end if;
end DiscreteStateSpace2;

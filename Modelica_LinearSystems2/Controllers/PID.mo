within Modelica_LinearSystems2.Controllers;
block PID "PID-controller in additive description form"
  extends Interfaces.PartialSampledBlock;
  import Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault;

  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
    annotation(Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
    annotation(Placement(transformation(extent={{100,-10},{120,10}})));

  parameter Types.PID_representation pidRep=Types.PID_representation.timeConstants
    "Type of PID representation";

  parameter Real k(min=0) = 1 "Gain of controller" annotation(Dialog(enable=pidRep==Types.PID_representation.timeConstants));
  parameter Modelica.SIunits.Time Ti(min=Modelica.Constants.small, start=0.5)
    "Time constant of Integrator block"
    annotation(Dialog(enable=pidRep==Types.PID_representation.timeConstants));
  parameter Modelica.SIunits.Time Td(min=0, start=0.1)
    "Time constant of Derivative block"
    annotation(Dialog(enable=pidRep==Types.PID_representation.timeConstants));
  parameter Real Nd(min=Modelica.Constants.small) = 10
    "The higher Nd, the more ideal the derivative block";

  parameter Real kp = 1 "P part parameter of gain representation"
    annotation(Dialog(enable=pidRep==Types.PID_representation.gains));
  parameter Real ki= 1 "I part parameter of gain representation"
    annotation(Dialog(enable=pidRep==Types.PID_representation.gains));
  parameter Real kd = 1 "D part parameter of gain representation"
    annotation(Dialog(enable=pidRep==Types.PID_representation.gains));

  parameter Real xi_start=0
    "Initial or guess value value for integrator output (= integrator state)"
    annotation (Dialog(group="Initialization"));
  parameter Real xd_start=0
    "Initial or guess value for state of derivative block"
    annotation (Dialog(group="Initialization"));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(
        enable = initType == Types.InitWithGlobalDefault.InitialOutput,
        group = "Initialization"));

  Sampler sampler(blockType=blockType)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  Modelica.Blocks.Math.Gain P(k=1)
    annotation (Placement(transformation(extent={{-30,30},{-10,50}}, rotation=0)));
  Derivative D(
    k= if pidRep==Types.PID_representation.timeConstants then Td else kd/kp,
    blockType=blockType,
    initType=if init == Types.Init.SteadyState or init == Types.Init.InitialOutput then
              InitWithGlobalDefault.SteadyState else if init == Types.Init.InitialState then
              InitWithGlobalDefault.InitialState else InitWithGlobalDefault.NoInit,
    x_start=xd_start,
    y_start=y_start,
    methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.StepExact,
    T=max([if pidRep==Types.PID_representation.timeConstants then Td/Nd else kd/kp/Nd,100*Modelica.Constants.eps]))
    annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

  Integrator I(
    k=if pidRep==Types.PID_representation.timeConstants then 1/Ti else kp/ki,
    y_start=xi_start,
    blockType=blockType,
    initType=if init == Types.Init.SteadyState then InitWithGlobalDefault.SteadyState else
              if init == Types.Init.InitialState then InitWithGlobalDefault.InitialState else
              InitWithGlobalDefault.NoInit,
    methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.StepExact)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

  Modelica.Blocks.Math.Add3 addPID
    annotation (Evaluate=true, Placement(transformation(
          extent={{10,-10},{30,10}},rotation=0)));
  Modelica.Blocks.Math.Gain gainPID(k=if pidRep==Types.PID_representation.timeConstants then k else kp)
    annotation (Placement(transformation(extent={{50,-10},
            {70,10}},          rotation=0)));

initial equation
  if init == Modelica_LinearSystems2.Controllers.Types.Init.InitialOutput then
    y = y_start;
  end if;

equation
  connect(P.y, addPID.u1) annotation (Line(
      points={{-9,40},{0,40},{0,8},{8,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gainPID.u, addPID.y) annotation (Line(
      points={{48,0},{31,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.u, u) annotation (Line(
      points={{-102,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(I.y, addPID.u2) annotation (Line(
      points={{-9,0},{8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(D.y, addPID.u3) annotation (Line(
      points={{-9,-50},{0,-50},{0,-8},{8,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.y, I.u) annotation (Line(
      points={{-79,0},{-32,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P.u, sampler.y) annotation (Line(
      points={{-32,40},{-60,40},{-60,0},{-79,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(D.u, sampler.y) annotation (Line(
      points={{-32,-50},{-60,-50},{-60,0},{-79,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gainPID.y, y) annotation (Line(
      points={{71,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    defaultComponentName="pID",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-80,50},{-73,-20},{30,60},{80,60}}, color={0,0,
              127}),
        Text(
          extent={{-46,-21},{80,-60}},
          lineColor={192,192,192},
          textString="PID"),
        Text(
          extent={{-102,85},{74,53}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Documentation(info="<html>
<p>
This is the text-book version of a PID-controller.
For a more practically useful PID-controller, use
block LimPID.
</p>

<p>
The PID block can be initialized in different
ways controlled by parameter <b>initType</b>. The possible
values of initType are defined in
<a href=\"modelica://Modelica.Blocks.Types.InitPID\">Modelica.Blocks.Types.InitPID</a>.
This type is identical to
<a href=\"modelica://Modelica.Blocks.Types.Init\">Types.Init</a>,
with the only exception that the additional option
<b>DoNotUse_InitialIntegratorState</b> is added for
backward compatibility reasons (= integrator is initialized with
InitialState whereas differential part is initialized with
NoInit which was the initialization in version 2.2 of the Modelica
standard library).
</p>

<p>
Based on the setting of initType, the integrator (I) and derivative (D)
blocks inside the PID controller are initialized according to the following table:
</p>

<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>initType</b></td>
      <td valign=\"top\"><b>I.initType</b></td>
      <td valign=\"top\"><b>D.initType</b></td></tr>

  <tr><td valign=\"top\"><b>NoInit</b></td>
      <td valign=\"top\">NoInit</td>
      <td valign=\"top\">NoInit</td></tr>

  <tr><td valign=\"top\"><b>SteadyState</b></td>
      <td valign=\"top\">SteadyState</td>
      <td valign=\"top\">SteadyState</td></tr>

  <tr><td valign=\"top\"><b>InitialState</b></td>
      <td valign=\"top\">InitialState</td>
      <td valign=\"top\">InitialState</td></tr>

  <tr><td valign=\"top\"><b>InitialOutput</b><br>
          and initial equation: y = y_start</td>
      <td valign=\"top\">NoInit</td>
      <td valign=\"top\">SteadyState</td></tr>

  <tr><td valign=\"top\"><b>DoNotUse_InitialIntegratorState</b></td>
      <td valign=\"top\">InitialState</td>
      <td valign=\"top\">NoInit</td></tr>
</table>

<p>
In many cases, the most useful initial condition is
<b>SteadyState</b> because initial transients are then no longer
present. If initType = InitPID.SteadyState, then in some
cases difficulties might occur. The reason is the
equation of the integrator:
</p>

<pre>
   <b>der</b>(y) = k*u;
</pre>

<p>
The steady state equation &quot;der(x)=0&quot; leads to the condition that the input u to the
integrator is zero. If the input u is already (directly or indirectly) defined
by another initial condition, then the initialization problem is <b>singular</b>
(has none or infinitely many solutions). This situation occurs often
for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and
since speed is both a state and a derivative, it is natural to
initialize it with zero. As sketched this is, however, not possible.
The solution is to not initialize u or the variable that is used
to compute u by an algebraic equation.
</p>
</html>"));
end PID;

within Modelica_LinearSystems2.Controllers;
block LimPID
  "P, PI, PD, and PID controller with limited output, anti-windup compensation and setpoint weighting"
  extends Modelica.Blocks.Interfaces.SVcontrol;
  extends Interfaces.PartialSampledBlock;

  import Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault;
  import Modelica.Blocks.Types.SimpleController;

  output Real controlError=u_s - u_m "Control error (set point - measurement)";

  parameter Types.PID_representation pidRepresentation=Types.PID_representation.timeConstants
    "Type of PID representation";
  parameter Modelica.Blocks.Types.SimpleController controllerType=Modelica.Blocks.Types.SimpleController.PID
    "Type of controller";
  parameter Real k(min=0) = 1 "Gain of controller" annotation(Dialog(enable=pidRepresentation==Types.PID_representation.timeConstants));
  parameter Modelica.SIunits.Time Ti(min=Modelica.Constants.small, start=0.5)
    "Time constant of Integrator block"
    annotation(Dialog(enable=pidRepresentation==Types.PID_representation.timeConstants and
                        (controllerType==SimpleController.PI or
                         controllerType==SimpleController.PID)));
  parameter Modelica.SIunits.Time Td(min=0, start=0.1)
    "Time constant of Derivative block"
    annotation(Dialog(enable=pidRepresentation==Types.PID_representation.timeConstants and (controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID)));
  parameter Real kp=1 "P part parameter of gain representation" annotation(Dialog(enable=pidRepresentation==Types.PID_representation.gains));
  parameter Real ki=1 "I part parameter of gain representation" annotation(Dialog(enable=pidRepresentation==Types.PID_representation.gains and   (controllerType==SimpleController.PI or  controllerType==SimpleController.PID)));
  parameter Real kd=1 "D part parameter of gain representation" annotation(Dialog(enable=pidRepresentation==Types.PID_representation.gains and (controllerType==SimpleController.PD or controllerType==SimpleController.PID)));
  parameter Real yMax(start=1) "Upper limit of output";
  parameter Real yMin=-yMax "Lower limit of output";
  parameter Real wp(min=0) = 1 "Set-point weight for Proportional block (0..1)";
  parameter Real wd(min=0) = 0 "Set-point weight for Derivative block (0..1)"
    annotation(Dialog(enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Real Ni(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation"
     annotation(Dialog(enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));
  parameter Real Nd(min=100*Modelica.Constants.eps) = 10
    "The higher Nd, the more ideal the derivative block"
    annotation(Dialog(enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Boolean limitsAtInit=true
    "= false, if limits are ignored during initializiation"
    annotation(Evaluate=true, Dialog(group="Initialization",
                       enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));
  parameter Real xi_start=0
    "Initial or guess value value for integrator output (= integrator state)"
    annotation (Dialog(group="Initialization",
                enable=controllerType==SimpleController.PI or
                       controllerType==SimpleController.PID));
  parameter Real xd_start=0
    "Initial or guess value for state of derivative block"
    annotation (Dialog(group="Initialization",
                         enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(
      enable = initType == Types.InitWithGlobalDefault.InitialOutput,
      group = "Initialization"));

  Sampler sampler_s(blockType=blockType)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  Sampler sampler_m(blockType=blockType)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-90})));
  Modelica.Blocks.Math.Add addP(k1=wp, k2=-1)
    annotation (Placement(transformation(extent={{-60,40},{-40,60}}, rotation=
           0)));
  Modelica.Blocks.Math.Add addD(k1=wd, k2=-1) if
                                        with_D
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}},
          rotation=0)));
  Modelica.Blocks.Math.Add3 addI(k2=-1) if with_I
    annotation (Evaluate=true, Placement(
        transformation(extent={{-60,-60},{-40,-40}}, rotation=0)));
  Modelica.Blocks.Math.Gain P(k=1)
    annotation (Placement(transformation(extent={{-30,40},{-10,
            60}},     rotation=0)));
  Derivative D(k=if pidRepresentation==Types.PID_representation.timeConstants then Td else kd/kp, T=max([Td/Nd,1.e-14]),
    blockType=blockType,
    initType=if init==Types.Init.SteadyState or
                init==Types.Init.InitialOutput then InitWithGlobalDefault.SteadyState else
             if init==Types.Init.InitialState then InitWithGlobalDefault.InitialState else
                InitWithGlobalDefault.NoInit,
    x_start=xd_start,
    y_start=y_start,
    methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.StepExact) if
                          with_D
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  Integrator I(k=if pidRepresentation==Types.PID_representation.timeConstants then 1/Ti else kp/ki, y_start=xi_start,
    blockType=blockType,
    initType=if init==Types.Init.SteadyState then
                InitWithGlobalDefault.SteadyState else
             if init==Types.Init.InitialState then
                InitWithGlobalDefault.InitialState else InitWithGlobalDefault.NoInit,
    methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.StepExact) if with_I
    annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

  Modelica.Blocks.Sources.Constant Izero(k=0) if not with_I
    annotation (Placement(transformation(extent={{20,-54.5},{10,-44.5}},
          rotation=0)));
protected
  parameter Boolean with_I = controllerType==SimpleController.PI or
                             controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
  parameter Boolean with_D = controllerType==SimpleController.PD or
                             controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
public
  Modelica.Blocks.Sources.Constant Dzero(k=0) if not with_D
    annotation (Placement(transformation(extent={{-20,19.5},{-10,29.5}},
          rotation=0)));
  Modelica.Blocks.Math.Add3 addPID
    annotation (Evaluate=true, Placement(transformation(
          extent={{10,-10},{30,10}},rotation=0)));
  Modelica.Blocks.Math.Gain gainPID(k=if pidRepresentation==Types.PID_representation.timeConstants then k else kp)
    annotation (Placement(transformation(extent={{40,-10},
            {60,10}}, rotation=0)));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1) if with_I
    annotation (Evaluate=true, Placement(transformation(
        origin={80,-40},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Math.Gain gainTrack(k=1/(k*Ni)) if with_I
    annotation (Placement(transformation(extent={{40,-80},{20,-60}}, rotation=0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(
    uMax=yMax,
    uMin=yMin,
    limitsAtInit=limitsAtInit)
    annotation (Placement(transformation(extent={{70,-10},{90,10}}, rotation=0)));
  UnitDelay unitDelay(blockType=blockType) if with_I
    annotation (Placement(transformation(extent={{70,-80},{50,-60}})));

initial equation
  if init ==Modelica_LinearSystems2.Controllers.Types.Init.InitialOutput and controllerType<>SimpleController.P then
    y = y_start;
  end if;

equation
  connect(addP.u1, sampler_s.y) annotation (Line(
      points={{-62,56},{-75,56},{-75,0},{-79,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addI.u1, sampler_s.y) annotation (Line(
      points={{-62,-42},{-75,-42},{-75,0},{-79,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addD.u1, sampler_s.y) annotation (Line(
      points={{-62,6},{-75,6},{-75,0},{-79,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler_m.y, addI.u2) annotation (Line(
      points={{6.73556e-016,-79},{6.73556e-016,-75},{-70,-75},{-70,-50},{-62,
          -50}},
      color={0,0,127},
      smooth=Smooth.None,
      thickness=0.5));
  connect(sampler_m.y, addD.u2) annotation (Line(
      points={{6.73556e-016,-79},{0,-79},{0,-75},{-70,-75},{-70,-6},{-62,-6}},
      color={0,0,127},
      thickness=0.5,
      smooth=Smooth.None));
  connect(sampler_m.y, addP.u2) annotation (Line(
      points={{6.73556e-016,-79},{6.73556e-016,-75},{-70,-75},{-70,44},{-62,44}},
      color={0,0,127},
      thickness=0.5,
      smooth=Smooth.None));

  connect(sampler_s.u, u_s) annotation (Line(
      points={{-102,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler_m.u, u_m) annotation (Line(
      points={{-7.34788e-016,-102},{-7.34788e-016,-98.5},{0,-98.5},{0,-120}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P.u, addP.y) annotation (Line(
      points={{-32,50},{-39,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(D.u, addD.y) annotation (Line(
      points={{-32,0},{-39,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(I.u, addI.y) annotation (Line(
      points={{-32,-50},{-39,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P.y, addPID.u1) annotation (Line(
      points={{-9,50},{0,50},{0,8},{8,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(D.y, addPID.u2) annotation (Line(
      points={{-9,0},{8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Dzero.y, addPID.u2) annotation (Line(
      points={{-9.5,24.5},{-5,24.5},{-5,0},{8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(I.y, addPID.u3) annotation (Line(
      points={{-9,-50},{-1,-50},{-1,-8},{8,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Izero.y, addPID.u3) annotation (Line(
      points={{9.5,-49.5},{-0.5,-49.5},{-0.5,-8},{8,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gainPID.y, addSat.u2) annotation (Line(points={{61,0},{65,0},{65,-20},
          {74,-20},{74,-28}},      color={0,0,127}));
  connect(gainPID.y, limiter.u)
    annotation (Line(points={{61,0},{72,0},{68,0}},
                                             color={0,0,127}));
  connect(gainPID.u, addPID.y) annotation (Line(
      points={{38,0},{31,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter.y, y) annotation (Line(
      points={{91,0},{96,0},{96,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u1, limiter.y) annotation (Line(
      points={{86,-28},{86,-20},{95,-20},{95,0},{91,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.u, addSat.y) annotation (Line(
      points={{72,-70},{80,-70},{80,-51}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gainTrack.y, addI.u3) annotation (Line(
      points={{19,-70},{-67,-70},{-67,-58},{-62,-58}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.y, gainTrack.u) annotation (Line(
      points={{49,-70},{42,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    defaultComponentName="PID",
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
          extent={{-35,-22},{79,-60}},
          lineColor={192,192,192},
          textString="PID"),
        Text(
          extent={{-85,93},{91,61}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Documentation(info="<html>
<p>
Via parameter <b>controllerType</b> either <b>P</b>, <b>PI</b>, <b>PD</b>,
or <b>PID</b> can be selected. If, e.g., PI is selected, all components belonging to the
D-part are removed from the block (via conditional declarations).
The example model
<a href=\"modelica://Modelica.Blocks.Examples.PID_Controller\">Modelica.Blocks.Examples.PID_Controller</a>
demonstrates the usage of this controller.
Several practical aspects of PID controller design are incorporated
according to chapter 3 of the book:
</p>

<dl>
<dt>&Aring;str&ouml;m K.J., and H&auml;gglund T. (1995):</dt>
<dd> <b>PID Controllers: Theory, Design, and Tuning</b>.
     Instrument Society of America, 2nd edition, ISBN: 1556175167.
     </dd>
</dl>

<p>
Besides the additive <b>proportional, integral</b> and <b>derivative</b>
part of this controller, the following features are present:
</p>
<ul>
<li> The output of this controller is limited. If the controller is
     in its limits, anti-windup compensation is activated to drive
     the integrator state to zero. </li>
<li> The high-frequency gain of the derivative part is limited
     to avoid excessive amplification of measurement noise.</li>
<li> Setpoint weighting is present, which allows to weight
     the setpoint in the proportional and the derivative part
     independently from the measurement. The controller will respond
     to load disturbances and measurement noise independently of this setting
     (parameters wp, wd). However, setpoint changes will depend on this
     setting. For example, it is useful to set the setpoint weight wd
     for the derivative part to zero, if steps may occur in the
     setpoint signal.</li>
</ul>

<p>
The parameters of the controller can be manually adjusted by performing
simulations of the closed loop system (= controller + plant connected
together) and using the following strategy:
</p>

<ol>
<li> Set very large limits, e.g., yMax = Modelica.Constants.inf</li>
<li> Select a <b>P</b>-controller and manually enlarge parameter <b>k</b>
     (the total gain of the controller) until the closed-loop response
     cannot be improved any more.</li>
<li> Select a <b>PI</b>-controller and manually adjust parameters
     <b>k</b> and <b>Ti</b> (the time constant of the integrator).
     The first value of Ti can be selected, such that it is in the
     order of the time constant of the oscillations occurring with
     the P-controller. If, e.g., vibrations in the order of T=10 ms
     occur in the previous step, start with Ti=0.01 s.</li>
<li> If you want to make the reaction of the control loop faster
     (but probably less robust against disturbances and measurement noise)
     select a <b>PID</b>-Controller and manually adjust parameters
     <b>k</b>, <b>Ti</b>, <b>Td</b> (time constant of derivative block).</li>
<li> Set the limits yMax and yMin according to your specification.</li>
<li> Perform simulations such that the output of the PID controller
     goes in its limits. Tune <b>Ni</b> (Ni*Ti is the time constant of
     the anti-windup compensation) such that the input to the limiter
     block (= limiter.u) goes quickly enough back to its limits.
     If Ni is decreased, this happens faster. If Ni=infinity, the
     anti-windup compensation is switched off and the controller works bad.</li>
</ol>

<h4>Initialization</h4>

<p>
This block can be initialized in different
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

<blockquote><pre>
<b>der</b>(y) = k*u;
</pre></blockquote>

<p>
The steady state equation &quot;der(x)=0&quot; leads to the condition that the input u to the
integrator is zero. If the input u is already (directly or indirectly) defined
by another initial condition, then the initialization problem is <b>singular</b>
(has none or infinitely many solutions). This situation occurs often
for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and
since speed is both a state and a derivative, it is natural to
initialize it with zero. As sketched this is, however, not possible.
The solution is to not initialize u_m or the variable that is used
to compute u_m by an algebraic equation.
</p>

<p>
If parameter <b>limitAtInit</b> = <b>false</b>, the limits at the
output of this controller block are removed from the initialization problem which
leads to a much simpler equation system. After initialization has been
performed, it is checked via an assert whether the output is in the
defined limits. For backward compatibility reasons
<b>limitAtInit</b> = <b>true</b>. In most cases it is best
to use <b>limitAtInit</b> = <b>false</b>.
</p>
</html>"));
end LimPID;

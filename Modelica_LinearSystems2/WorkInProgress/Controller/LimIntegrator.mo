within Modelica_LinearSystems2.WorkInProgress.Controller;
block LimIntegrator
  "Output is the integral limited to upper or lower limit (continiuous or discrete)"
  extends Interfaces.PartialSampledBlock;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault;
  import Modelica_LinearSystems2.Controller.Interfaces;

  parameter Real k=1 "Integrator gain";
  parameter Boolean withDelay=false
    "True, if the output is delayed by one sample period (only if discrete)";

  parameter Real y_start=0 "Initial or guess value of output (=state)" annotation(Dialog(tab="Advanced options"));

  parameter Boolean limitsAtInit=true
    "= false, if limits are ignored during initialization (i.e., y=u)";

public
  Modelica.Blocks.Interfaces.RealInput limit1
    "Connector of Real input signal used as maximum of input u"
    annotation (Placement(transformation(extent={{-140,60},{-100,100}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput limit2
    "Connector of Real input signal used as minimum of input u"
    annotation (Placement(transformation(extent={{-140,-100},{-100,-60}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signals of block"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signals of block"
    annotation (Placement(transformation(extent={{100,-10},{120,10}}, rotation=0)));
  Modelica_LinearSystems2.WorkInProgress.Controller.Integrator integrator(
    k=k,
    blockType=blockType,
    methodType=methodType,
    initType=initType,
    y_start=y_start,
    limitsAtInit=limitsAtInit,
    sampleFactor=sampleFactor,
    withDelay=false)
    annotation (Placement(transformation(extent={{10,-10},{30,10}})));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1)
    annotation (Evaluate=true, Placement(transformation(
        origin={32,-28},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_LinearSystems2.Controller.Internal.Add add
    annotation (Placement(transformation(extent={{-20,10},{0,-10}})));
  Modelica_LinearSystems2.Controller.Sampler sampler1(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-84,-10},{-64,10}})));
  Modelica_LinearSystems2.Controller.Sampler sampler2(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-80,-78},{-60,-58}})));
  Modelica_LinearSystems2.Controller.Sampler sampler3(sampleFactor=sampleFactor,
      blockType=blockType) annotation (Placement(
        transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-30,-50})));
  Modelica.Blocks.Math.Gain gain(k=1/k)
    annotation (Placement(transformation(extent={{12,-38},{-8,-18}})));

  Modelica_LinearSystems2.Controller.UnitDelay unitDelay(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-20,-38},{-40,-18}})));
  Modelica.Blocks.Math.Abs abs1
    annotation (Placement(transformation(extent={{-52,-50},{-72,-30}})));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold=0)
    annotation (Placement(transformation(extent={{-82,-50},{-102,-30}})));
  Modelica.Blocks.Logical.Switch switch1
    annotation (Placement(transformation(extent={{-50,10},{-30,30}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-84,24},{-64,44}})));
equation
  connect(variableLimiter.y, y) annotation (Line(
      points={{81,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableLimiter.u, integrator.y) annotation (Line(
      points={{58,0},{31,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u2, integrator.y) annotation (Line(
      points={{44,-22},{48,-22},{48,0},{31,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u1, variableLimiter.y) annotation (Line(
      points={{44,-34},{90,-34},{90,0},{81,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, integrator.u) annotation (Line(
      points={{-1,0},{8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler1.u, u) annotation (Line(
      points={{-86,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler2.u, limit2) annotation (Line(
      points={{-82,-68},{-90,-68},{-90,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, addSat.y) annotation (Line(
      points={{14,-28},{21,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, sampler3.u) annotation (Line(
      points={{-9,-28},{-12,-28},{-12,-50},{-18,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator.limit2, limit2) annotation (Line(
      points={{14,12},{14,56},{-100,56},{-100,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.y, add.u1) annotation (Line(
      points={{-41,-28},{-44,-28},{-44,-12},{-10,-12},{-10,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.u, gain.y) annotation (Line(
      points={{-18,-28},{-9,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit1, integrator.limit1) annotation (Line(
      points={{-120,80},{26,80},{26,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit1, variableLimiter.limit1) annotation (Line(
      points={{-120,80},{40,80},{40,8},{58,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit2, variableLimiter.limit2) annotation (Line(
      points={{-120,-80},{52,-80},{52,-8},{58,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(abs1.u, addSat.y) annotation (Line(
      points={{-50,-40},{-44,-40},{-44,-72},{20,-72},{20,-28},{21,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(greaterThreshold.u, abs1.y) annotation (Line(
      points={{-80,-40},{-73,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(greaterThreshold.y, switch1.u2) annotation (Line(
      points={{-103,-40},{-130,-40},{-130,20},{-52,20}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(switch1.u1, const.y) annotation (Line(
      points={{-52,28},{-56,28},{-56,34},{-63,34}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch1.u3, sampler1.y) annotation (Line(
      points={{-52,12},{-56,12},{-56,0},{-63,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch1.y, add.u2) annotation (Line(
      points={{-29,20},{-22,20},{-22,0},{-18,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <em>integrator</em>:
</p>
<pre>
          k
     y = --- * u
          s
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
<p>
It is not possible to initialize a continuous integrator in steady state.
For this reason, option \"initType = SteadyState\" is ignored for
a continuous integrator block and
interpreted as \"initType = InitialState\".
</p>
</html>"), Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
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
        Text(
          extent={{0,-10},{60,-70}},
          textColor={192,192,192},
          textString="I"),
        Text(
          extent={{-150,-150},{150,-110}},
          textColor={0,0,0},
          textString="k=%k"),
        Line(points={{-80,-80},{80,80}}, color={0,0,127}),
        Text(
          extent={{-94,78},{88,46}},
          textColor={0,0,0},
          textString="%sampleFactor")}));
end LimIntegrator;

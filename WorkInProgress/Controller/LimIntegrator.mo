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

  parameter Real y_start=0 "Initial or guess value of output (=state)"                                                               annotation(Dialog(tab="Advanced options"));

  parameter Boolean limitsAtInit=true
    "= false, if limits are ignored during initializiation (i.e., y=u)";

public
  Modelica.Blocks.Interfaces.RealInput limit1
    "Connector of Real input signal used as maximum of input u"
                              annotation (Placement(transformation(extent={{-140,60},
            {-100,100}},          rotation=0)));
  Modelica.Blocks.Interfaces.RealInput limit2
    "Connector of Real input signal used as minimum of input u"
                              annotation (Placement(transformation(extent={{-140,
            -100},{-100,-60}},      rotation=0)));

  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signals of block"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signals of block"
    annotation (Placement(transformation(extent={{100,-10},{120,10}}, rotation=
            0)));
  Modelica_LinearSystems2.WorkInProgress.Controller.Integrator
                                                integrator(
    k=k,
    blockType=blockType,
    methodType=methodType,
    initType=initType,
    y_start=y_start,
    limitsAtInit=limitsAtInit,
    sampleFactor=sampleFactor,
    withDelay=false)
    annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter(limitsAtInit=
        limitsAtInit)
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1)
    annotation (Evaluate=true, Placement(transformation(
        origin={32,-28},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_LinearSystems2.Controller.Internal.Add add
    annotation (Placement(transformation(extent={{-50,8},{-30,-12}})));
  Modelica_LinearSystems2.Controller.Sampler sampler1(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-94,-10},{-74,10}})));
  Modelica_LinearSystems2.Controller.Sampler sampler2(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
  Modelica_LinearSystems2.Controller.Sampler sampler3(sampleFactor=sampleFactor,
      blockType=blockType)                            annotation (Placement(
        transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-36,-64})));
  Modelica.Blocks.Math.Gain gain(k=1/k)
    annotation (Placement(transformation(extent={{12,-38},{-8,-18}})));

  Modelica_LinearSystems2.Controller.UnitDelay unitDelay(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-20,-38},{-40,-18}})));
  Modelica.Blocks.Math.Abs abs1
    annotation (Placement(transformation(extent={{-52,-46},{-72,-26}})));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold=0)
    annotation (Placement(transformation(extent={{-86,-46},{-106,-26}})));
  Modelica.Blocks.Logical.Switch switch1
    annotation (Placement(transformation(extent={{-56,14},{-36,34}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-132,26},{-112,46}})));
equation
  connect(variableLimiter.y, y) annotation (Line(
      points={{81,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableLimiter.u, integrator.y) annotation (Line(
      points={{58,0},{26,0},{26,-2},{11,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u2, integrator.y) annotation (Line(
      points={{44,-22},{40,-22},{40,-2},{11,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u1, variableLimiter.y) annotation (Line(
      points={{44,-34},{90,-34},{90,0},{81,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, integrator.u) annotation (Line(
      points={{-31,-2},{-12,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler1.u, u) annotation (Line(
      points={{-96,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler2.u, limit2) annotation (Line(
      points={{-82,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, addSat.y) annotation (Line(
      points={{14,-28},{21,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, sampler3.u) annotation (Line(
      points={{-9,-28},{-24,-28},{-24,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator.limit2, limit2) annotation (Line(
      points={{-6,10},{-6,40},{-100,40},{-100,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.y, add.u1) annotation (Line(
      points={{-41,-28},{-40,-28},{-40,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.u, gain.y) annotation (Line(
      points={{-18,-28},{-9,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit1, integrator.limit1) annotation (Line(
      points={{-120,80},{-64,80},{-64,82},{6,82},{6,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit1, variableLimiter.limit1) annotation (Line(
      points={{-120,80},{40,80},{40,8},{58,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limit2, variableLimiter.limit2) annotation (Line(
      points={{-120,-80},{-44,-80},{-44,-82},{50,-82},{50,-8},{58,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(abs1.u, addSat.y) annotation (Line(
      points={{-50,-36},{-21,-36},{-21,-28},{21,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(greaterThreshold.u, abs1.y) annotation (Line(
      points={{-84,-36},{-73,-36}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(greaterThreshold.y, switch1.u2) annotation (Line(
      points={{-107,-36},{-130,-36},{-130,24},{-58,24}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(switch1.u1, const.y) annotation (Line(
      points={{-58,32},{-92,32},{-92,36},{-111,36}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch1.u3, sampler1.y) annotation (Line(
      points={{-58,16},{-72,16},{-72,0},{-73,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch1.y, add.u2) annotation (Line(
      points={{-35,24},{-32,24},{-32,-2},{-48,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Window(
      x=0.29,
      y=0.05,
      width=0.53,
      height=0.54),
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>integrator</i>:
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
It is not possible to initalize a continuous integrator in steady state.
For this reason, option \"initType = SteadyState\" is ignored for
a continuous integrator block and
interpreted as \"initType = InitialState\".
</p>
</html>
"), Icon(coordinateSystem(
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
          lineColor={192,192,192},
          textString="I"),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="k=%k"),
        Line(points={{-80,-80},{80,80}}, color={0,0,127}),
        Text(
          extent={{-94,78},{88,46}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics));
end LimIntegrator;

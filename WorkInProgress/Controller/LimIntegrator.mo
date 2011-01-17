within Modelica_LinearSystems2.WorkInProgress.Controller;
block LimIntegrator
  "Output is the integral limited to upper or lower limit (continiuous or discrete)"
  import Modelica_LinearSystems2;

  import Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault;
  import Modelica_LinearSystems2.Controller.Interfaces;

  extends Interfaces.PartialSampledBlock;

  extends Interfaces.PartialBlockIcon;

  parameter Real k=1 "Integrator gain";
  parameter Boolean withDelay=false
    "= true, if the output is delayed by one sample period (only if discrete)";

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
    withDelay=withDelay,
    blockType=blockType,
    methodType=methodType,
    initType=initType,
    y_start=y_start,
    limitsAtInit=limitsAtInit,
    sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter(limitsAtInit=
        limitsAtInit)
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1)
    annotation (Evaluate=true, Placement(transformation(
        origin={20,-30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_LinearSystems2.Controller.Internal.Add add
    annotation (Placement(transformation(extent={{-70,10},{-50,-10}})));
  Modelica_LinearSystems2.Controller.Sampler sampler(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-100,70},{-80,90}})));
  Modelica_LinearSystems2.Controller.Sampler sampler1(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  Modelica_LinearSystems2.Controller.Sampler sampler2(blockType=blockType,
      sampleFactor=sampleFactor)
    annotation (Placement(transformation(extent={{-100,-90},{-80,-70}})));
  Modelica_LinearSystems2.Controller.Sampler sampler3(sampleFactor=sampleFactor,
      blockType=blockType)                            annotation (Placement(
        transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-40,-30})));
  Modelica.Blocks.Math.Gain gain(k=1/k)
    annotation (Placement(transformation(extent={{0,-40},{-20,-20}})));

equation
  connect(variableLimiter.y, y) annotation (Line(
      points={{81,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableLimiter.u, integrator.y) annotation (Line(
      points={{58,0},{-9,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u2, integrator.y) annotation (Line(
      points={{32,-24},{40,-24},{40,0},{-9,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.u1, variableLimiter.y) annotation (Line(
      points={{32,-36},{90,-36},{90,0},{81,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, integrator.u) annotation (Line(
      points={{-51,0},{-32,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler1.y, add.u2) annotation (Line(
      points={{-79,0},{-68,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler1.u, u) annotation (Line(
      points={{-102,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.y, variableLimiter.limit1) annotation (Line(
      points={{-79,80},{50,80},{50,8},{58,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.u, limit1) annotation (Line(
      points={{-102,80},{-120,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler2.y, variableLimiter.limit2) annotation (Line(
      points={{-79,-80},{50,-80},{50,-8},{58,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler2.u, limit2) annotation (Line(
      points={{-102,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, addSat.y) annotation (Line(
      points={{2,-30},{9,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, sampler3.u) annotation (Line(
      points={{-21,-30},{-28,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator.limit2, limit2) annotation (Line(
      points={{-26,12},{-26,40},{-100,40},{-100,-80},{-120,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator.limit1, sampler.y) annotation (Line(
      points={{-14,12},{-14,80},{-79,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler3.y, add.u1) annotation (Line(
      points={{-51,-30},{-60,-30},{-60,-8}},
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

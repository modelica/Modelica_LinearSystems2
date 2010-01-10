within Modelica_LinearSystems2.Controller.Internal;
block DiscreteInterpolator
  "Increasing the sampling frequency with linear interpolation"
  extends Interfaces.PartialBlockIcon;
  parameter Integer outputSampleFactor(min=1)=1
    "Output sample time = outputSampleFactor * sampleClock.sampleTime" 
     annotation (Dialog(enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));
  parameter Integer inputSampleFactor(min=1)=1
    "Input sample time = inputSampleFactor * outputSampleFactor * sampleClock.sampleTime"
     annotation (Dialog(enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block" 
    annotation (extent=[-140, -20; -100, 20]);
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block" 
    annotation (extent=[100, -10; 120, 10]);
protected
  outer SampleClock sampleClock "Global options";
  Boolean outputSampleTrigger "= true, if output sample time";
  Integer outputTicks(start=0, fixed=true);
  Integer ticks(start=1, fixed=true);
  Real pre_u;
  Boolean sampleIn;
equation
  when sampleClock.sampleTrigger then
      outputTicks = if pre(outputTicks) < outputSampleFactor then pre(outputTicks) + 1 else 1;
  end when;
  outputSampleTrigger = sampleClock.sampleTrigger and outputTicks >= outputSampleFactor;

  when outputSampleTrigger then
     ticks = if pre(ticks) < inputSampleFactor then pre(ticks) + 1 else 1;
     sampleIn = ticks == 1;
     y = if ticks == 1 then (u - pre(u))/inputSampleFactor + pre(u) else 
                             pre(pre_u) + ticks/inputSampleFactor*(pre(u) - pre(pre_u));
  end when;

  when {initial(), sampleIn} then
      pre_u = pre(u);
  end when;

initial equation
  u = pre(u);
  pre_u = pre(pre_u);
  u = y;

  annotation (
    Window(
      x=0.27,
      y=0.1,
      width=0.57,
      height=0.75),
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-30,78},{-30,-90}}, color={192,192,192}),
        Polygon(
          points={{-30,92},{-38,70},{-22,70},{-30,90},{-30,92}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-42,-78},{82,-78}}, color={192,192,192}),
        Polygon(
          points={{90,-78},{68,-70},{68,-86},{90,-78}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-72,34},{-58,-68}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-46,34},{-84,34},{-66,68},{-46,34}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-30,-12},{20,40},{88,52}},
          color={0,0,0},
          smooth=Smooth.None),
        Rectangle(extent={{-36,-6},{-24,-18}}, lineColor={0,0,0}),
        Rectangle(extent={{14,46},{26,34}}, lineColor={0,0,0}),
        Rectangle(extent={{82,58},{94,46}}, lineColor={0,0,0}),
        Ellipse(
          extent={{-20,10},{-8,0}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-4,26},{8,16}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{34,48},{46,38}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{58,52},{70,42}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
<p>
</p>
</html>
"));

end DiscreteInterpolator;

within Modelica_LinearSystems2.Controller;
block SampleClock
  "Global options for blocks of Controller library (in particular sample clock)"

  parameter Types.BlockType blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Continuous
    "Type of Sampled blocks"
    annotation (
      Evaluate=true,
      choices(__Dymola_radioButtons=true, choice=Modelica_LinearSystems2.Controller.Types.BlockType.Continuous
        "Continuous",
        choice=Modelica_LinearSystems2.Controller.Types.BlockType.Discrete
        "Discrete"));
  parameter Modelica_LinearSystems2.Types.Method methodType=
      Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method for discrete blocks";
  parameter Modelica.SIunits.Time sampleTime = 1
    "Base sample time for discrete blocks";
  parameter Types.Init initType=Modelica_LinearSystems2.Controller.Types.Init.SteadyState
    "Type of initialization of Sampled blocks" annotation(Evaluate=true);
  output Boolean sampleTrigger "Trigger that is true at every sampleTime"
                                                                         annotation(HideResult=true);
equation
  if blockType == Types.BlockType.Continuous then
     sampleTrigger = false;
  else
     sampleTrigger = sample(sampleTime, sampleTime);
  end if;
  annotation (
    defaultComponentName="sampleClock",
    defaultComponentPrefixes="inner",
    missingInnerMessage="A \"sampleClock\" component is not defined. A default
sampleClock with blockType = Continuous will be used. If this is not desired,
drag Modelica_LinearSystems2.Controller.SampleClock into the top level of your model.",
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics={
        Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Sphere,
          fillColor={255,255,255}),
        Line(points={{-100,0},{-36,0}}, color={0,0,0}),
        Line(points={{36,0},{100,0}}, color={0,0,0}),
        Line(points={{-35,0},{30,35}}, color={0,0,0}),
        Text(
          extent={{-150,-110},{150,-140}},
          lineColor={0,0,0},
          textString="%sampleTime s"),
        Ellipse(
          extent={{-25,-10},{-45,10}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,-20},{150,-50}},
          lineColor={0,0,0},
          textString="Continuous",
          visible=blockType==Modelica_LinearSystems2.Controller.Types.BlockType.Continuous),
        Text(
          extent={{-150,-20},{150,-50}},
          lineColor={0,0,0},
          textString="Discrete",
          visible=blockType==Modelica_LinearSystems2.Controller.Types.BlockType.Discrete),
        Text(
          extent={{-150,140},{150,100}},
          lineColor={0,0,255},
          fillColor={169,199,255},
          fillPattern=FillPattern.Solid,
          textString="%name"),
        Ellipse(
          extent={{45,-10},{25,10}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<html>
<p>
Global block that defines options for all components of the
Controller library that are on the same or on a lower level
as the sampleClock component. In particular it is defined whether
the blocks shall be used by default in a continuous or a
discrete representation. In the latter case, the default
discretization method and the base sample time is defined.
The sample time of a block is an integer multiple of the base sample
time defined in the SampleClock component.
</p>
</html>"));
end SampleClock;

within Modelica_LinearSystems2.Controller.Templates;
partial model TwoDegreesOfFreedomControlSISO
  "Template of a controller with two structural degrees of freedom and an inverse plant model in forward path"
  parameter Integer l = 1 "number of measurable outputs";
  parameter Boolean additionalMeasurableOutputs = true;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-140,
            -100},{140,100}}), graphics={Text(
          extent={{4,10},{36,4}},
          lineColor={0,0,255},
          textString="controller")}),     Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-140,-100},{140,100}})),
    Documentation(info="<html>
<p>
Template of a controller with two structural degrees of freedom and an inverse plant model in forward path.
The functionality of such contorl system structures is discribed in [1]
<p>

<A name=\"References\"><B><FONT SIZE=\"+1\">References</FONT></B></A> <PRE>
  [1] Looye, G. et al, \"Nonlinear inverse moldes for control\",
      Proceedings Modelica Conference 2005, pp. 267-279, 2005.
</PRE>



</html>"));
  Modelica.Blocks.Math.Feedback feedback[plant.l]
    annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
  Modelica_LinearSystems2.Controller.Internal.Add2 add
    annotation (Placement(transformation(extent={{50,-20},{70,0}})));
  Modelica_LinearSystems2.Controller.Filter filter
    annotation (Placement(transformation(extent={{-100,10},{-80,30}})));
  replaceable Modelica_LinearSystems2.Controller.Interfaces.PartialSISO
    controller constrainedby Interfaces.PartialSISO
    annotation (Placement(transformation(extent={{10,-20},{30,0}})));
  replaceable Plant_SISO plant(l=l,
      additionalMeasurableOutputs=additionalMeasurableOutputs) constrainedby
    PlantTemplate_SISO
    annotation (Placement(transformation(extent={{90,-20},{110,0}})));
  Modelica.Blocks.Math.InverseBlockConstraints forwardControlModel
    annotation (Placement(transformation(extent={{-68,6},{-14,34}})));
  replaceable Plant_SISO plant2(l=l, additionalMeasurableOutputs=
        additionalMeasurableOutputs) constrainedby PlantTemplate_SISO
    annotation (Placement(transformation(extent={{-30,9},{-50,30}})));
equation
  connect(controller.u, feedback[1].y)  annotation (Line(
      points={{8,-10},{-1,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, feedback.u2) annotation (Line(
      points={{100,-21},{100,-40},{-10,-40},{-10,-18}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant2.y, forwardControlModel.u2)         annotation (Line(
      points={{-51,19.5},{-62,19.5},{-62,20},{-62.6,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant2.u, forwardControlModel.y2)         annotation (Line(
      points={{-28,19.5},{-26,19.5},{-26,20},{-18.05,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, plant2.ym)     annotation (Line(
      points={{-18,-10},{-40,-10},{-40,7.95}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter.y, forwardControlModel.u1)     annotation (Line(
      points={{-79,20},{-70.7,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(controller.y, add.u2[1]) annotation (Line(
      points={{31,-10},{52,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1[1], forwardControlModel.y1)
                                      annotation (Line(
      points={{60,-2},{60,20},{-12.65,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y[1], plant.u) annotation (Line(
      points={{69,-10},{88,-10}},
      color={0,0,127},
      smooth=Smooth.None));
end TwoDegreesOfFreedomControlSISO;

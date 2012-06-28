within Modelica_LinearSystems2.Controller.Templates;
partial model TwoDOFinverseModelController
  "Template for a controller with two structural degrees of freedom and an inverse plant model in forward path"
  parameter Boolean additionalMeasurableOutputs = true
    "Enable additional output vector of dimension l";
  parameter Integer l=1 "Number of measurable outputs"
    annotation (Dialog(enable=additionalMeasurableOutputs));

  Modelica.Blocks.Math.Feedback feedback[plant.l]
    annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
  Controller.Internal.Add2 add
    annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
  Controller.Filter filter
    annotation (Placement(transformation(extent={{-90,10},{-70,30}})));
  replaceable Controller.Interfaces.PartialSISO controller
    constrainedby Controller.Interfaces.PartialSISO
    annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
  replaceable Controller.Templates.Internal.Plant_SISO plant(l=l,
      additionalMeasurableOutputs=additionalMeasurableOutputs) constrainedby
    Templates.Internal.PlantTemplate_SISO
    annotation (Placement(transformation(extent={{70,-30},{90,-10}})));
  Modelica.Blocks.Math.InverseBlockConstraints inverseModel
    annotation (Placement(transformation(extent={{-58,6},{-4,34}})));
  replaceable Internal.Plant_SISO plant_inv(l=l,
      additionalMeasurableOutputs=additionalMeasurableOutputs) constrainedby
    Templates.Internal.PlantTemplate_SISO
    annotation (Placement(transformation(extent={{-20,10},{-40,30}})));
equation
  connect(controller.u, feedback[1].y)  annotation (Line(
      points={{8,-20},{-1,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, feedback.u2) annotation (Line(
      points={{80,-31},{80,-50},{-10,-50},{-10,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant_inv.y, inverseModel.u2)             annotation (Line(
      points={{-41,20},{-52.6,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant_inv.u, inverseModel.y2)             annotation (Line(
      points={{-18,20},{-8.05,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, plant_inv.ym)  annotation (Line(
      points={{-18,-20},{-30,-20},{-30,9}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter.y, inverseModel.u1)            annotation (Line(
      points={{-69,20},{-60.7,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(controller.y, add.u2[1]) annotation (Line(
      points={{31,-20},{42,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1[1], inverseModel.y1) annotation (Line(
      points={{50,-12},{50,20},{-2.65,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y[1], plant.u) annotation (Line(
      points={{59,-20},{68,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}),       graphics={Text(
          extent={{-58,42},{-2,36}},
          lineColor={0,0,0},
          textString="inverted system")}),Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<p>
Template of a controller with two structural degrees of freedom and an inverse plant model in forward path.
The functionality of such contorl system structures is discribed in [1]
</p>

<h4><a name=\"References\">R</a>eferences</h4>
<pre>
  [1] Looye, G. et al, \"Nonlinear inverse moldes for control\",
      Proceedings Modelica Conference 2005, pp. 267-279, 2005.
</pre>
</html>"));
end TwoDOFinverseModelController;

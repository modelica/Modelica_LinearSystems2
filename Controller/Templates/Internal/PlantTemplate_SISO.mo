within Modelica_LinearSystems2.Controller.Templates.Internal;
partial model PlantTemplate_SISO "SISO plant template"
  parameter Boolean additionalMeasurableOutputs = true
    "Enable additional output vector of dimension l";
  parameter Integer l=1 "Number of measurable outputs"
    annotation (Dialog(enable=additionalMeasurableOutputs));

  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput ym[l] if additionalMeasurableOutputs
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-110})));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics));
end PlantTemplate_SISO;

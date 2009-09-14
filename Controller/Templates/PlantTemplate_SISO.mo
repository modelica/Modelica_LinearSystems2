within Modelica_LinearSystems2.Controller.Templates;
partial model PlantTemplate_SISO "SISO plant template"
  parameter Integer l = 1 "number of measurable outputs";
  parameter Boolean additionalMeasurableOutputs = true;

  Modelica.Blocks.Interfaces.RealOutput y 
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput u 
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics));
  Modelica.Blocks.Interfaces.RealOutput ym[l] if additionalMeasurableOutputs 
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-110})));
end PlantTemplate_SISO;

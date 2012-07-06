within Modelica_LinearSystems2.Controller.Templates;
partial model PlantTemplate "MIMO plant template"
  parameter Integer n=1 "system order";
  parameter Integer m=1 "number of inputs";
  parameter Boolean additionalMeasurableOutputs = true
    "Enable additional output vector of dimension l";
  parameter Integer l=1 "Number of measurable outputs"
    annotation (Dialog(enable=additionalMeasurableOutputs));

  Modelica.Blocks.Interfaces.RealOutput y[n]
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput u[m]
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput ym[l] if additionalMeasurableOutputs
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-110})));
end PlantTemplate;

within Modelica_LinearSystems2.Controller.Templates;
model Plant_SISO "SISO dummy plant model which can be used in templates"
  extends PlantTemplate_SISO;
equation
  connect(u, ym[1]) annotation (Line(
      points={{-120,0},{0,0},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, y) annotation (Line(
      points={{-120,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}), graphics));
end Plant_SISO;

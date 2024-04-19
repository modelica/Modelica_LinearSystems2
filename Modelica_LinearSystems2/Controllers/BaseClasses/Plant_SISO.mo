within Modelica_LinearSystems2.Controllers.BaseClasses;
model Plant_SISO "SISO dummy plant model which can be used in templates"
  extends PartialPlantSISO;
equation
  connect(u, ym[1]) annotation (Line(
      points={{-120,0},{0,0},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, y) annotation (Line(
      points={{-120,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
end Plant_SISO;

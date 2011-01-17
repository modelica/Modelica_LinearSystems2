within Modelica_LinearSystems2.Controller.Examples.Components;
block TwoPoint "Block with two-point behaviour"
  extends Modelica.Blocks.Interfaces.SISO;
  parameter Real b(min=0)=1;

equation
  y = if u<0 then -b else b;

  annotation (Icon(graphics={
        Line(
          points={{-4,70},{-4,-86}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-90,0},{80,0}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-74,-36},{-4,-36},{-4,32},{66,32}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5)}));
end TwoPoint;

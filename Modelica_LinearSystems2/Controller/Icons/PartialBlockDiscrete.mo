within Modelica_LinearSystems2.Controller.Icons;
partial block PartialBlockDiscrete "Icon for a discrete block"

  annotation (
    Icon(
      graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          fillColor={213,255,190},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Text(
          extent={{-150,150},{150,110}},
          textString="%name",
          lineColor={0,0,255}),
        Rectangle(
          extent={{-100,100},{100,-100}},
          fillColor={230,230,255},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None,
          lineColor={0,0,0})}),
    Documentation(info="<html>
<p>This partial class is intended to design a <em>default icon for a discrete block</em>.</p>
</html>"));
end PartialBlockDiscrete;

within Modelica_LinearSystems2.Controller.Icons;
partial block PartialBlockIcon "Icon for a discrete/continuous block"

protected
  Boolean cont=true;

  annotation (
    Icon(
      coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
        graphics={
        Rectangle(
          visible=not cont,
          extent={{-100,100},{100,-100}},
          fillColor={213,255,190},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None),
        Rectangle(
          visible=cont,
          extent={{-100,100},{100,-100}},
          fillColor={230,230,255},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Text(
          extent={{-150,140},{150,100}},
          lineColor={0,0,255},
          fillColor={169,199,255},
          fillPattern=FillPattern.Solid,
          textString="%name")}), Documentation(info="<html>
<p>This partial class is intended to design a <i>default icon for a discrete or continuous block</i>. The background color of this icon depends on the boolean parameter <code>cont</code>.</p>
</html>"));

end PartialBlockIcon;

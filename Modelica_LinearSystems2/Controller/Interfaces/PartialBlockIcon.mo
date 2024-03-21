within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialBlockIcon
  "Obsolete model: Basic graphical layout of discrete/continuous block"
  extends Modelica.Icons.ObsoleteModel;

  annotation (
    obsolete = "Obsolete function - use Modelica_LinearSystems2.Controller.Icons.PartialBlockIcon instead",
    Documentation(info="<html>
<p>This model is obsolete. Use <a href=\"Modelica://Modelica_LinearSystems2.Controller.Icons.PartialBlockIcon\">PartialBlockIcon</a> instead.</p>
</html>"),
    Icon(
      graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          fillColor={230,230,255},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Text(
          extent={{-150,150},{150,110}},
          textString="%name",
          textColor={0,0,255})}));
end PartialBlockIcon;

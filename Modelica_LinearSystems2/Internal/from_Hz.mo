within Modelica_LinearSystems2.Internal;
function from_Hz "Obsolete class: Convert from Hz to rad/s"
  extends Modelica.Units.Icons.Conversion;
  extends Modelica.Icons.ObsoleteModel;

  input Modelica.Units.SI.Frequency f "frequency";
  output Modelica.Units.SI.AngularVelocity w "angular velocity";

algorithm
  w := 2*Modelica.Constants.pi*f;
  annotation (
    obsolete = "Obsolete function - use Modelica.Units.Conversions.from_Hz instead",
    Icon(graphics={
        Text(
          extent={{-20,100},{-100,20}},
          textColor={0,0,0},
          textString="Hz"),
        Text(
          extent={{100,-20},{20,-100}},
          textColor={0,0,0},
          textString="rad/s")}), Documentation(info="<html>
<p>This model is obsolete. Use <a href=\"modelica://Modelica.Units.Conversions.from_Hz\">Modelica.Units.Conversions.from_Hz</a> instead.</p>
</html>"));
end from_Hz;

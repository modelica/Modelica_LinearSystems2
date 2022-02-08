within Modelica_LinearSystems2.Internal;
function to_Hz "Obsolete class: Convert from rad/s to Hz"
  extends Modelica.Units.Icons.Conversion;
  extends Modelica.Icons.ObsoleteModel;

  input Modelica.Units.SI.AngularVelocity w "angular velocity";
  output Modelica.Units.SI.Frequency f "frequency";
algorithm
  f := w/(2*Modelica.Constants.pi);
  annotation (
    obsolete = "Obsolete function - use Modelica.Units.Conversions.to_Hz instead",
    Icon(graphics={
        Text(
          extent={{-20,100},{-100,20}},
          textColor={0,0,0},
          textString="rad/s"),
        Text(
          extent={{100,-20},{20,-100}},
          textColor={0,0,0},
          textString="Hz")}), Documentation(info="<html>
<p>This model is obsolete. Use <a href=\"Modelica://Modelica.Units.Conversions.to_Hz\">Modelica.Units.Conversions.to_Hz</a> instead.</p>
</html>"));
end to_Hz;

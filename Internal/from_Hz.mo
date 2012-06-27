within Modelica_LinearSystems2.Internal;
function from_Hz "Obsolete class: Convert from Hz to rad/s"
  extends Modelica.SIunits.Conversions.ConversionIcon;
  extends Modelica.Icons.ObsoleteModel;

  input Modelica.SIunits.Frequency f "frequency";
  output Modelica.SIunits.AngularVelocity w "angular velocity";

algorithm
  w := 2*Modelica.Constants.pi*f;
  annotation (Icon(graphics={Text(
          extent={{-20,100},{-100,20}},
          lineColor={0,0,0},
          textString="K"), Text(
          extent={{100,-20},{20,-100}},
          lineColor={0,0,0},
          textString="°C")}), Documentation(info="<html>
<p>This model is obsolete. Use <a href=\"Modelica://Modelica.SIunits.Conversions.from_Hz\">Modelica.SIunits.Conversions.from_Hz</a> instead.</p>
</html>"));
end from_Hz;

within Modelica_LinearSystems2.Internal;
function to_Hz "Obsolete class: Convert from rad/s to Hz"
  extends Modelica.SIunits.Conversions.ConversionIcon;
  extends Modelica.Icons.ObsoleteModel;

  input Modelica.SIunits.AngularVelocity w "angular velocity";
  output Modelica.SIunits.Frequency f "frequency";
algorithm
  f := w/(2*Modelica.Constants.pi);
  annotation (Icon(graphics={Text(
          extent={{-20,100},{-100,20}},
          lineColor={0,0,0},
          textString="rad/s"), Text(
          extent={{100,-20},{20,-100}},
          lineColor={0,0,0},
          textString="Hz")}), Documentation(info="<html>
<p>This model is obsolete. Use <a href=\"Modelica://Modelica.SIunits.Conversions.to_Hz\">Modelica.SIunits.Conversions.to_Hz</a> instead.</p>
</html>"));
end to_Hz;

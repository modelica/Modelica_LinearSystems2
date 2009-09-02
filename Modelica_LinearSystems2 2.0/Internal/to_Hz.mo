within Modelica_LinearSystems2.Internal;
function to_Hz "Convert from rad/s to Hz"
  extends Modelica.SIunits.Conversions.ConversionIcon;
  input Modelica.SIunits.AngularVelocity w "angular velocity";
  output Modelica.SIunits.Frequency f "frequency";
  annotation (Icon(graphics={Text(
          extent={{-20,100},{-100,20}},
          lineColor={0,0,0},
          textString="K"), Text(
          extent={{100,-20},{20,-100}},
          lineColor={0,0,0},
          textString="°C")}));
algorithm
  f := w/(2*Modelica.Constants.pi);
end to_Hz;

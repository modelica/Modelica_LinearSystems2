within Modelica_LinearSystems2.Internal;
function from_Hz "Convert from Hz to rad/s"
  extends Modelica.SIunits.Conversions.ConversionIcon;
  input Modelica.SIunits.Frequency f "frequency";
  output Modelica.SIunits.AngularVelocity w "angular velocity";

  annotation (Icon(graphics={Text(
          extent={{-20,100},{-100,20}},
          lineColor={0,0,0},
          textString="K"), Text(
          extent={{100,-20},{20,-100}},
          lineColor={0,0,0},
          textString="°C")}));
algorithm
  w := 2*Modelica.Constants.pi*f;
end from_Hz;

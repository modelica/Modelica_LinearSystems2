within Modelica_LinearSystems2.Records;
record SetParameter "Set value of one parameter from a translated model"

  String Name "Name of parameter" annotation (Dialog);
  Real Value "Value of parameter" annotation (Dialog);

  annotation (
    Dialog(
      __Dymola_importDsin(
        button="select" "Select the model parameters to be included in this table",
        onlyStart=true,
        fields(
          Name=initialName,
          Value=initialValue.value))),
  Icon(graphics={
        Rectangle(
          extent={{-100,-30},{100,-90}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-40,-30},{-40,80}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{20,28},{-20,0},{20,-28}},
          color={0,0,0},
          smooth=Smooth.None,
          origin={-40,-10},
          rotation=90)}),
    Documentation(info="<html>
<p>
Set a&nbsp;name of whatever parameter and its value. The record can be vectorized to
set more parameters and values at once.
</p>
</html>"));

end SetParameter;

within Modelica_LinearSystems2.Records;
record ParameterVariation
  "Define variation of one parameter in a given range and optionally select the parameter from a translated model"

  String Name "Name of parameter" annotation (Dialog);
  Modelica_LinearSystems2.Utilities.Types.Grid grid=Modelica_LinearSystems2.Utilities.Types.Grid.OneValue "Definition of parameter grid" annotation (Dialog);

  Real Value=0 "Value of parameter" annotation (Dialog(group="if grid = OneValue"));
  Real Min=-1e100 "Minimum value of parameter"
    annotation(Dialog(group="if grid = Equidistant or Logarithmic"));
  Real Max=1e100 "Maximum value of parameter"
    annotation(Dialog(group="if grid = Equidistant or Logarithmic"));
  Integer nPoints(min=2) = 11
    "Number of parameter values in the range Min .. Max" annotation(Dialog(group="if grid = Equidistant or Logarithmic"));

  annotation (
     Dialog(
       __Dymola_label="Model parameters",
       __Dymola_importDsin(
         button="Select" "Select the model parameters to be included in this table",
         onlyStart=true,
    fields(Name=initialName,
           Value=initialValue.value,
           Min=initialValue.minimum,
           Max=initialValue.maximum,
           Unit=initialValue.unit))),
  Icon(graphics={
        Rectangle(
          extent={{-100,-30},{100,-90}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-100,-30},{-100,44}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{100,-32},{100,44}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-100,40},{100,40}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-60,68},{-100,40},{-60,12}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{60,68},{100,40},{60,12}},
          color={0,0,0},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p>
Set a&nbsp;name of whatever parameter and values needed for parameter variation.
The record can be vectorized to set more parameters and values at once.
</p>
</html>"));

end ParameterVariation;

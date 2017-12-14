within Modelica_LinearSystems2.WorkInProgress.Internal;
record ModelParameters
  "Model parameters to be selected after translation of a model"
  String parName="" "Name of parameter" annotation (Dialog);
  Real parValue=1 "Value of parameter" annotation (Dialog(
      __Dymola_absoluteWidth=20));
  Real parMin = 1 "Minimum value of parameter" annotation (Dialog);
  Real parMax = 6 "Maximum value of parameter" annotation (Dialog);
  Integer nVar(min=2) = 10 "Number of variations parameter (min=2)" annotation (Dialog);
  String parUnit="" "Unit of parameter" annotation (Dialog);
//  annotation(Dialog(__Dymola_importDsin(onlyStart=true,
//    fields(name=initialName))));
// __Dymola_importDsin(
// onlyStart=true,
// fields(
// name=initialName,
// Value=initialValue.value,
// min=initialValue.minimum,
// max=initialValue.maximum,
// fixed=initialValue.fixed,
// unit=initialValue.unit,
// ValueType=initialValue.Type,
// Description=initialValue.description))
  annotation (Icon(graphics={
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
          smooth=Smooth.None)}));
end ModelParameters;

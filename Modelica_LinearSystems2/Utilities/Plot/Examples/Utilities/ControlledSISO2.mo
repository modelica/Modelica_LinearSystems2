within Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities;
model ControlledSISO2 "Controlled SISO system with parameters different to ControlledSISO1"
  extends ControlledSISO1(zerosAndPoles(
    system(
      d2=[20,101; 22,122],
      d1=fill(0.0, 0),
      n1={4},
      n2=[2,2])));
  annotation (
    Documentation(
      info="<html>
<p>
Utility model in order to demonstrate the plotting of a root locus:
</p>

<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/ControlledSISO1a.png\"/>
</div>

<p>
This model consists of a linear, time invariant single-input, single-output plant \"zerosAndPoles\"
</p>
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/ControlledSISO2c.png\"/>
</div>
<p>
that is controlled by a P-Controller with a constant gain k. The pole/zero pattern of the plant
is shown in the next figure:
</p>

<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/ControlledSISO2b.png\"/>
</div>
</html>"));
end ControlledSISO2;

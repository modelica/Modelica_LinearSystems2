within Modelica_LinearSystems2.Utilities.Plot.Records;
record DefaultDiagram
  "Diagram with no curves (might be used as default diagram)"
  extends Diagram(curve=fill(Curve(
                             x=fill(0.0,0),y=fill(0.0,0)),0));

  annotation (Documentation(info="<html>
<p>
This record defines a default diagram without curves. In several cases a default is needed in order
that no translation error occurs and then this record can be used.
</p>
</html>"));
end DefaultDiagram;

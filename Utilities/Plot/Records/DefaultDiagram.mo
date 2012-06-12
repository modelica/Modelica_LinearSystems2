within Modelica_LinearSystems2.Utilities.Plot.Records;
record DefaultDiagram
  "Diagram with no curves (might be used as default diagram)"
  extends Diagram(curve=fill(Curve(
                             x=fill(0.0,0),y=fill(0.0,0)),0));
end DefaultDiagram;

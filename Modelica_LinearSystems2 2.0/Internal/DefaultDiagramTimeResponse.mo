within Modelica_LinearSystems2.Internal;
record DefaultDiagramTimeResponse "Default diagram for a time response plot"
  import Modelica;

   extends Modelica.Icons.Record;
   extends Modelica_LinearSystems2.Utilities.Plot.Records.DefaultDiagram(
                      heading="time response",
                      heightRatio=0.4,
                      legend=false,
                      xLabel="time [s]",
                      yLabel="system output");
  annotation (Documentation(info="<html>
<p>
This record contains the default diagram options for pole/zero plots.
</p>
</html>"));
end DefaultDiagramTimeResponse;

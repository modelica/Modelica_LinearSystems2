within Modelica_LinearSystems2.Internal;
record DefaultDiagramPolesAndZeros "Default diagram for polesAndZeros plot"
  import Modelica;

   extends Modelica.Icons.Record;
   extends Modelica_LinearSystems2.Utilities.Plot.Records.DefaultDiagram(
                      heading="Poles (x) and invariant zeros (o)",
                      heightRatio=0.6,
                      legend=false,
                      xLabel="Real part",
                      yLabel="Imaginary part");
  annotation (Documentation(info="<html>
<p>
This record contains the default diagram options for pole/zero plots.
</p>
</html>"));
end DefaultDiagramPolesAndZeros;

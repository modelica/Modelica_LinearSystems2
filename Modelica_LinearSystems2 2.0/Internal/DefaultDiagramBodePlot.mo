within Modelica_LinearSystems2.Internal;
record DefaultDiagramBodePlot "Default diagram for Bode plot"
  import Modelica;

   extends Modelica.Icons.Record;
   extends Modelica_LinearSystems2.Utilities.Plot.Records.DefaultDiagram(
                      heading="Bode-Diagram",
                      heightRatio=0.4,
                      legend=false,
                      xLabel="Frequency [Hz]",
                      yLabel="",
                      logX=true,
                      logY=true);
  annotation (Documentation(info="<html>
<p>This record contains the default diagram options for magnitude in bode plots. </p>
</html>"));
end DefaultDiagramBodePlot;

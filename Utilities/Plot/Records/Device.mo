within Modelica_LinearSystems2.Utilities.Plot.Records;
record Device "Properties of a device"
  extends Modelica.Icons.Record;

  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm xTopLeft=0
    "Horizontal position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog);
  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm yTopLeft=0
    "Vertical position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog);
  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm diagramWidth=140
    "Width of diagram" annotation(Dialog);
  Modelica_LinearSystems2.Utilities.Plot.Types.ImageResolution_dpi windowResolution=96
    "[dpi] Image resolution in window if applicable (e.g. unscaled window)"   annotation(Dialog);

  Boolean autoLineColor = true
    "If automatic line properties: distinguish curves by color otherwise by line style"
       annotation(Dialog,choices(__Dymola_checkBox=true));

  annotation (Documentation(info="<html>
<p>
With this record the placement of a diagram in the drawing window is defined,
as well as the diagram width (the height to width ratio of a diagram is defined
with record
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Diagram\">Modelica_LinearSystems2.Utilities.Plot.Records.Diagram</a>.
Furthermore the window resolution is defined here in order to compute the number of pixels from length definitions in
mm.
</p>
</html>"));
end Device;

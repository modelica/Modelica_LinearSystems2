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
    "if automatic line properties: distinguish curves by color otherwise by line style"
       annotation(Dialog,choices(__Dymola_checkBox=true));
end Device;

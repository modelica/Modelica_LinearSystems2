within Modelica_LinearSystems2.Utilities.Plot.Types;
type DrawingUnit_mm "Drawing unit specifing the nominal size in [mm]"
   extends Modelica.Icons.TypeReal(final quantity="Length", final unit="mm");
  annotation (Documentation(info="<html>
<p>
All size information for plotting, such as width or height
of a window or the thickness of a line, are defined by type <b>DrawingUnit_mm</b>. The
recommended interpretation is that the DrawingUnit is the unscaled size in a document or
on printer in [mm]. For example, if the width of a diagram is 120,
and the diagram is pasted into a Word or PowerPoint document,
then the width of the diagram in the document is 120 mm.
</p>
</html>"));
end DrawingUnit_mm;

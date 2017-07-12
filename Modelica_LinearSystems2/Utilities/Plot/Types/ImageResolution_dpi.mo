within Modelica_LinearSystems2.Utilities.Plot.Types;
type ImageResolution_dpi
  "Resolution of image in pixel per inch (screen) or dots per inch (printer)"
   extends Modelica.Icons.TypeReal(final quantity="ImageResolution");

  annotation (Documentation(info="<html>
<p>\"ImageResolution_dpi\" defines the mapping of a length coordinate to the resolution of the output device. The resolution [dpi] is defined as \"dots-per-inch\" and therefore a length L_mm defined in [mm] is mapped to a length L_dot in dots (or pixel) with the formula: </p>

<pre>   L_dot = round(ImageResolution_dpi/25.4 * L_mm)
</pre>

<p>where function round(..) rounds to the nearest integer. Typical values are \"96 dpi\" (for screen) or \"600 dpi\" for printer. For example if an \"ImageResolution = 96 dpi\" shall be used for a screen, then 1 mm is mapped to 4 pixel. </p>
</html>"));
end ImageResolution_dpi;

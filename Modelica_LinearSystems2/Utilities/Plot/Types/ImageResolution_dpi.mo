within Modelica_LinearSystems2.Utilities.Plot.Types;
type ImageResolution_dpi
  "Resolution of image in pixel per inch (screen) or dots per inch (printer)"
   extends Modelica.Icons.TypeReal(final quantity="ImageResolution");

  annotation (Documentation(info="<html>
<p>&quot;ImageResolution_dpi&quot; defines the mapping of a length coordinate to the resolution of the output device. The resolution [dpi] is defined as &quot;dots-per-inch&quot; and therefore a length L_mm defined in [mm] is mapped to a length L_dot in dots (or pixel) with the formula: </p>

<pre>   L_dot = round(ImageResolution_dpi/25.4 * L_mm)
</pre>

<p>where function round(..) rounds to the nearest integer. Typical values are &quot;96 dpi&quot; (for screen) or &quot;600 dpi&quot; for printer. For example if an &quot;ImageResolution = 96 dpi&quot; shall be used for a screen, then 1&nbsp;mm is mapped to 4 pixel. </p>
</html>"));
end ImageResolution_dpi;

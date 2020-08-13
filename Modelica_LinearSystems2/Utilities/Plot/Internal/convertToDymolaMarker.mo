within Modelica_LinearSystems2.Utilities.Plot.Internal;
function convertToDymolaMarker
  "Convert from marker type of Modelica_LinearSystems2 to Dymola marker"
  import Symbol = Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol;

  input Symbol symbol
    "Enumeration value of point symbol used in Modelica_LinearSystems2";
  output MarkerStyle result
    "Enumeration value of point symbol used in Dymola";
algorithm
  result :=
    if symbol == Symbol.None         then MarkerStyle.None else
    if symbol == Symbol.Cross        then MarkerStyle.Cross else
    if symbol == Symbol.Circle       then MarkerStyle.Circle else
    if symbol == Symbol.Square       then MarkerStyle.Square else
    if symbol == Symbol.FilledSquare then MarkerStyle.FilledSquare else
    if symbol == Symbol.TriangleDown then MarkerStyle.TriangleDown else
    if symbol == Symbol.TriangleUp   then MarkerStyle.TriangleUp else
    if symbol == Symbol.Diamond      then MarkerStyle.Diamond else
    if symbol == Symbol.Dot          then MarkerStyle.Dot else
    MarkerStyle.None;

  annotation (Inline=true, Documentation(info="<html>
<p>
This function converts the point symbol enumeration as used in the Modelica_LinearSystems2 library
to the enumeration used in Dymola.
</p>
</html>"));
end convertToDymolaMarker;

within Modelica_LinearSystems2.Utilities.Plot;
package Internal "Internal functions, that should not be utilized by a user"
  extends Modelica.Icons.Package;
  function convertToDymolaMarker
    "Convert from marker type of Modelica_LinearSystems2 to Dymola marker"
    import Symbol = Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol;
     input Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol symbol;
     output MarkerStyle result;
  algorithm
     result :=if symbol == Symbol.None then
                 MarkerStyle.None else
              if symbol == Symbol.Cross then
                 MarkerStyle.Cross else
              if symbol == Symbol.Circle then
                 MarkerStyle.Circle else
              if symbol == Symbol.Square then
                 MarkerStyle.Square else
              if symbol == Symbol.FilledSquare then
                 MarkerStyle.FilledSquare else
              if symbol == Symbol.TriangleDown then
                 MarkerStyle.TriangleDown else
              if symbol == Symbol.TriangleUp then
                 MarkerStyle.TriangleUp else
              if symbol == Symbol.Diamond then
                 MarkerStyle.Diamond else MarkerStyle.None;

     annotation(Inline=true);
  end convertToDymolaMarker;

  function convertToDymolaPattern
    "Convert from pattern type of Modelica_LinearSystems2 to Dymola pattern"
    import Pattern = Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern;
     input Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern pattern;
     output LinePattern result;
  algorithm
     result:=if pattern == Pattern.None then
                 LinePattern.None else
              if pattern == Pattern.Solid then
                 LinePattern.Solid else
              if pattern == Pattern.Dash then
                 LinePattern.Dash else
              if pattern == Pattern.Dot then
                 LinePattern.Dot else
              if pattern == Pattern.DashDot then
                 LinePattern.DashDot else
              if pattern == Pattern.DashDotDot then
                 LinePattern.DashDotDot else LinePattern.Solid;

     annotation(Inline=true);
  end convertToDymolaPattern;
end Internal;

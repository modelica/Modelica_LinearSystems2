within Modelica_LinearSystems2.Utilities.Plot.Internal;
function convertToDymolaPattern
  "Convert from pattern type of Modelica_LinearSystems2 to Dymola pattern"
  extends Modelica.Icons.Function;

  import Pattern = Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern;

  input Pattern pattern
    "Enumeration value of line pattern used in Modelica_LinearSystems2";
  output LinePattern result
    "Enumeration value of line pattern used in Dymola";
algorithm
  result :=
    if pattern == Pattern.None       then LinePattern.None else
    if pattern == Pattern.Solid      then LinePattern.Solid else
    if pattern == Pattern.Dash       then LinePattern.Dash else
    if pattern == Pattern.Dot        then LinePattern.Dot else
    if pattern == Pattern.DashDot    then LinePattern.DashDot else
    if pattern == Pattern.DashDotDot then LinePattern.DashDotDot else
    LinePattern.Solid;

  annotation (Inline=true, Documentation(info="<html>
<p>
This function converts the line pattern enumeration as used in the Modelica_LinearSystems2 library
to the enumeration used in Dymola.
</p>
</html>"));
end convertToDymolaPattern;

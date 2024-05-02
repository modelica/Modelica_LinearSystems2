within Modelica_LinearSystems2.Internal;
function printHTML_end "Print end of html file"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams;

  input String fileName="??????.html"
    "File on which the html end should be written";
algorithm
  Streams.print("</html>", fileName);
end printHTML_end;

within Modelica_LinearSystems2.Internal;
function printHTML_start "Print start of html file"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;

  input String fileName="??????.html"
    "File on which the html start should be written";
algorithm
   Modelica.Utilities.Files.removeFile(fileName);
   // Following doesn't work in Dymola
   //print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", fileName);
   print("<html>", fileName);
   print("<style type=\"text/css\">", fileName);
   print("* { font-size: 10pt; font-family: Arial,sans-serif; }", fileName);
   print("</style>", fileName);
end printHTML_start;

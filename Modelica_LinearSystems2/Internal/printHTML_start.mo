within Modelica_LinearSystems2.Internal;
function printHTML_start "Print start of html file"

  import Modelica.Utilities.Files;
  import Modelica.Utilities.Streams;

  input String fileName="??????.html"
    "File on which the html start should be written";
algorithm
   Files.removeFile(fileName);
   // Following doesn't work in Dymola
   //Streams.print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", fileName);
   Streams.print("<html>", fileName);
   Streams.print("<style type=\"text/css\">", fileName);
   Streams.print("* { font-size: 10pt; font-family: Arial,sans-serif; }", fileName);
   Streams.print("</style>", fileName);
end printHTML_start;

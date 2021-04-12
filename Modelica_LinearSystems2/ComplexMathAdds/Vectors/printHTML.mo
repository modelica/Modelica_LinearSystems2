within Modelica_LinearSystems2.ComplexMathAdds.Vectors;
function printHTML
  "Print complex vector as HTML in sorted form (vector is expected to have pure real and/or conjugate complex values (so poles and/or zeros)"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Files;
  import Modelica.Utilities.Streams;
  import Modelica.Utilities.Streams.print;

  input Complex c[:] "Complex vector to be printed";
  input String heading="Zeros" "Heading above the table";
  input String name="zeros" "Heading of value column";
  input Boolean sort=true
    "= true, if values are sorted, otherwise no sorting";
  input Boolean ascending = true
    "= true if ascending order, otherwise descending order for sorting";
  input Boolean sortFrequency=true
    "= true, if sorting is first for imaginary then for real value, otherwise sorting is for absolute value";
protected
  Integer nc = size(c,1);
  Complex cSorted[size(c,1)];
  Integer cIndex[size(c,1)];
  String tempFile = "TemporaryForPrint.html";
  Integer nReal;

  function printTable
    "Print the table with eigenvalues in html format on file"
    extends Modelica.Icons.Function;

    import Modelica.Utilities.Strings;
    import Modelica.Utilities.Streams.print;
    import Complex;

    input Complex systemZeros[:];
    input Integer nReal;
    input String tempFile;
  protected
    Integer nz=size(systemZeros, 1);

    String number;
    Real timeConstant;
    Real freq;
    Real damp;
    Integer i;

  algorithm
    for i in 1:nReal loop
      // Build eigenvalue number

      number := String(
              i,
              minimumLength=7,
              leftJustified=false);
      timeConstant := if abs(systemZeros[i].re) > 10*Modelica.Constants.eps
         then 1/abs(systemZeros[i].re) else 1/(10*Modelica.Constants.eps);

      print(
        "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
         + number + " </td>" + "\n  <td> &nbsp; " + String(systemZeros[i].re,
        format="14.4e") + " </td>" + "\n  <td> &nbsp; " + String(
        timeConstant, format="9.4f") + " </td>" +
        "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
        "\n  <td style=\"text-align:center\"> &nbsp; --- </td>\n</tr>",
        tempFile);

    end for;

    for i in nReal + 1:2:nz loop
      number := String(i) + "/" + String(i + 1);
      number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;

      // Determine frequency and number of corresponding zero
      (freq,damp) := Modelica_LinearSystems2.ComplexMathAdds.frequency(systemZeros[i]);

      print(
        "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
         + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
         + String(systemZeros[i].re, format="14.4e") + " &plusmn; " +
        String(systemZeros[i].im, format="12.4e") + "j </td>" +
        "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
        "\n  <td style=\"text-align:left\"> &nbsp; " + String(freq, format=
        "9.4f") + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
         + String(damp, format="9.4f") + " </td>\n</tr>", tempFile);

    end for;

    print("</table>\n", tempFile);
  end printTable;
algorithm
  if size(c,1) < 1 then
     return;
  end if;

  // Sort complex vector
  if sort then
    (cSorted,cIndex) := Modelica.ComplexMath.Vectors.sort(c,ascending,sortFrequency);
  else
    cSorted :=c;
    cIndex :=1:nc;
  end if;

  // Number of real zeros
  nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(cSorted);

  // Remove temporary file, if it exists
  Files.removeFile(tempFile);

  // Print heading
  // Following doesn't work in Dymola
  //Streams.print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", tempFile);
  print("<html>", tempFile);
  print("<style type=\"text/css\">", tempFile);
  print("* { font-size: 10pt; font-family: Arial,sans-serif; }", tempFile);
  print("</style>", tempFile);

  print("<table style=\"background-color:rgb(100, 100, 100); margin:20px 0 20px 20px;\" "
         + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", tempFile);
  print("<caption>" + heading + "</caption>", tempFile);
  print("<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">" +
        "\n  <td> number </td>\n  <td>" + name + "</td>\n  <td> time constant [s] </td>" +
        "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n</tr>", tempFile);

  // Print values of complex vector
  printTable(cSorted, nReal, tempFile);

  // Print end of file
  print("</html>", tempFile);

  // Read file into output window and remove temporary file
  Streams.readFile(tempFile);
  Files.removeFile(tempFile);
end printHTML;

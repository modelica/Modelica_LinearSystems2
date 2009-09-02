within Modelica_LinearSystems2.Internal.Streams;
function stateSpaceString_html
  import Modelica;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica.Utilities.Strings;

  input StateSpace ss
    "State space system to be transformed in a String representation";
  input Integer significantDigits=12
    "Number of significant digits that are shown";
  input String name="ss" "Independent variable name used for printing";
  output String s="";
protected
  String space=Strings.repeat(5);
  String space2=Strings.repeat(3);
  String space3=Strings.repeat(8);
  Integer nx=size(ss.A, 1);
  Integer nu=size(ss.B, 2);
  Integer ny=size(ss.C, 1);
  Integer sizeD=size(ss.D, 2);
  Integer stringMaxLength;
  Boolean namesExist=false;

algorithm
//Checking if name arrays are empty
  for i in 1:nx loop
    namesExist := namesExist or (ss.xNames[i] <> "");
  end for;

  for i in 1:ny loop
    namesExist := namesExist or (ss.yNames[i] <> "");
  end for;

  for i in 1:nu loop
    namesExist := namesExist or (ss.uNames[i] <> "");
  end for;

  if namesExist then
    Modelica.Utilities.Streams.print("namesExist == true");
  else
    Modelica.Utilities.Streams.print("namesExist == false");
  end if;

  stringMaxLength := max(size(ss.xNames, 1), min(size(ss.yNames, 1),
    11));
//, size(ss.uNames,1)));
                                                               //, max(size(ss.yNames), size(ss.uNames))

  if nx == 0 and sizeD == 0 then
    s := name + ".A = []<br>  " + name + ".B = []<br>   " + name + ".C = [] <br>   " + name + ".D = []";
  else
    s := "<br>" + name + ".A = <br>";

//Horizontal
// Two alternatives when printing state names
    if namesExist == false then
//    if ss.xNames == fill(" ", size(A, 2)) then

      s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
        "x1 ";
    else
      s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
        ss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[
        1]), 11)) + " " + Strings.substring(
            ss.xNames[1],
            1,
            min(Strings.length(ss.xNames[1]), 11));
    end if;

    for i in 2:nx loop

//Two alternatives when printing state names

      if namesExist == false then
        s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
           + String(i - 1))) + "x" + String(i) + " ";
      else
        s := s + " " + Strings.repeat(significantDigits + 11 - min(
          Strings.length(ss.xNames[i - 1]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11));

      end if;

    end for;
    s := s + "<br>";

    for i in 1:nx loop
//Vertical
//Two alternatives when printing state names
      if namesExist == false then
        s := s + space + "x" + String(i) + " ";
      else
        s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
          ss.xNames[i]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11)) + " ";

      end if;

      for j in 1:nx loop
        if ss.A[i, j] >= 0 then
          s := s + " ";
        end if;
        s := s + String(ss.A[i, j], significantDigits=significantDigits) +
          Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
          ss.A[i, j]), significantDigits=significantDigits)));
      end for;
      s := s + "<br>";
    end for;
//--------------------------------------------------------------------------------------------------------------------------------------------------
    s := s + "<br>" + name + ".B = <br>";
 //Horizontal
// Two alternatives when printing state names
    if namesExist == false then
      s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
        "u1 ";
    else
      s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
        ss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[
        1]), 11)) + " " + Strings.substring(
            ss.uNames[1],
            1,
            min(Strings.length(ss.uNames[1]), 11));
    end if;

    for i in 2:nu loop
//Two alternatives when printing state names
      if namesExist == false then
        s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
           + String(i - 1))) + "u" + String(i) + " ";
      else
        s := s + " " + Strings.repeat(significantDigits + 11 - min(
          Strings.length(ss.uNames[i - 1]), 11)) + Strings.substring(
              ss.uNames[i],
              1,
              min(Strings.length(ss.uNames[i]), 11));
      end if;
    end for;
    s := s + "<br>";
    for i in 1:nx loop

//Vertical
//Two alternatives when printing state names
      if namesExist == false then
        s := s + space + "x" + String(i) + " ";
      else
        s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
          ss.xNames[i]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11)) + " ";

      end if;

      for j in 1:nu loop
        if ss.B[i, j] >= 0 then
          s := s + " ";
        end if;
        s := s + String(ss.B[i, j], significantDigits=significantDigits) +
          Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
          ss.B[i, j]), significantDigits=significantDigits)));
      end for;
      s := s + "<br>";
    end for;

//--------------------------------------------------------------------------------------------------------------------------------------------------
    s := s + "<br>" + name + ".C = <br>";
 //Horizontal
// Two alternatives when printing state names
    if namesExist == false then
      s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
        "x1 ";
    else
      s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
        ss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[
        1]), 11)) + " " + Strings.substring(
            ss.xNames[1],
            1,
            min(Strings.length(ss.xNames[1]), 11));
    end if;

    for i in 2:nx loop
//Two alternatives when printing state names
      if namesExist == false then
        s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
           + String(i - 1))) + "x" + String(i) + " ";
      else
        s := s + " " + Strings.repeat(significantDigits + 11 - min(
          Strings.length(ss.xNames[i - 1]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11));
      end if;
    end for;
    s := s + "<br>";

    for i in 1:nu loop
//Vertical
//Two alternatives when printing state names
      if namesExist == false then
        s := s + space + "y" + String(i) + " ";
      else
        s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
          ss.yNames[i]), 11)) + Strings.substring(
              ss.yNames[i],
              1,
              min(Strings.length(ss.yNames[i]), 11)) + " ";

      end if;

      for j in 1:nx loop
        if ss.C[i, j] >= 0 then
          s := s + " ";
        end if;
        s := s + String(ss.C[i, j], significantDigits=significantDigits) +
          Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
          ss.C[i, j]), significantDigits=significantDigits)));
      end for;
      s := s + "<br>";
    end for;
//--------------------------------------------------------------------------------------------------------------------------------------------------
    s := s + "<br>" + name + ".D = <br>";
 //Horizontal
// Two alternatives when printing state names
    if namesExist == false then
      s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
        "u1 ";
    else
      s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
        ss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[
        1]), 11)) + " " + Strings.substring(
            ss.uNames[1],
            1,
            min(Strings.length(ss.uNames[1]), 11));
    end if;

    for i in 2:nu loop
//Two alternatives when printing state names
      if namesExist == false then
        s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
           + String(i - 1))) + "u" + String(i) + " ";
      else
        s := s + " " + Strings.repeat(significantDigits + 11 - min(
          Strings.length(ss.uNames[i - 1]), 11)) + Strings.substring(
              ss.uNames[i],
              1,
              min(Strings.length(ss.uNames[i]), 11));
      end if;
    end for;
    s := s + "<br>";
    for i in 1:ny loop
//Vertical
//Two alternatives when printing state names
      if namesExist == false then
        s := s + space + "y" + String(i) + " ";
      else
        s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
          ss.yNames[i]), 11)) + Strings.substring(
              ss.yNames[i],
              1,
              min(Strings.length(ss.yNames[i]), 11)) + " ";

      end if;

      for j in 1:nu loop
        if ss.D[i, j] >= 0 then
          s := s + " ";
        end if;
        s := s + String(ss.D[i, j], significantDigits=significantDigits) +
          Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
          ss.D[i, j]), significantDigits=significantDigits)));
      end for;
      s := s + "<br>";
    end for;

  end if;

end stateSpaceString_html;

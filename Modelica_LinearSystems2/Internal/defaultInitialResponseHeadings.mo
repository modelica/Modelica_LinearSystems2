within Modelica_LinearSystems2.Internal;
function defaultInitialResponseHeadings
  "Compute default headings for initial response plot"
  extends Modelica.Icons.Function;

  input Integer sizeC1;
  output String heading[sizeC1,1];
algorithm
  heading[:, :] := matrix({"Initial response (Out " + String(i) + ")" for i in
    1:sizeC1});
end defaultInitialResponseHeadings;

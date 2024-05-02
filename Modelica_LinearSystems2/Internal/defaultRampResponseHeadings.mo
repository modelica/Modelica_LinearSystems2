within Modelica_LinearSystems2.Internal;
function defaultRampResponseHeadings
  "Compute default headings for ramp response plot"
  extends Modelica.Icons.Function;

  input Integer sizeB2;
  input Integer sizeC1;
  output String heading[sizeC1,sizeB2];
algorithm
  heading[:, :] := {{"Ramp response[" + String(i) + "," + String(j) + "]" for j in
        1:sizeB2} for i in 1:sizeC1};
end defaultRampResponseHeadings;

within Modelica_LinearSystems2.Internal;
function defaultStepResponseHeadings
  "Compute default headings for step response plot"

  input Integer sizeB2;
  input Integer sizeC1;
  output String heading[sizeC1,sizeB2];
algorithm
  heading[:, :] := {{"Step response[" + String(i) + "," + String(j) + "]" for j in
        1:sizeB2} for i in 1:sizeC1};
end defaultStepResponseHeadings;

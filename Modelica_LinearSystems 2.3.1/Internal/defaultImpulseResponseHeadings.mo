within Modelica_LinearSystems2.Internal;
function defaultImpulseResponseHeadings
  "Compute default headings for impulse response plot"

  input Integer sizeB2;
  input Integer sizeC1;
  output String heading[sizeC1,sizeB2];
algorithm
  heading[:, :] := {{"Impulse response[" + String(i) + "," + String(j) + "]"
    for j in 1:sizeB2} for i in 1:sizeC1};
end defaultImpulseResponseHeadings;

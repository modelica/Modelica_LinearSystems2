within Modelica_LinearSystems2.Internal;
function defaultColumnLabels "Compute default column labels for plot"
  import Modelica_LinearSystems2.StateSpace;

  input Integer sizeC1;
  output String columnLabels[sizeC1 + 1];

algorithm
  columnLabels := cat(
    1,
    {"time [t]"},
    {"Out[" + String(i) + "]" for i in 1:sizeC1}) "Column labels";
end defaultColumnLabels;

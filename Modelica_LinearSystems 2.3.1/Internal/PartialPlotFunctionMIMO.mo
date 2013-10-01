within Modelica_LinearSystems2.Internal;
partial function PartialPlotFunctionMIMO
  "Interface of a plot function for MIMO systems"

  import Modelica_LinearSystems2.Utilities.Plot;

  input Plot.Records.Diagram defaultDiagram "Default diagram layout"               annotation(Dialog);
  input Plot.Records.Device device=Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);
end PartialPlotFunctionMIMO;

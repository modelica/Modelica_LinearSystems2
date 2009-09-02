within Modelica_LinearSystems2.Internal;
partial function PartialPlotFunction "Interface of a plot function"

  import Modelica_LinearSystems2.Utilities.Plot;

  input Plot.Records.Diagram defaultDiagram "Default diagram layout" annotation(Dialog);
  input Plot.Records.Device device=Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);
end PartialPlotFunction;

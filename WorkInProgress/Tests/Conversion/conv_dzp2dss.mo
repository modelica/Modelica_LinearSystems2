within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function conv_dzp2dss

  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.DiscreteTransferFunction;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  input Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;
  input Integer n=0;
protected
  DiscreteStateSpace dss1=DiscreteStateSpace(dzp);
  Modelica_LinearSystems2.DiscreteTransferFunction dtf1=
                                DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss1);
  Modelica_LinearSystems2.DiscreteTransferFunction dtf2=
                                DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp);
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp2=
                             Modelica_LinearSystems2.DiscreteZerosAndPoles(
                                                   dtf2);
  DiscreteStateSpace dss2=DiscreteStateSpace(dzp2);
  Modelica_LinearSystems2.DiscreteTransferFunction dtf3=
                                DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss2);

algorithm
  if n>0 then
  Modelica.Utilities.Streams.print("\nNr. " + String(n));
  end if;

  Modelica.Utilities.Streams.print("\ndzp = " + String(dzp));
  Modelica.Utilities.Streams.print("dtf1 = " + String(dtf1));
  Modelica.Utilities.Streams.print("dtf2 = " + String(dtf2));
  Modelica.Utilities.Streams.print("\ndzp2 = " + String(dzp2));
  Modelica.Utilities.Streams.print("dtf3 = " + String(dtf3));
  Modelica.Utilities.Streams.print("dtf2 = " + String(dtf2));

end conv_dzp2dss;

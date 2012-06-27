within Modelica_LinearSystems2.WorkInProgress.Tests.temp;
partial function baseFunction
  extends Modelica.Icons.Function;

  input Real u1[2];
  input Real u2[:];

  output Real y1[size(u1,1)];
  output Real y2[size(u2,1)];

end baseFunction;

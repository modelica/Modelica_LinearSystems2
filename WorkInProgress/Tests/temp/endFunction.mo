within Modelica_LinearSystems2.WorkInProgress.Tests.temp;
function endFunction

  input Modelica_LinearSystems2.WorkInProgress.Tests.temp.baseFunction f;
  input Real u1[2];
  input Real u2[2];

  output Real y1[size(u1,1)];
  output Real y2[size(u2,1)];

algorithm
  (y1,y2):=f(u1,u2);

end endFunction;

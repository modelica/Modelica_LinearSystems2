within Modelica_LinearSystems2.Internal;
function timeResponseSamples
  "Estimate reasonable discretisation sample time and simulation time span for time response plot"
  import Modelica_LinearSystems2;
  import Complex;

  input Modelica_LinearSystems2.StateSpace sc;
  output Real dt "Sample time";
  output Real tSpan "Time span";
protected
  Complex eig[size(sc.A, 1)];
  Real realp[size(sc.A, 1)];
  Real sorted[size(sc.A, 1)];
  Real indices[size(sc.A, 1)];
  Integer i;
algorithm
  if size(sc.A,1) > 0 then
     eig := Modelica_LinearSystems2.StateSpace.Analysis.eigenValues(sc);
     for i in 1:size(eig, 1) loop
       realp[i] := eig[i].re;
     end for;
     (sorted,indices) := Modelica.Math.Vectors.sort(realp);

     // Estimate simulation time span
     if sorted[end] < 0 then
       tSpan := -5/sorted[end];
     elseif sorted[end] > 0 then
       tSpan := 15/sorted[end];
     elseif sorted[end] == 0 then
       tSpan := 15000;
     end if;
     // Curb simulation time span to maximal 15000s
     if tSpan > 15000 then
       tSpan := 15000;
     end if;

     // set sample time:
     dt := tSpan/1000;
  else
    tSpan := 1.0;
    dt := 0.01;
  end if;
end timeResponseSamples;

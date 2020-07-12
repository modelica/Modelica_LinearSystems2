within Modelica_LinearSystems2.ComplexMathAdds;
function maxElement "Return maximum element of complex vector"
  extends Modelica.Icons.Function;

  input Complex v[:] "Vector";
  output Real result "Element of v with largest absolute value";
  output Integer index "v[index] has the largest absolute value";

protected
  Real absv_i;
algorithm
  if size(v,1) > 0 then
     result := Modelica.ComplexMath.abs(v[1]);
     index  := 1;
     for i in 2:size(v,1) loop
        absv_i :=Modelica.ComplexMath.abs(v[i]);
        if absv_i > result then
           result := absv_i;
           index := i;
        end if;
     end for;
  else
     result := 0;
     index  := 0;
  end if;
  annotation (Documentation(info="<html>
</html>"));
end maxElement;

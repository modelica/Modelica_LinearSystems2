within Modelica_LinearSystems2.Math.Vectors;
encapsulated function find "find element in vector"
  import Modelica;
  import Modelica_LinearSystems2.Math;
  extends Modelica.Icons.Function;
  input Integer s "Search for s";
  input Integer v[:] "Vector";
  output Integer result
    "v[result] = s (first occurrence of s); result=0, if not found";
protected
  Integer i;
algorithm
  result :=0;
  i :=1;
  while i <= size(v,1) loop
     if v[i] == s then
        result :=i;
        i :=size(v, 1) + 1;
     else
        i :=i + 1;
     end if;
  end while;
end find;

within Modelica_LinearSystems2.Math.Vectors;
encapsulated function find "Find element in vector"
  extends Modelica.Icons.Function;
  import Modelica;

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

  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Vectors.find instead");
end find;

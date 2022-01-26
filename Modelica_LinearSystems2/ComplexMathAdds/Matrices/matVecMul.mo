within Modelica_LinearSystems2.ComplexMathAdds.Matrices;
function matVecMul "Multiply a complex matrices with a complex vector"
  extends Modelica.Icons.Function;

  input Complex m[:,:] "Complex matrix";
  input Complex vi[size(m, 2)] "Complex vector";
  output Complex vo[size(m, 1)] "= m*vi";
//   Real Rem[size(m, 1),size(m, 2)]=Re(m);
//   Real Imm[size(m, 1),size(m, 2)]=Im(m);
//   Real Revi[size(vi, 1)]=Re(vi);
//   Real Imvi[size(vi, 1)]=Im(vi);
protected
  Integer l1;

algorithm
//  vo := Rem*Revi - Imm*Imvi + j*(Rem*Imvi + Imm*Revi);

//    for l1 in 1:size(m, 1) loop
//      vo[l1] := Complex(0);
//      for l2 in 1:size(m, 2) loop
//        vo[l1] := vo[l1] + m[l1, l2]*vi[l2];
//      end for;//l2
//    end for;  //l1

for l1 in 1:size(m, 1) loop
  vo[l1] := ComplexMathAdds.Vectors.multiply(m[l1,:],vi);
end for;

end matVecMul;

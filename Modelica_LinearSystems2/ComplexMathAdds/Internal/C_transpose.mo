within Modelica_LinearSystems2.ComplexMathAdds.Internal;
function C_transpose "Computes the transposed matrix of a complex matrix"
  extends Modelica.Icons.Function;
  import Modelica.ComplexMath.j;

  input Complex C[:,:];
  output Complex CT[size(C, 2),size(C, 1)];
protected
  Integer l1;
  Integer l2;
algorithm
//  CT := Complex(1)*transpose(C[:,:].re)-j*transpose(C[:,:].im);// too slow
  for l1 in 1:size(C, 1) loop
    for l2 in 1:size(C, 2) loop
//      CT[l2, l1] := Re(C[l1, l2]) - j*Im(C[l1, l2]);
      CT[l2, l1] := Modelica.ComplexMath.conj(C[l1,l2]);
    end for;
  end for;

end C_transpose;

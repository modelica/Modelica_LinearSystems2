within Modelica_LinearSystems2.ComplexMathAdds.Matrices;
function matMatMul "Multiply two complex matrices"
  extends Modelica.Icons.Function;

  input Complex m1[:,:] "Complex matrix 1";
  input Complex m2[size(m1, 2),:] "Complex matrix 2";
  output Complex m3[size(m1, 1),size(m2, 2)] "= m1*m2";
//   Real m1Re[size(m1,1),size(m1,2)]=Re(m1);
//   Real m1Im[size(m1,1),size(m1,2)]=Im(m1);
//   Real m2Re[size(m2,1),size(m2,2)]=Re(m2);
//   Real m2Im[size(m2,1),size(m2,2)]=Im(m2);
protected
  Integer l1;
  Integer l2;
algorithm
//  m3 := m1[:, :].re*m2[:, :].re - m1[:, :].im*m2[:, :].im + j*(m1[:, :].re*m2[:,:].im + m1[:, :].im*m2[:, :].re);
//  m3 :=m1Re*m2Re - m1Im*m2Im + j*(m1Re*m2Im + m1Im*m2Re);

 for l1 in 1:size(m1,1) loop
   for l2 in 1:size(m2,2) loop
     m3[l1,l2] := ComplexMathAdds.Vectors.multiply(m1[l1,:],m2[:,l2]);
   end for;
   end for;

end matMatMul;

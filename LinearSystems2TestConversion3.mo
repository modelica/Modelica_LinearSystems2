within ;
package LinearSystems2TestConversion3
  package Math
    function complexNumerics
      import Modelica_LinearSystems2.Math.Complex;

    protected
      Complex j = Modelica_LinearSystems2.Math.Complex.j();
      Complex c1=2+3*j;
      Complex c2=3+4*j;
      Complex c3;
      Complex cv[3] "Vector";
      Complex cm[3,2] "Matrix";
    algorithm
      c3 := Complex(2);
      cv := {Complex(2), Complex(1,7), Complex(3,-3)};
      cm := [cv, Modelica_LinearSystems2.Math.Complex.Vectors.reverse(cv)];
      c3 := Modelica_LinearSystems2.Math.Complex.'constructor'(9, -4);
      Complex.'-'.negate(c1);
      Modelica_LinearSystems2.Math.Complex.'-'.subtract(c1, c2);
      Complex.'+'(c1, c2);
      Modelica_LinearSystems2.Math.Complex.'*'(c1, c2);
      Complex.'/'(c1, c2);
      Modelica_LinearSystems2.Math.Complex.'=='(c1, c2);
      Complex.'String'(c3);
      Complex.'abs'(c3);
      Complex.'sqrt'(c3);
      Modelica_LinearSystems2.Math.Complex.'max'(cv);
      Complex.exp(c1);
      Complex.log(c1);
      Complex.sin(c1);
      Complex.cos(c1);
      Complex.arg(c1);
      Complex.conj(c1);
      Modelica_LinearSystems2.Math.Complex.real(c1);
      Modelica_LinearSystems2.Math.Complex.imag(c1);
      Modelica_LinearSystems2.Math.Complex.eigenValues(diagonal({2,3,6}));
      Complex.eigenVectors(diagonal({2,3,6}));
      Modelica_LinearSystems2.Math.Complex.frequency(c1);
      Complex.Vectors.print("c1", c=cv);
      Modelica_LinearSystems2.Math.Complex.Vectors.printHTML(cv);
      Modelica_LinearSystems2.Math.Complex.Vectors.length(cv);
      Complex.Vectors.norm(cv);
      Complex.Vectors.normalize(cv);
      Complex.Vectors.sortComplex(cv);
      Complex.Vectors.multiply(cv,cv);
      Complex.Vectors.reverse(cv);
      Modelica_LinearSystems2.Math.Complex.Matrices.print(cm);
      Complex.Matrices.matMatMul(cm, cm);
      Complex.Matrices.matVecMul(cm, cv);
    end complexNumerics;

    function readMatrixGainTest

    protected
      Real K0[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(m=1,n=1);
      Real K1[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain("filename.mat",m=1,n=1);
      Real K2[1,1] = Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(matrixName="testK",m=1,n=1);
    algorithm
    end readMatrixGainTest;
  end Math;

  package Types

    model Issue13 "Conversion test concerning github issue #13"
      import Modelica_LinearSystems2.Types;
      parameter Modelica_LinearSystems2.Types.Grid grid = Modelica_LinearSystems2.Types.Grid.Equidistant;
      parameter Types.AnalogFilter analogFilter = Modelica_LinearSystems2.Types.AnalogFilter.Chebyshev;
      parameter Modelica_LinearSystems2.Types.FilterType filterType = Types.FilterType.HighPass;
      parameter Modelica_LinearSystems2.Types.Method method = Modelica_LinearSystems2.Types.Method.Trapezoidal;
      parameter Modelica_LinearSystems2.Types.StaircaseMethod staircaseMethod = Modelica_LinearSystems2.Types.StaircaseMethod.QR;
      parameter Types.TimeResponse timeResponse = Types.TimeResponse.Impulse;
      parameter Modelica_LinearSystems2.Types.Window window = Modelica_LinearSystems2.Types.Window.Bartlett;

    end Issue13;
  end Types;

  package Streams

    function readSystemDimensionTest

    protected
      Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension();
      Integer xuy2[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension2();
      Integer xuy3[3] = Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension();
    algorithm

    end readSystemDimensionTest;
  end Streams;
  annotation (uses(Modelica_LinearSystems2(version="2.4.0")));
end LinearSystems2TestConversion3;

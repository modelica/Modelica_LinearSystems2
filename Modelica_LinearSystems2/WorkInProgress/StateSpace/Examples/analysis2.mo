within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function analysis2 "Example to check controllability of a state space system"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Complex;
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;

  input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
      A=[-6.0,0,2,-13.0; -9.0,-12.0,-3.25,19.25; 0,0,1,34.0; 0,0,-34.0,-31.0]/6,
      B=[1; 1; 1; 1],
      C=identity(4),
      D=zeros(4, 1),
      xNames={"x1","x2","x3","x4"},
      uNames={"u"},
      yNames={"y1","y2","y3","y4"});
  input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
      Modelica_LinearSystems2.Internal.AnalyseOptions(
      plotEigenValues=true,
      plotInvariantZeros=true,
      plotStepResponse=true,
      plotFrequencyResponse=false,
      printEigenValues=true,
      printEigenValueProperties=true,
      printInvariantZeros=true,
      printControllability=false,
      printObservability=false,
      headingEigenValues="Eigenvalues",
      headingInvariantzeros="Invariant zeros",
      headingStepResponse="Step response",
      headingFrequencyResponse="Frequency response");

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
   annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
                                                                           annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                     caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"  annotation(Dialog(group="system data definition",enable = systemOnFile));

  output Boolean ok;

protected
  StateSpace ss=if systemOnFile then Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else ssi;
  StateSpace ssModal=ss;
  Complex ceigvec[4,4] "complex eigen values of the system";

  Complex T[4,4]=fill(Complex(0),4,4);
  Complex Vs[4,4];
  Real VsRe[4,4];
  Real Lambda[4,4];
  Real eigvec[4,4] "eigen values of the system";
  Real eigvecRe[4,4]=fill(0, 4, 4) "eigen values of the system";
  Real eigvecIm[4,4]=fill(0, 4, 4) "eigen values of the system";
  Complex eigval[4]=fill(Complex(0), 4)
  annotation (Documentation(info="<html>
This example shows the usage of <b>function Modelica_LinearSystems2.StateSpace.Design.assignPolesMI</b> which is
to design pole assignment controllers for state space systems with multiple input.
</html>"));
  Complex j = Modelica.ComplexMath.j;
  Complex Tp[2,2] = sqrt(0.5)*[1+0*j, -j; 1+0*j, j];
  Integer i;

  Real x0[size(ss.A, 1)]=ones(size(ss.A, 1)) "Initial state vector";
  Real z0[size(ss.A, 1)]=ones(size(ss.A, 1)) "Initial state vector";

algorithm
  ok := false;
  StateSpace.Analysis.analysis(ss, fileName="analysis.html", analyseOptions=analyseOptions, description="Description of the system");
  (eigvec,eigval) := Modelica_LinearSystems2.StateSpace.Analysis.eigenVectors(ss, false);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("eval", eigval);
  Modelica_LinearSystems2.Math.Matrices.printMatrix(eigvec, 6, "leftEigenVectors");

  i := 1;
  while i <= 4 loop
    if abs(Im(eigval[i])) > 1e-12 then
      eigvecRe[:, i] := eigvec[:, i];
      eigvecIm[:, i] := eigvec[:, i + 1];
      eigvecRe[:, i + 1] := eigvec[:, i];
      eigvecIm[:, i + 1] := -eigvec[:, i + 1];
      T[i,i] := Tp[1,1];
      T[i,i+1] := Tp[1,2];
      T[i+1,i] := Tp[2,1];
      T[i+1,i+1] := Tp[2,2];
      i := i + 2;
    else
      eigvecRe[:, i] := eigvec[:, i];
      eigvecIm[:, i] := fill(0, 4);
      T[i,i] := Complex(1);
      i := i + 1;
    end if;
  end while;
  ceigvec := Complex(1)*eigvecRe + j*eigvecIm;
  Vs := ceigvec*T;
  VsRe := Re(Vs);
  Lambda:= Modelica.Math.Matrices.inv(VsRe)*ss.A*VsRe;

  ssModal.A := Lambda;
  ssModal.B := Modelica.Math.Matrices.solve2(VsRe,ss.B);
//  ssModal.C := ss.C*VsRe;
  ssModal.C := ss.C;
  ssModal.D := ss.D;
  z0 := x0;//Modelica.Math.Matrices.solve(VsRe,x0);

  Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(ceigvec), 6, "cevRe");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(ceigvec), 6, "cevIm");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(T), 6, "TRe");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(T), 6, "TIm");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(Vs), 6, "VsRe");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(Vs), 6, "VsIm");
  Modelica_LinearSystems2.Math.Matrices.printMatrix(Lambda, 6, "Lam");

  Modelica_LinearSystems2.StateSpace.Plot.initialResponse(ss=ss, tSpan=4, x0=x0, defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Initial response system states"),subPlots=false);
  Modelica_LinearSystems2.StateSpace.Plot.initialResponse(ss=ssModal, tSpan=4, x0=z0, defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Initial response modal states"),subPlots=false);

  ok := true;
  annotation (
    Documentation(info="<html>
This example shows the usage of <b>function Modelica_LinearSystems2.StateSpace.Analysis.isControllable</b> which is
to check whether a system is controllable or not.
</html>"));
end analysis2;

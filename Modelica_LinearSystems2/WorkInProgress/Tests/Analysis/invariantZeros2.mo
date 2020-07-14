within Modelica_LinearSystems2.WorkInProgress.Tests.Analysis;
function invariantZeros2
  "Example to compute the invariant zeros of a state space system"
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;
  import Complex;

  output Complex iz[:];
  output Boolean ok;

protected
  Real A[5,5] = [17.0,   24.0,    1.0,    8.0,   15.0;
                 23.0,    5.0,    7.0,   14.0,   16.0;
                  4.0,    6.0,   13.0,   20.0,   22.0;
                 10.0,   12.0,   19.0,   21.0,    3.0;
                 11.0,   18.0,   25.0,    2.0,    9.0];

  Real B[5,2] = [-1.0,   -4.0;
                  4.0,    9.0;
                 -9.0,  -16.0;
                 16.0,   25.0;
                -25.0,  -36.0];

  Real C[2,5] = [1, 0, 1, 0, 0;
                 0, 1, 0, 1, 1];

  Real D[2,2] = [0, 0;
                 0, 0];

  StateSpace ss=StateSpace(
      A=A,
      B=B,
      C=C,
      D=D);

algorithm
  iz := StateSpace.Analysis.invariantZeros(ss);
  ok := true;

end invariantZeros2;

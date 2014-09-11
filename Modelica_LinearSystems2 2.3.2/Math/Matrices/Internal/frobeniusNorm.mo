within Modelica_LinearSystems2.Math.Matrices.Internal;
function frobeniusNorm "Return the Frobenius norm of a matrix"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Input matrix";
  output Real result=if min(size(A))>0 then sqrt(sum(A.*A)) else -1e100
    "frobenius norm of matrix A";

algorithm
  annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end frobeniusNorm;

within Modelica_LinearSystems2.Internal;
function eigenValuesFromLinearization
  "Return eigen values from a set of linearization data (assumes files dslinX.mat in current directory)"
  extends Modelica.Icons.Function;

  input Real is[:] "File indices X of the dslinX.mat files";
  input Integer nx "Number of states";
  input Integer nu "Number of inputs";
  input Integer ny "Number of outputs";
  input Boolean reorder=false
    "True, if eigen values shall be reordered so that they are closest to the previous ones";
  output Real Re[nx,size(is,1)]
    "Re[nx,np] Real values of eigenvalues Re[j,i], where i are the different parameter values and j the eigenvalue numbers";
  output Real Im[nx,size(is,1)]
    "Im[nx,np] Imaginary values of eigenvalues Im[j,i], where i are the different parameter values and j the eigenvalue numbers";
protected
  Integer np = size(is,1) "Number of linearization files";
  String fileName="dslin";
  String fileName2;
  Real A[nx,nx] "A-matrix of one linearization point";
  Real ABCD[nx+ny,nx+nu] "Linearization matrices A,B,C,D";
  Integer info;
  Real Re_i[nx];
  Real Im_i[nx];
  Real Re_old[nx];
  Real Im_old[nx];
  Real Ones[nx] = ones(nx);
  Real diff[nx];
  Integer idiff[nx];
  Boolean compare[nx];
  Integer r;
  Boolean first;
  String file="log_rootLocus.txt";
algorithm
  // Read all matrices from file, compute eigenvalues and store in output arrays
  for i in 1:np loop
    // Read matrix A from file i
    fileName2 := fileName+String(i)+".mat";
    ABCD :=Modelica.Utilities.Streams.readRealMatrix(
      fileName2, "ABCD", nx + ny, nx + nu);
    A :=ABCD[1:nx, 1:nx];

    // Compute eigen values of A
    (Re_i, Im_i, info) := Modelica.Math.Matrices.LAPACK.dgeev_eigenValues(A);
     assert(info == 0, "Calculating the eigen values with function LAPACK function dgeev_eigenValues\n"+
                       "failed, since the numerical algorithm did not converge.");

    // Delete linearization file
    Modelica.Utilities.Files.removeFile(fileName2);

    // Copy eigen values to result arrays
    if not reorder or i == 1 then
      Re[:,i] :=Re_i;
      Im[:,i] :=Im_i;
    else
      compare :=fill(true, nx);
      Re_old :=Re[:, i - 1];
      Im_old :=Im[:, i - 1];

      /* Determine difference from new to old eigenvalues.
         The eigen values are inspected in the order from the smallest
         to the largest difference
      */
      diff :=(Re_i - Re_old).^2 + (Im_i - Im_old).^2;
      (,idiff) :=Modelica.Math.Vectors.sort(diff);

      for j in idiff loop
        // Find first the right order of the pure real eigen values
        // Compute squared distance of eigen value j to previous eigen values
        diff :=(Re_i[j]*Ones - Re_old).^2 + (Im_i[j]*Ones - Im_old).^2;

        // Determine index of minimum
        first :=true;
        for k in 1:nx loop
          if compare[k] then
            if first then
              first :=false;
              r :=k;
            elseif diff[k] < diff[r] then
              r :=k;
            end if;
          end if;
        end for;

        // Copy eigen value
        Re[r,i] :=Re_i[j];
        Im[r,i] :=Im_i[j];
        compare[r] :=false;
      end for;
    end if;
  end for;

  annotation (__Dymola_translate=true);
end eigenValuesFromLinearization;

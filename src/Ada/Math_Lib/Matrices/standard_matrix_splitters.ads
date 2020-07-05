with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;

package Standard_Matrix_Splitters is

-- DESCRIPTION :
--   Procedures are defined to split a matrix of complex numbers into
--   columns with the real parts and columns with the imaginary parts.
--   The mirror of the split is the merge.
--   For matrix-vector multiplication, matrices are splitted in rows.

  procedure Complex_Parts
              ( mat : in Standard_Complex_Matrices.Matrix;
                rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in rvv the real parts of the elements in mat
  --   and in ivv the imaginary parts of the elements in mat.
  --   The vector representation of the matrix mat is column wise.

  -- REQUIRED :
  --   Space is allocated for the vectors rvv and ivv, with ranges:
  --   rvv'range = ivv'range = mat'range(2), and for all k in rvv'range:
  --   rvv(k)'range = ivv(k)'range = mat'range(1).

  -- ON ENTRY :
  --   mat      a matrix with complex numbers;
  --   rvv      space allocated for the real parts of the numbers
  --            in the columns of mat;
  --   ivv      space allocated for the imaginary parts of the numbers
  --            in the columns of mat.

  -- ON RETURN :
  --   rvv      real parts of the columns of mat;
  --   ivv      imaginary parts of the columns of mat.

  procedure Complex_Merge
              ( rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec;
                mat : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given in rvv and ivv are the real and imaginary parts of the
  --   columns of a complex matrix, with values defined on return.
  --   Mirrors the Complex_Parts procedure.

  -- REQUIRED :
  --   rvv'range = ivv'range = mat'range(2), and for all k in rvv'range:
  --   rvv(k)'range = ivv(k)'range = mat'range(1).

  -- ON ENTRY :
  --   rvv      columns with the real parts of complex numbers;
  --   ivv      columns with the imaginary parts of complex numbers.

  -- ON RETURN :
  --   mat      matrix of complex numbers with real parts from rvv
  --            and imaginary parts from ivv.

  procedure Split_Rows
              ( A : in Standard_Complex_Matrices.Link_to_Matrix;
                rArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                iArows : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Splits the rows of A in vectors of real and imaginary parts of
  --   the complex numbers on the rows of A.

  -- ON ENTRY :
  --   A        a matrix of complex numbers;
  --   rArows   space allocated for the real parts of all numbers of A;
  --   iArows   space allocated for the imaginary parts of all numbers of A.

  -- ON RETURN :
  --   rArows   rArows(i)(j) contains the real part of A(i,j);
  --   iArows   iArows(i)(j) contains the imaginary part of A(i,j).

  procedure Split_Rows
              ( vm : in Standard_Complex_VecMats.VecMat;
                rv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Splits the rows of the matrix vm(k) into vectors of real and
  --   imaginary parts of the complex numbers on the rows of vm(k),
  --   for k in rv'range = iv'range.

  -- REQUIRED :
  --   rv'range fits within vm'range and 
  --   rv(k)'range = vm(k)'range(1), for all k in rv'range, and
  --   rv(k)(j)'range = vm(k)'range(2), for all j in rv(k)'range;
  --   iv has the same dimensions as rv.

  -- ON ENTRY :
  --   vm       a vector of complex matrices;
  --   rv       rv(k) has space allocated for the real parts 
  --            of all numbers of vm(k);
  --   iv       iv(k) has space allocated for the imaginary parts 
  --            of all numbers of vm(k).

  -- ON RETURN :
  --   rv       rv(k) stores the real parts of vm(k), for k in rv'range;
  --   iv       iv(k) stores the imaginary parts of vm(k), for k in iv'range.

end Standard_Matrix_Splitters;

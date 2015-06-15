with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;

package Matrix_to_VecVec_Conversions is

-- DESCRIPTION :
--   The mat2vv functions convert matrices to vectors of vectors.
--   In the vector of vectors representation, a matrix is stored
--   as a vector of columns, because many linear algebra operations
--   are oriented towards columns.

  function mat2vv ( A : Standard_Complex_Matrices.Matrix )
                  return Standard_Complex_VecVecs.VecVec;
  function mat2vv ( A : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Complex_VecVecs.VecVec;
  function mat2vv ( A : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Complex_VecVecs.VecVec;
  function mat2vv ( A : Multprec_Complex_Matrices.Matrix )
                  return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of vectors, with as content in the vectors
  --   the columns of the matrix A.

end Matrix_to_VecVec_Conversions;

with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Dense_Vector_Series;
with Standard_Dense_Matrix_Series;

package Standard_Interpolating_Series is

-- DESCRIPTION :
--   Via interpolation at random points we can solve linear systems
--   with matrix series, even if several of the leading coefficients
--   of the matrix coefficients are singular matrices.

  function Eval ( v : Standard_Dense_Vector_Series.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector obtained by evaluating v at t.

  function Eval ( m : Standard_Dense_Matrix_Series.Matrix;
                  t : Complex_Number )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix obtained by evaluating m at t.

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               d : integer32; verbose : boolean := true ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix series up to degree d <= m.deg
  --   has full rank when evaluated at a random complex value.

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               verbose : boolean := true ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 if for all degrees d in 0..m.deg the matrix has no
  --   full rank when evaluated at a random complex value, otherwise
  --   returns the smallest d for which the evaluation of m at a
  --   random complex value yields a full rank matrix.
  --   If the matrix series m does not have full rank, then the solution
  --   of a linear system with m as a coefficient matrix may not work.

  function Sample ( v : Standard_Dense_Vector_Series.Vector;
                    t : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Samples the vector series m at the points in t,
  --   returning the evaluated vectors v(t).
  --   The range of the returned result is t'range.

  function Sample ( m : Standard_Dense_Matrix_Series.Matrix;
                    t : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Samples the matrix series m at the points in t,
  --   returning the evaluated matrices m(t).
  --   The range of the returned result is t'range.

  function Solve_Linear_Systems
             ( m : Standard_Complex_VecMats.VecMat;
               v : Standard_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec;
            
  -- DESCRIPTION :
  --   Solves the linear systems m(i)*x = v(i) and returns
  --   the solutions x in the vector of vectors on return.

  -- REQUIRED : m'range = v'range and all matrices in m and
  --   vectors in v have the same dimension.

  function Residuals
             ( m : Standard_Complex_VecMats.VecMat;
               v,x : Standard_Complex_VecVecs.VecVec )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of residuals for the solutions in x
  --   of the systems defined by the matrices in m and the
  --   right hand side vectors in v.

  function Transpose ( x : Standard_Complex_VecVecs.VecVec )
                     return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   If the range of x is 0..deg and the range of each x(i)
  --   is 1..dim, then the vector of vectors on return has
  --   the range 1..dim and each vector in the vector of vectors
  --   on return has the range 1..deg+1.

  function Vandermonde_Matrix
             ( t : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Vandermonde matrix defined by the points in t.
  --   The range of the square matrix on return is 1..dim,
  --   where dim equals the number of points in t.

  function Solve_Interpolation_Systems
             ( v : Standard_Complex_Matrices.Matrix;
               f : Standard_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Solves the interpolation linear systems defined by the
  --   Vandermonde matrix in v and the right hand side vectors in f.

  function Construct ( x : Standard_Complex_VecVecs.VecVec )
                     return Standard_Dense_Vector_Series.Vector;

  -- DESCRIPTION :
  --   Takes the solutions of the interpolation systems in x
  --   and constructs the vector series.

end Standard_Interpolating_Series;

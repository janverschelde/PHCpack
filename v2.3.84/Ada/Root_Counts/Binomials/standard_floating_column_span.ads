with text_io;                           use text_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Standard_Floating_Column_Span is

-- DESCRIPTION :
--   Applies the Gram-Schmidt orthogonalization to decide whether
--   a vector belongs to the column span of a matrix.

  function In_Span ( v : Standard_Integer_VecVecs.VecVec;
                     x : Standard_Integer_Vectors.Vector )
                   return boolean;

  -- DESCRIPTION :
  --   This function is a wrapper to the floating-point versions.

  function In_Span ( v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector )
                   return boolean;
  function In_Span ( v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector;
                     tol : double_float ) return boolean;
  function In_Span ( file : file_type;
                     v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector;
                     tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x can be written as a linear combination
  --   of the vectors in v.

  -- REQUIRED : the vectors in v are linearly independent.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   v        a sequence of vectors;
  --   x        one vector of the same range as the vectors in v;
  --   tol      tolerance to decide on the residual of the
  --            least squares problem (the default is 1.0E-8).

  -- ON RETURN :
  --   true if x can be written as a linear combination of the
  --   vectors in v (with respect to the given tolerance),
  --   and false otherwise.

end Standard_Floating_Column_Span;

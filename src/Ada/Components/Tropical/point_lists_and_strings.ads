with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;

package Point_Lists_and_Strings is

-- DESCRIPTION :
--   This package provides functions to convert a point configuration given
--   by a matrix with in its columns the integer coordinates of points
--   into a list of tuples for export to Python.
--   Also routines to parse a Python list of tuples into a matrix of points
--   are provided.

  function convert ( A : Standard_Integer64_Matrices.Matrix ) return string;
  function convert ( A : Multprec_Integer_Matrices.Matrix ) return string;

  -- DESCRIPTION :
  --   Returns a list of tuples, suitable for parsing by Python.
  --   The number of items in the string representation of the Python list
  --   equals the number of columns of the matrix A.
  --   Every column defines one tuple in the list on return.

end Point_Lists_and_Strings;

with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with Lists_of_Integer_Vectors;

package Point_Lists_and_Strings is

-- DESCRIPTION :
--   This package provides functions to convert a point configuration given
--   by a matrix with in its columns the integer coordinates of points
--   into a list of tuples for export to Python.  As an example of a matrix
--   format of a configuration of 5 points in 3-space, consider:
--
--    -9 -3 8 -7 -4
--    -1 -8 9 -7 0
--    3 -4 -1 -3 9
--
--   The corresponding string representation is
--   [(-9, -1, 3), (-3, -8, -4), (8, 9, -1), (-7, -7, -3), (-4, 0, 9)].
--
--   Also routines to parse a Python list of tuples into a matrix of points
--   are provided.

  function write ( A : Standard_Integer64_Matrices.Matrix ) return string;
  function write ( A : Multprec_Integer_Matrices.Matrix ) return string;

  -- DESCRIPTION :
  --   Returns a list of tuples, suitable for parsing by Python.
  --   The number of items in the string representation of the Python list
  --   equals the number of columns of the matrix A.
  --   Every column defines one tuple in the list on return.

  function write ( A : Lists_of_Integer_Vectors.List ) return string;

  -- DESCRIPTION :
  --   Given in A the list of points, writes the string representation
  --   as define when A is the matrix with the points in its columns.

  procedure Extract_Dimensions ( s : in string; rows,cols : out integer32 );

  -- DESCRIPTION :
  --   Extracts the dimensions, number of rows and columsn from the
  --   string representation of the point configuration.
  --   The number of columns equals the number of opening round brackets.
  --   When parsing the first point, the number of rows equals
  --   the number of commas plus one.

  function parse ( s : string; rows,cols : integer32 )
                 return Standard_Integer64_Matrices.Matrix;
  function parse ( s : string; rows,cols : integer32 )
                 return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Parses the string and returns a matrix of points,
  --   of dimensions 1..rows, 1..cols.

  -- REQUIRED :
  --   rows and cols are computed via Extract_Dimensions(s,rows,cols).

  function parse ( n : integer32; s : string )
                 return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Extracts the numbers of the string representation s of the point
  --   with n integer coordinates as stored as a Python tuple.

end Point_Lists_and_Strings;

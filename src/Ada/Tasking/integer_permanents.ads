with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Boolean_Matrices;

package Integer_Permanents is

-- DESCRIPTION :
--   Reference code to develop a multitasked permanent computation.
--   For Boolean matrices, the permanent computation is a pure
--   counting problem of contributing permutations.

  procedure Permanent
              ( row : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 );

  -- DESCRIPTION :
  --   Row expansion for the permanent of an integer matrix.

  -- ON ENTRY :
  --   row      current row index, initialize with mat'first(1);
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   per      initialize with zero.

  -- ON RETURN :
  --   per      the permanent of the matrix mat.

  procedure Permanent
              ( file : in file_type; row : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 );

  -- DESCRIPTION :
  --   Row expansion for the permanent of an integer matrix.
  --   For every factor which contributes to the permanent,
  --   writes the column selection and the value of the factor.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   row      current row index, initialize with mat'first(1);
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   per      initialize with zero.

  -- ON RETURN :
  --   per      the permanent of the matrix mat.

  procedure Permanent
              ( row : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 );

  -- DESCRIPTION :
  --   Row expansion for the permanent of a Boolean matrix.

  -- ON ENTRY :
  --   row      current row index, initialize with mat'first(1);
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   per      initialize with zero.

  -- ON RETURN :
  --   per      the permanent of the matrix mat.

  procedure Permanent
              ( file : in file_type; row : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 );

  -- DESCRIPTION :
  --   Row expansion for the permanent of a Boolean matrix.
  --   For every factor which contributes to the permanent,
  --   writes the column selection.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   row      current row index, initialize with mat'first(1);
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   per      initialize with zero.

  -- ON RETURN :
  --   per      the permanent of the matrix mat.

  function Permanent ( mat : Standard_Integer_Matrices.Matrix )
                     return integer64;

  -- DESCRIPTION :
  --   Returns the permanent of the matrix mat.

  function Permanent ( file : file_type;
                       mat : Standard_Integer_Matrices.Matrix )
                     return integer64;

  -- DESCRIPTION :
  --   Returns the permanent of the matrix mat.
  --   Writes every contributing factor to file.

  function Permanent ( mat : Boolean_Matrices.Matrix )
                     return integer64;

  -- DESCRIPTION :
  --   Returns the permanent of the matrix mat.

  function Permanent ( file : file_type;
                       mat : Boolean_Matrices.Matrix )
                     return integer64;

  -- DESCRIPTION :
  --   Returns the permanent of the matrix mat.
  --   Writes every contributing factor to file.

end Integer_Permanents;

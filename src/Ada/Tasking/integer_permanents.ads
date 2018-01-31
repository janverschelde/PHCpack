with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
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

  function Number_of_Start_Columns
             ( dim,nbrows : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the maximum number of start columns for the first number
  --   of rows in nbrows of a matrix of dimension dim.

  procedure Start_Columns
              ( row,nbrows : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                idx : in out integer32;
                stc : in out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the column indices for the first rows in the
  --   row expansion for the permanent of an integer matrix.

  -- ON ENTRY :
  --   row      current row index, initialize with mat'first(1);
  --   nbrows   number of rows to consider;
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   idx      initialize index to zero;
  --   stc      must have the range 
  --            1..Number_of_Start_Columns(mat'last(1),nbrows).

  -- ON RETURN :
  --   idx      last index in stc of start column indices;
  --   stc      column indices for the first number of rows nbrows 
  --            in the row expansion for the permanent of mat.

  procedure Start_Columns
              ( row,nbrows : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                idx : in out integer32;
                stc : in out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the column indices for the first rows in the
  --   row expansion for the permanent of a Boolean matrix.

  -- ON ENTRY :
  --   row      current row index, initialize with mat'first(1);
  --   nbrows   number of rows to consider;
  --   mat      a square matrix;
  --   cols     selected column indices of the matrix,
  --            the range of cols must be mat'range(2);
  --   cnts     counts the number of times to select each column,
  --            of range mat'range(2) and initialized to ones;
  --   idx      initialize index to zero;
  --   stc      must have the range 
  --            1..Number_of_Start_Columns(mat'last(1),nbrows).

  -- ON RETURN :
  --   idx      last index in stc of start column indices;
  --   stc      column indices for the first number of rows nbrows 
  --            in the row expansion for the permanent of mat.

  function Start_Permanent
             ( row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Standard_Integer_Matrices.Matrix )
             return integer64;

  -- DESCRIPTION :
  --   Given the start columns for the first rows in stc,
  --   returns the permanent of the integer matrix mat.
  --   The value of row is the next row after the first rows.

  function Start_Permanent
             ( file : file_type; row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Standard_Integer_Matrices.Matrix )
             return integer64;

  -- DESCRIPTION :
  --   Given the start columns for the first rows in stc,
  --   returns the permanent of the integer matrix mat.
  --   The value of row is the next row after the first rows.
  --   The contributing permutations are written to file.

  function Start_Permanent
             ( row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Boolean_Matrices.Matrix )
             return integer64;

  -- DESCRIPTION :
  --   Given the start columns for the first rows in stc,
  --   returns the permanent of the Boolean matrix mat.
  --   The value of row is the next row after the first rows.

  function Start_Permanent
             ( file : file_type; row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Boolean_Matrices.Matrix )
             return integer64;

  -- DESCRIPTION :
  --   Given the start columns for the first rows in stc,
  --   returns the permanent of the Boolean matrix mat.
  --   The value of row is the next row after the first rows.
  --   The contributing permutations are written to file.

end Integer_Permanents;

with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Boolean_Matrices;
with Standard_Integer_Matrices;

package Multitasking_Integer_Permanents is

-- DESCRIPTION :
--   Exports procedures to compute the permanent of an integer matrix
--   with many tasks.

  procedure Initialize
              ( mat : in Boolean_Matrices.Matrix;
                nbrows : in integer32; size : out integer32 );

  -- DESCRIPTION :
  --   Initializes the static columns queue with jobs,
  --   for the multitasked row expansion of a Boolean matrix.

  -- ON ENTRY :
  --   mat      a square Boolean matrix;
  --   nbrows   the number of first rows in mat.

  -- ON RETURN :
  --   size     the number of columns stored in the queue.

  procedure Initialize
              ( mat : in Standard_Integer_Matrices.Matrix;
                nbrows : in integer32; size : out integer32 );

  -- DESCRIPTION :
  --   Initializes the static columns queue with jobs,
  --   for the multitasked row expansion of an integer matrix.

  -- ON ENTRY :
  --   mat      a square integer matrix;
  --   nbrows   the number of first rows in mat.

  -- ON RETURN :
  --   size     the number of columns stored in the queue.

  function Permanent ( mat : Boolean_Matrices.Matrix;
                       nbrows,ntasks : integer32; verbose : boolean )
                     return integer64;

  -- DESCRIPTION :
  --   Applies row expansion to compute the permanent of a square
  --   Boolean matrix in mat, with multitasking.

  -- ON ENTRY :
  --   mat      a square Boolean matrix;
  --   nbrows   number of rows to define the jobs;
  --   ntasks   number of tasks;
  --   verbose  if true, then the tasks write output to screen,
  --            otherwise, the workers remain silent.

  -- ON RETURN :
  --   The permanent of the matrix mat.

  function Permanent ( mat : Standard_Integer_Matrices.Matrix;
                       nbrows,ntasks : integer32; verbose : boolean )
                     return integer64;

  -- DESCRIPTION :
  --   Applies row expansion to compute the permanent of a square
  --   integer matrix in mat, with multitasking.

  -- ON ENTRY :
  --   mat      a square integer matrix;
  --   nbrows   number of rows to define the jobs;
  --   ntasks   number of tasks;
  --   verbose  if true, then the tasks write output to screen,
  --            otherwise, the workers remain silent.

  -- ON RETURN :
  --   The permanent of the matrix mat.

end Multitasking_Integer_Permanents;

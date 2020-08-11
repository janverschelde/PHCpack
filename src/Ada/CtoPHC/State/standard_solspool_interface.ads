with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Standard_SolsPool_Interface is

-- DESCRIPTION :
--   The functions below work on a pool of solution lists in double precision.

  function Standard_SolsPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the pool of solution lists in double precision
  --   to hold a given number of lists.

  -- ON ENTRY :
  --   a       in a[0] is the number of lists to be stored in the pool;
  --   vrblvl  is the verbose level.

  function Standard_SolsPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of solution lists currently stored
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the size of the solutions pool.

  function Standard_SolsPool_Length
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the length of a solution list in the pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the list in the pool;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the length of the solution list,
  --           with the index given in a[0].

  function Standard_SolsPool_Dimension
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of a solution list in the pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the list in the pool;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the dimension of the solution list,
  --           with the index given in a[0].

  function Standard_SolsPool_Add
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a solution to a list in the pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the list in the pool;
  --   b       integer attributes of the solution;
  --   c       double attributes of the solution;
  --   vrblvl  is the verbose level.

  function Standard_SolsPool_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a solution from a list in the pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the list in the pool;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       integer attributes of the solution;
  --   c       double attributes of the solution.

end Standard_SolsPool_Interface;

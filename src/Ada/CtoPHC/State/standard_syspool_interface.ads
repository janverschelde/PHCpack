with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Standard_SysPool_Interface is

-- DESCRIPTION :
--   The functions below interface to the pool of polynomial systems
--   with coefficients in double precision.

  function Standard_SysPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the pool of systems in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of systems the pool can hold;
  --   vrblvl  is the verbose level.

  function Standard_SysPool_Read
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system in double precision
  --   and places it at the given position in the pool of systems
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the system in the pool;
  --   vrblvl  is the verbose level.

  function Standard_SysPool_Write
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the system at the given position in the pool of
  --   systems in double precision to the defined output.

  -- ON ENTRY :
  --   a       in a[0] is the index of the system in the pool;
  --   vrblvl  is the verbose level.

  function Standard_SysPool_from_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a copy of a system in double precision
  --   and places it at the given position in the pool of systems
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index in the pool of systems
  --           where to place the copy of the system in the container;
  --   vrblvl  is the verbose level.

  function Standard_SysPool_into_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   The polynomials in the system at the given position in the pool
  --   of systems in double precision are set into the container
  --   to store a polynomial system in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index in the pool of systems
  --           where to get the system into the container;
  --   vrblvl  is the verbose level.

  function Standard_SysPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the pool of systems in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the size of the systems pool.

  function Standard_SysPool_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the pool of systems in double precision.
  --   The verbose level is given in vrblvl.

end Standard_SysPool_Interface;

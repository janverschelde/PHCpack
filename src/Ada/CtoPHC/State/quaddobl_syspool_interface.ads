with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package QuadDobl_SysPool_Interface is

-- DESCRIPTION :
--   The functions below interface to the pool of polynomial systems
--   with coefficients in quad double precision.

  function QuadDobl_SysPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the pool of systems in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of systems the pool can hold;
  --   vrblvl  is the verbose level.

  function QuadDobl_SysPool_from_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a copy of a system in quad double precision
  --   and places it at the given position in the pool of systems
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index in the pool of systems
  --           where to place the copy of the system in the container;
  --   vrblvl  is the verbose level.

  function QuadDobl_SysPool_into_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   The polynomials in the system at the given position in the pool
  --   of systems in quad double precision are set into the container
  --   to store a polynomial system in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index in the pool of systems
  --           where to get the system into the container;
  --   vrblvl  is the verbose level.

  function QuadDobl_SysPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the pool of systems in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the size of the systems pool.

  function QuadDobl_SysPool_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the pool of systems in quad double precision.
  --   The verbose level is given in vrblvl.

end QuadDobl_SysPool_Interface;

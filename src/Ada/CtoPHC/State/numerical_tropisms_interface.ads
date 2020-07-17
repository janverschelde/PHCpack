with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Numerical_Tropisms_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   numerically computed tropisms.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function Standard_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the container for numerical tropisms computed
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] are the number of tropisms, and
  --           in a[1] is the length of each tropism;
  --   b       the winding numbers, as many as the value of a[0]
  --   c       all floating-point coefficients,
  --           as many as a[0]*(a[1]+1), for the coordinates of the
  --           tropisms and the errors on the numerical tropisms,
  --           or as many as 2*a[0]*(a[1]+1) in double double precision,
  --           or as many as 4*a[0]*(a[1]+1) in quad double precision;
  --   vrblvl  the verbose level.

  function Store_Standard_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function Store_DoblDobl_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function Store_QuadDobl_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores in the container a tropism computed in 
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the tropism, and
  --           in a[1] the index of the tropism;
  --   b       the winding number;
  --   c       all the floating-point coefficients, with the coordinates
  --           of the tropisms and its error,
  --           as many as a[0] + 1 in double precision, or
  --           as many as 2*a[0] + 2 in double double precision, or
  --           as many as 4*a[0] + 4 in quad double precision;
  --   vrblvl  the verbose level.

  function Standard_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of tropisms computed 
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of tropisms in the container.

  function Standard_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the tropisms computed
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the dimension of the tropisms in the container.

  function Standard_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves from the container one tropism computed 
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the tropism, and
  --           in a[1] is the index of the tropism;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the winding number of the tropism;
  --   c       the coordinates of the tropism and its error,
  --           as many as a[0] + 1 in double precision, or
  --           as many as 2*a[0] + 2 in double double precision, or
  --           as many as 4*a[0] + 4 in quad double precision.

  function Standard_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves from the container all tropisms computed
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       in a[0] are the number of tropisms, and
  --           in a[1] the length of each tropism;
  --   b       the winding numbers for each tropism;
  --   c       the coefficients, coordinates and errors, 
  --           as many as a[0]*(a[1]+1) in double precision, or
  --           as many as 2*a[0]*(a[1]+1) in double double precision, or
  --           as many as 4*a[0]*(a[1]+1) in quad double precision.

  function Standard_Clear ( vrblvl : integer32 := 0 ) return integer32;
  function DoblDobl_Clear ( vrblvl : integer32 := 0 ) return integer32;
  function QuadDobl_Clear ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates data for numerical tropisms computed
  --   in standard double, double double, or quad double precision.
  --   The verbose level is given in vrblvl.

end Numerical_Tropisms_Interface;

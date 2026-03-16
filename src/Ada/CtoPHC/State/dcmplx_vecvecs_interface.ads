with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package DCMPLX_VecVecs_Interface is

-- DESCRIPTION :
--   Defines functions to access the container of arrays of vectors
--   of complex vectors.

  function DCMPLX_VecVecs_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the container with the number of arrays
  --   and the number of vectors in each array.

  -- ON ENTRY :
  --   a       in a[0] is the number of arrays;
  --   b       array with a[0] integers,
  --           b[i] contains the number of vectors in the (i+1)-th array;
  --   vrblvl  is the verbose level.

  function DCMPLX_VecVecs_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns in a[0] the number of arrays in the container.
  --   The verbose level is given by vrblvl.

  function DCMPLX_VecVecs_Get_Size_Array
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns in b[0] the number of vectors in the k-th array,
  --   where k is given by a[0].
  --   The verbose level is given by vrblvl.

  function DCMPLX_VecVecs_Get_Size_Vector
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns in b[0] the size of the i-th vector in the k-th array,
  --   where k is given by a[0] and i by a[1].
  --   The verbose level is given by vrblvl.

  function DCMPLX_VecVecs_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the vector with coordinates given in c
  --   to the proper place in the container.

  -- ON ENTRY :
  --   a       dimension of the vector is in a[0];
  --   b       defines the i-vector of the k-th array
  --           as k = b[0] and i = b[1];
  --   c       as many doubles as the value of 2*a[0];
  --   vrblvl  is the verbose level.

  function DCMPLX_VecVecs_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns in c the coordinates of the vector
  --   as defined by the input values of a and b.

  -- ON ENTRY :
  --   a       dimension of the vector is in a[0];
  --   b       defines the i-vector of the k-th array
  --           as k = b[0] and i = b[1];
  --   c       has space for as many doubles as the value of 2*a[0];
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the i-th vector in the k-th array.

  function DCMPLX_VecVecs_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates space occupied by the vectors in the container.

end DCMPLX_VecVecs_Interface;

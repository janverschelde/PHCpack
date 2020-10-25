with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;

package TripDobl_Vector_Splitters is

-- DESCRIPTION :
--   Provides memory allocators for triple double complex vectors.

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return TripDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the complex coefficients
  --   of a series truncated to degree deg.

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return TripDobl_Complex_VecVecs.VecVec;
  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return TripDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of a vector 
  --   of series, all truncated to degree deg.
  --   The vector on return has range 1..dim.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return TripDobl_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return TripDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim.

end TripDobl_Vector_Splitters;

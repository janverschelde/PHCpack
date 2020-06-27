with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Standard_Vector_Splitters is

-- DESCRIPTION :
--   A vector splitter separates real and imaginary parts of complex vectors,
--   for more efficient computations with the vectors.

-- FUNCTIONS WITH MEMORY ALLOCATIONS :
--   The functions below are easy to use for correctness checks,
--   as they automatically allocate memory, good for one time use.
--   But memory allocations are not appropriate for work space
--   and separating real and imaginary parts is better done in
--   one single loop, instead of two separate loops.

  function Real_Part ( x : Standard_Complex_Vectors.Vector )
                     return Standard_Floating_Vectors.Vector;
  function Real_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector;
  function Real_Part ( x : Standard_Complex_VecVecs.Link_to_VecVec )
                     return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector (or vector of vectors) of the real parts
  --   of the complex vector x.  Memory is allocated.

  function Imag_Part ( x : Standard_Complex_Vectors.Vector )
                     return Standard_Floating_Vectors.Vector;
  function Imag_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector;
  function Imag_Part ( x : Standard_Complex_VecVecs.Link_to_VecVec )
                     return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector (or vector of vectors) of the imaginary parts
  --   of the complex vector x.  Memory is allocated.

  function Make_Complex
             ( rpx,ipx : Standard_Floating_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Make_Complex
             ( rpx,ipx : Standard_Floating_Vectors.Link_to_Vector )
             return Standard_Complex_Vectors.Link_to_Vector;
  function Make_Complex
             ( rpx,ipx : Standard_Floating_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec;
  function Make_Complex
             ( rpx,ipx : Standard_Floating_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the vector (of vectors) of complex numbers,
  --   with real and imaginary parts given in the vectors rpx and ipx.
  --   Memory is allocated for the resulting vectors.

  procedure Split_Complex
              ( x : in Standard_Complex_Vectors.Vector;
                rpx,ipx : out Standard_Floating_Vectors.Vector );
  procedure Split_Complex
              ( x : in Standard_Complex_Vectors.Vector;
                rpx,ipx : out Standard_Floating_Vectors.Link_to_Vector );
  procedure Split_Complex
              ( x : in Standard_Complex_Vectors.Link_to_Vector;
                rpx,ipx : out Standard_Floating_Vectors.Link_to_Vector );
  procedure Split_Complex
              ( x : in Standard_Complex_VecVecs.VecVec;
                rpx,ipx : out Standard_Floating_VecVecs.VecVec );
  procedure Split_Complex
              ( x : in Standard_Complex_VecVecs.Link_to_VecVec;
                rpx,ipx : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Splits the complex vector (of vectors) x into two real vectors,
  --   with its real and imaginary parts of the complex numbers in x.
  --   Memory is allocated for the resulting vectors.

-- MEMORY ALLOCATORS :

  function Allocate_Floating_Coefficients
             ( deg : integer32 )
             return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the real or imaginary parts of
  --   complex coefficients of a series truncated to degree deg.

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return Standard_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the complex coefficients
  --   of a series truncated to degree deg.

  function Allocate_Floating_Coefficients
             ( dim,deg : integer32 )
             return Standard_Floating_VecVecs.VecVec;
  function Allocate_Floating_Coefficients
             ( dim,deg : integer32 )
             return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the real or imaginary parts of
  --   coefficients of a vector of series, all truncated to degree deg.
  --   The vector on return has range 1..dim.

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec;
  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return Standard_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of a vector 
  --   of series, all truncated to degree deg.
  --   The vector on return has range 1..dim.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Floating_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim.

-- PROCEDURES TO PART AND MERGE :

  procedure Complex_Parts
              ( x : in Standard_Complex_Vectors.Link_to_Vector;
                rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( x : in Standard_Complex_Vectors.Vector;
                rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( x : in Standard_Complex_VecVecs.VecVec;
                rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( x : in Standard_Complex_VecVecs.Link_to_VecVec;
                rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two real vectors,
  --   with its real and imaginary parts of the complex numbers in x.

  -- REQUIRED :
  --   The vectors rpx and ipx are completely allocated
  --   and have the same ranges as x.

  procedure Complex_Parts
              ( deg : in integer32;
                x : in Standard_Complex_Vectors.Vector;
                rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in Standard_Complex_VecVecs.VecVec;
                rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in Standard_Complex_VecVecs.Link_to_VecVec;
                rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two real vectors,
  --   with its real and imaginary parts of the complex numbers in x,
  --   with the vectors x(k) splitted up to degree deg.

  -- REQUIRED :
  --   The vectors rpx and ipx are allocated, at least up to deg,
  --   and have the same ranges as x, at least up to the index deg.

  procedure Complex_Merge
             ( rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector;
               cvx : in Standard_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
             ( rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector;
               cvx : out Standard_Complex_Vectors.Vector );
  procedure Complex_Merge
             ( rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec;
               cvx : in Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
             ( rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec;
               cvx : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts in rpx and ipx
  --   into the vector (of vectors) of complex numbers.

  -- REQUIRED :
  --   The vector cvx is completely allocated
  --   and has the same ranges as rpx and ipx.

  procedure Complex_Merge
             ( deg : in integer32;
               rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector;
               cvx : in Standard_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
             ( deg : in integer32;
               rpx,ipx : in Standard_Floating_Vectors.Link_to_Vector;
               cvx : out Standard_Complex_Vectors.Vector );
  procedure Complex_Merge
             ( deg : in integer32;
               rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec;
               cvx : in Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
             ( deg : in integer32;
               rpx,ipx : in Standard_Floating_VecVecs.Link_to_VecVec;
               cvx : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts in rpx and ipx
  --   into the vector (of vectors) of complex numbers,
  --   with the vectors cvx(k) merged to degree deg.

  -- REQUIRED :
  --   The vector cvx is allocated, at least up to the index deg,
  --   and has the same ranges as rpx and ipx, at least up to deg.

end Standard_Vector_Splitters;

with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_VecVecs;

package Standard_Complex_VecVecVecs is

-- DESCRIPTION :
--   A VecVecVec is a three dimensional data structure of vectors 
--   of vectors of vectors of complex numbers in double precision.

  type VecVecVec is
    array ( integer32 range <> ) of Standard_Complex_VecVecs.Link_to_VecVec;
 
  type Link_to_VecVecVec is access VecVecVec;

  type VecVecVec_Array is array ( integer32 range <> ) of Link_to_VecVecVec;

  procedure Allocate ( v : out Link_to_VecVecVec;
                       d1first,d1last : in integer32;
                       d2first,d2last : in integer32;
                       d3first,d3last : in integer32 );

  -- DESCRIPTION :
  --   Allocates a three dimensional matrix of complex numbers,
  --   with leading dimension in the range d1first..d1last,
  --   the second dimension has the range d2first..d2last, and
  --   the third dimension has the range d3first..d3last.
  --   All numbers in v are initialized to zero.

  procedure Copy ( v_from : in Link_to_VecVecVec;
                   v_to : out Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Copies the v_from to the v_to,
  --   after the deallocation of v_to.

  procedure Clear ( v : in out VecVecVec );
  procedure Clear ( v : in out Link_to_VecVecVec );
  procedure Clear ( v : in out VecVecVec_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by v.

end Standard_Complex_VecVecVecs;

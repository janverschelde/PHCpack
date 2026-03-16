with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Double_VecVecs_Container is

-- DESCRIPTION :
--   Defines operations to store an array of vectors
--   of floating-point vectors in double precision.
--   The verbose level of all procedures and functions
--   is defined by the value of vrblvl.

  procedure Initialize ( n : in integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Allocates space for n arrays of floating-point vectors of vectors.

  procedure Initialize ( m : in Standard_Integer_Vectors.Vector;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Allocates space for m(k) vectors of vectors as the k-th component
  --   in an array of vectors of vectors, of range equal to m'range.

  procedure Store_Link
              ( v : in Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Stores the given pointer v to the arrays.

  procedure Store_Copy ( v : in Standard_Floating_VecVecs.Array_of_VecVecs;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Makes a deep copy of v and stores the arrays.

  procedure Store_Link ( k : in integer32;
                         v : in Standard_Floating_VecVecs.Link_to_VecVec;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Stores the given pointer v as the k-th component of the arrays,
  --   where k must be in range 1..n.

  procedure Store_Copy ( k : in integer32;
                         v : in Standard_Floating_VecVecs.VecVec;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Stores a deep copy of v as the k-th component of the arrays,
  --   where k must be in range 1..n.

  procedure Store_Link ( k,i : in integer32;
                         v : in Standard_Floating_Vectors.Link_to_Vector;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Stores the given pointer v as the i-th component
  --   of the k-th vector of vectors,
  --   where k must be in range 1..n and i must be in the range
  --   of the k-th vector of vectors.

  procedure Store_Copy ( k,i : in integer32;
                         v : in Standard_Floating_Vectors.Vector;
                         vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Stores a deep copy of v as the i-th component
  --   of the k-th vector of vectors,
  --   where k must be in range 1..n and i must be in the range
  --   of the k-th vector of vectors.

  function Size ( vrblvl : in integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of vectors of vectors stored.

  function Size ( k : integer32; vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the k-th vector of vectors stored.

  -- REQUIRED : k in range 1..n, where n equals Size.

  function Size ( k,i : integer32;
                  vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the i-th vector of the k-th vector of vectors.

  -- REQUIRED : k in range 1..n, where n equals Size,
  --   and i is in the range of Size(k).

  function Get ( vrblvl : in integer32 := 0 )
               return Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns the pointer to the arrays of vectors stored.

  function Get ( k : integer32; vrblvl : integer32 := 0 )
               return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns a pointer to the k-th vector of vectors stored.

  -- REQUIRED : k in range 1..n, where n equals Size.

  function Get ( k,i : integer32; vrblvl : integer32 := 0 )
               return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns a pointer to the i-th vector
  --   in the k-th vector of vectors stored.

  -- REQUIRED : k in range 1..n, where n equals Size,
  --   and i in range 1..Size(k).

  procedure Clear ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Deallocates everything in the container.

end Double_VecVecs_Container;

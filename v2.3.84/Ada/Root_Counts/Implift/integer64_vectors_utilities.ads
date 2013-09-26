with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;         use Standard_Integer64_Vectors;
with Standard_Integer64_Transformations; use Standard_Integer64_Transformations;

package Integer64_Vectors_Utilities is

-- DESCRIPTION :
--   This package offers utilities for transforming integer vectors.

  function Pivot ( v : Vector ) return integer32;
  function Pivot ( v : Link_to_Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the first nonzero entry out of v.
  --   If all entries of v are zero, then v'last+1 is returned.

  function  Reduce ( v : Vector; i : integer32 ) return Vector;
  procedure Reduce ( v : in out Link_to_Vector; i : in integer32 );
  function  Reduce ( v : Link_to_Vector; i : integer32 ) return Link_to_Vector;

  -- DESCRIPTION :
  --   The i-th component will be deleted out of the vector.

  function  Insert ( v : Vector; i : integer32; a : integer64 ) return Vector;
  procedure Insert ( v : in out Link_to_Vector;
                     i : in integer32; a : in integer64 );
  function  Insert ( v : Link_to_Vector;
                     i : integer32; a : integer64 ) return Link_to_Vector;

  -- DESCRIPTION :
  --   The i-th component will be inserted, using the value a.

  function Insert_and_Transform
             ( v : Vector; i : integer32; a : integer64; t : Transfo )
             return Vector;
  procedure Insert_and_Transform
             ( v : in out Link_to_Vector; i : in integer32; a : integer64;
               t : in Transfo );
  function Insert_and_Transform
             ( v : Link_to_Vector; i : integer32; a : integer64; t : Transfo )
             return Link_to_Vector;

  -- DESCRIPTION :
  --   Inserts the i-th component in the vector v,
  --   using the value a, and transforms the vector, applying t.

end Integer64_Vectors_Utilities;

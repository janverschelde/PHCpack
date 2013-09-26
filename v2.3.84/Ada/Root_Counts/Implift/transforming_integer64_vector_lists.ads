with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;         use Standard_Integer64_Vectors;
with Lists_of_Integer64_Vectors;         use Lists_of_Integer64_Vectors;
with Standard_Integer64_Transformations; use Standard_Integer64_Transformations;

package Transforming_Integer64_Vector_Lists is

-- DESCRIPTION :
--   Unimodular transformations of lists of standard 64-bit integer vectors.

  procedure Shift ( L : in out List; v : in Vector );
  procedure Shift ( L : in out List; v : in Link_to_Vector );

  function  Shift ( L : List; v : Vector ) return List;
  function  Shift ( L : List; v : Link_to_Vector ) return List;

  -- DESCRIPTION :
  --   The list will be shifted: Shift(l,v) = { y-v | Is_In(l,y) }

  function "*"( L : List; t : Transfo ) return List;
  function "*"( t : Transfo; L : List ) return List;

  -- DESCRIPTION :
  --   Returns the transformed list of points.

  procedure Apply ( L : in out List; t : in Transfo );

  -- DESCRIPTION :
  --   Applies the transformation t to the list L.

  function  Reduce ( L : List; i : integer32 ) return List;
  procedure Reduce ( L : in out List; i : in integer32 );

  -- DESCRIPTION :
  --   Returns a list of vectors where the i-th component has been deleted.

  function  Insert ( L : List; i : integer32; a : integer64 ) return List;
  procedure Insert ( L : in out List; i : in integer32; a : integer64 );

  -- DESCRIPTION :
  --   Returns a list of vectors where the i-th component has been inserted,
  --   for all d in L: d(i) = a.

  function  Transform_and_Reduce ( t : Transfo; i : integer32; L : List )
                                 return List;
  procedure Transform_and_Reduce ( t : in Transfo; i : in integer32;
                                   L : in out List );

  -- DESCRIPTION :
  --   Transforms the list L and deletes the i-th component
  --   of every element in the transformed list.

  function Insert_and_Transform
             ( L : List; i : integer32; a : integer64; t : Transfo )
             return List;
  procedure Insert_and_Transform
             ( L : in out List; i : in integer32; a : in integer64;
               t : in Transfo );

  -- DESCRIPTION :
  --   Inserts the i-th component of every element in the list L,
  --   using the value a, and transforms the list, applying t.

end Transforming_Integer64_Vector_Lists;
